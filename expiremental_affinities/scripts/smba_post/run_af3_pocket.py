"""Orchestrator for the AF3-pocket redocking experiment.

For each (protein, ligand) row in docking_results.csv with Ki or Kd affinity:
  1. AF3 protein+ligand inference         (alphafold/docking_inference.py)
  2. AF3 ligand prep (per-conf pdbqts)    (alphafold/prepare_ligand.py)
  3. Per-conf vina redocking              (defined below)

Stage 2 produces one Meeko-prepped pdbqt per SMBA conformation, pre-centered at
that conformation's AF3-aligned ligand COM. Stage 3 just reads each pdbqt and
docks against the matching receptor at the pdbqt's centroid.
"""

import csv
import shutil
import signal
import tarfile
import time
from argparse import ArgumentParser, Namespace
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from tqdm import tqdm

from alphafold.docking_inference import run_one as af3_run_one
from alphafold.prepare_ligand import run_one as ligand_prep_run_one
from utils import DATA_DIR, load_pairs
from varidock.runners.af3 import AF3Config
from varidock.stages import AF3Inference
from varidock.stages.vina_dock import VinaDocking, VinaDockingConfig
from varidock.types import DockingInput, PDBQT, PocketCenter

DOCKING_RESULTS_CSV = DATA_DIR.parent / "docking_results.csv"
KI_KD_TYPES = {"ki_nM", "kd_nM"}
NUM_CONFS = 11
MAX_DOCKING_WORKERS = 16
MAX_PROTEIN_WORKERS = 4
TAR = True

def load_args() -> Namespace:
    parser = ArgumentParser(
        description="Run AF3-pocket redocking pipeline for Ki/Kd pairs."
    )
    parser.add_argument("--write-only", action="store_true")
    parser.add_argument("--overwrite-inp", action="store_true")
    return parser.parse_args()


def filter_pairs() -> list[tuple[str, str, str]]:
    """Return (protein_id, ligand_id, smiles) for docking_results.csv rows with Ki/Kd."""
    smiles_by_pair = {(pid, lid): smi for pid, lid, smi in load_pairs()}

    pairs = []
    with open(DOCKING_RESULTS_CSV, newline="") as f:
        for row in csv.DictReader(f):
            if row["known_affinity_type"] not in KI_KD_TYPES:
                continue
            pid, lid = row["protein"], row["ligand"]
            smi = smiles_by_pair.get((pid, lid))
            if smi is None:
                print(f"Skipping {pid}+{lid}: no SMILES in known_structures.csv")
                continue
            pairs.append((pid, lid, smi))
    return pairs


def pdbqt_centroid(pdbqt_path) -> tuple[float, float, float]:
    coords = []
    for line in pdbqt_path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            coords.append(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
    c = np.array(coords).mean(axis=0)
    return float(c[0]), float(c[1]), float(c[2])



def is_complete(log_file: Path) -> bool:
    if not log_file.exists() or log_file.stat().st_size == 0:
        return False
    with open(log_file) as f:
        return f.read().rstrip().endswith("COMPLETED")



def gather_jobs(protein_id: str, ligand_id: str) -> list[tuple]:
    protein_dir = DATA_DIR / protein_id
    ligand_dir = protein_dir / "ligands_prepared_af3_pocket" / ligand_id
    receptor_dir = protein_dir / "protein_receptors"
    docking_dir = protein_dir / "docking_af3_pocket" / ligand_id

    jobs = []
    for i in range(NUM_CONFS):
        receptor = receptor_dir / f"protein_conf{i}.pdbqt"
        if not receptor.exists():
            continue

        ligand = ligand_dir / f"ligand_c{i}_p0.pdbqt"
        if not ligand.exists():
            continue

        log = docking_dir / f"out_c{i}_p0.log"
        if is_complete(log):
            continue

        cx, cy, cz = pdbqt_centroid(ligand)
        jobs.append(
            (
                receptor.resolve(),
                ligand.resolve(),
                PocketCenter(x=cx, y=cy, z=cz),
                i,
                0,
                docking_dir,
            )
        )

    return jobs


def dock_one(args: tuple) -> None:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    receptor_path, ligand_path, center, conf_idx, pocket_idx, out_dir = args

    try:
        out_dir.mkdir(exist_ok=True, parents=True)
        log_file = out_dir / f"out_c{conf_idx}_p{pocket_idx}.log"
        if is_complete(log_file):
            print(f"Skipping docking for c{conf_idx}_p{pocket_idx} - completed.")
            return

        t0 = time.time()

        config = VinaDockingConfig(output_dir=out_dir.resolve(), box_size=(30, 30, 30))
        VinaDocking(config).run(
            DockingInput(
                receptor=PDBQT(path=receptor_path),
                ligand=PDBQT(path=ligand_path),
                pocket_center=center,
                conf_index=conf_idx,
                pose_index=pocket_idx,
            )
        )
        print(f"  c{conf_idx}_p{pocket_idx}: {time.time() - t0:.1f}s", flush=True)

    except Exception as e:
        if "outside the grid box" in str(e):
            out_dir.mkdir(exist_ok=True, parents=True)
            with open(out_dir / f"out_c{conf_idx}_p{pocket_idx}.log", "w") as f:
                f.write("SKIPPED: ligand outside grid box\nCOMPLETED\n")
        else:
            raise


def redock_one(protein_id: str, ligand_id: str) -> None:
    """Dock the pre-centered per-conf AF3-pocket pdbqts against each conformation's receptor."""
    protein_dir = DATA_DIR / protein_id
    docking_dir = protein_dir / "docking_af3_pocket" / ligand_id

    if (protein_dir / "docking_af3_pocket.tar.gz").exists():
        return

    print(f"Processing {protein_id} + {ligand_id}")
    jobs = gather_jobs(protein_id, ligand_id)
    if not jobs:
        print(f"  {protein_id}: all done, tarring...")
    else:
        print(f"  {protein_id}: {len(jobs)} docking jobs")
        with ProcessPoolExecutor(max_workers=MAX_DOCKING_WORKERS) as pool:
            futures = [pool.submit(dock_one, job) for job in jobs]
            for fut in tqdm(
                as_completed(futures), total=len(futures), desc=f"  {protein_id}"
            ):
                fut.result()

    if docking_dir.exists() and TAR:
        tar_path = protein_dir / "docking_af3_pocket.tar.gz"
        with tarfile.open(tar_path, "w:gz") as tar:
            tar.add(docking_dir, arcname=f"docking_af3_pocket/{ligand_id}")
        shutil.rmtree(docking_dir)
        print(f"  Tarred {protein_id}")


def main():
    args = load_args()
    pairs = filter_pairs()
    print(f"Running AF3-pocket pipeline for {len(pairs)} Ki/Kd pairs")

    af3_cfg = AF3Config.from_config(
        singularity_args=("--nv",),
        script_args=("--gpu_device=0",),
    )
    af3_inference = AF3Inference(
        af3_cfg,
        jax_cache_dir="/mnt/scratch/hegazyab/.cache",
        write_only=args.write_only,
        overwrite_input=args.overwrite_inp,
    )

    ready_to_dock = []
    for pid, lid, smi in pairs:
        af3_run_one(af3_inference, pid, lid, smi)
        if not (
            DATA_DIR
            / pid
            / "inference_af3_lig"
            / lid
            / "af_output"
            / lid
            / f"{lid}_model.cif"
        ).exists():
            continue
        ligand_prep_run_one(pid, lid, smi)
        ready_to_dock.append((pid, lid))
        # redock_one(pid, lid)

    with ProcessPoolExecutor(max_workers=MAX_PROTEIN_WORKERS) as pool:
        futures = [pool.submit(redock_one, pid, lid) for pid, lid in ready_to_dock]
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Proteins"):
            fut.result()

if __name__ == "__main__":
    main()
