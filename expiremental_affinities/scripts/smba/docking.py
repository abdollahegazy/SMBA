from pathlib import Path
import tarfile
import shutil
import signal
import time
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from utils import load_pairs, DATA_DIR
from varidock.stages.vina_dock import VinaDocking, VinaDockingConfig
from varidock.types import DockingInput, PocketCenter, PDBQT


MAX_WORKERS = 1
NUM_CONFS = 11
TAR = False

def is_complete(log_file: Path) -> bool:
    if not log_file.exists() or log_file.stat().st_size == 0:
        return False
    with open(log_file) as f:
        return f.read().rstrip().endswith("COMPLETED")


def dock_one(args: tuple) -> None:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    receptor_path, ligand_path, center, conf_idx, pocket_idx, out_dir = args

    try:
        out_dir.mkdir(exist_ok=True, parents=True)
        log_file = out_dir / f"out_c{conf_idx}_p{pocket_idx}.log"
        if log_file.exists() and log_file.stat().st_size > 0:
            return

        t0 = time.time()

        config = VinaDockingConfig(output_dir=out_dir.resolve(),box_size=(30,30,30))
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


def gather_jobs(protein_id: str, ligand_id: str) -> list[tuple]:
    protein_dir = DATA_DIR / protein_id
    pocket_dir = protein_dir / "pockets"
    ligand_dir = protein_dir / "ligands_prepared" / ligand_id
    receptor_dir = protein_dir / "protein_receptors"
    docking_dir = protein_dir / "docking" / ligand_id

    jobs = []
    for i in range(NUM_CONFS):
        centers_file = pocket_dir / f"protein_conf{i}" / "centers.txt"
        if not centers_file.exists() or centers_file.stat().st_size == 0:
            continue

        receptor = receptor_dir / f"protein_conf{i}.pdbqt"
        if not receptor.exists():
            continue

        with open(centers_file) as f:
            for pocket, line in enumerate(f):
                x, y, z = map(float, line.strip().split())
                # x,y,z = -32.9272, 29.116, -6.38

                ligand = ligand_dir / f"ligand_c{i}_p{pocket}.pdbqt"

                if not ligand.exists():
                    continue

                log = docking_dir / f"out_c{i}_p{pocket}.log"
                if is_complete(log):
                    continue

                jobs.append(
                    (
                        receptor.resolve(),
                        ligand.resolve(),
                        PocketCenter(x=x, y=y, z=z),
                        i,
                        pocket,
                        docking_dir,
                    )
                )

    return jobs


def process_protein(protein_id: str, ligand_id: str) -> None:
    protein_dir = DATA_DIR / protein_id
    docking_dir = protein_dir / "docking" / ligand_id

    if (protein_dir / "docking.tar.gz").exists():
        return

    jobs = gather_jobs(protein_id, ligand_id)
    if not jobs:
        print(f"  {protein_id}: all done, tarring...")
    else:
        print(f"  {protein_id}: {len(jobs)} docking jobs")
        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
            futures = [pool.submit(dock_one, job) for job in jobs]
            for fut in tqdm(
                as_completed(futures), total=len(futures), desc=f"  {protein_id}"
            ):
                fut.result()

    if docking_dir.exists() and TAR:
        tar_path = protein_dir / "docking.tar.gz"
        with tarfile.open(tar_path, "w:gz") as tar:
            tar.add(docking_dir, arcname=f"docking/{ligand_id}")
        shutil.rmtree(docking_dir)
        print(f"  Tarred {protein_id}")


def main():
    pairs = load_pairs()
    for protein_id, ligand_id, _ in pairs:
        print(f"Processing {protein_id} + {ligand_id}")
        process_protein(protein_id, ligand_id)


if __name__ == "__main__":
    main()
