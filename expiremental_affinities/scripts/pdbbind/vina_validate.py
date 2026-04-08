"""
Run Vina on PDBbind crystal structures vs experimental dG.
1. ADFR prep protein -> pdbqt
2. Meeko prep ligand (centered at pocket COM) -> pdbqt
3. Vina dock
"""
import subprocess
import csv
import signal
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
# from rdkit import Chem
# 
# from varidock.stages.adfr_protein_receptor_prep import (
#     ADFRReceptorPrep,
#     ADFRReceptorPrepConfig,
# )
# from varidock.stages.meeko_ligand_prep import MeekoLigandPrep, MeekoLigandPrepConfig
# from varidock.stages.center_ligand_to_pocket import CenterLigand, CenterLigandConfig
from varidock.stages.vina_dock import VinaDocking, VinaDockingConfig
from varidock.types import (
    DockingInput,
    PocketCenter,
    PDBQT,
    # PDB,
    # Ligand,
    # LigandPrepInput,
)

PDBBIND = Path("/serviceberry/tank/abdolla/pdbbind_v2020/refined-set")
BENCHMARK = Path("../pdbbind/monomer_benchmark_small.csv")
OUT_DIR = Path("../data/pdbbind/vina_results")
MAX_WORKERS = 1


def get_ligand_com(mol2_file):
    xs, ys, zs = [], [], []
    with open(mol2_file) as f:
        in_atoms = False
        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atoms = True
                continue
            if line.startswith("@<TRIPOS>"):
                in_atoms = False
            if in_atoms:
                parts = line.split()
                if len(parts) >= 5:
                    xs.append(float(parts[2]))
                    ys.append(float(parts[3]))
                    zs.append(float(parts[4]))
    return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)

def prep_and_dock(row: dict) -> dict:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    pid = row["pdb_id"]
    entry_dir = PDBBIND / pid
    result_dir = OUT_DIR / pid
    result_dir.mkdir(exist_ok=True, parents=True)

    # Skip if already done
    log_file = result_dir / "out_c0_p0.log"
    if log_file.exists() and log_file.stat().st_size > 0:
        # print(f"{pid} already has Vina results, skipping")
        return {"pdb_id": pid, "status": "cached"}
    lig_x, lig_y, lig_z = get_ligand_com(entry_dir / f"{pid}_ligand.mol2")
    center = PocketCenter(x=lig_x, y=lig_y, z=lig_z)
    # center = PocketCenter(
    #     x=float(row["pocket_x"]),
    #     y=float(row["pocket_y"]),
    #     z=float(row["pocket_z"]),
    # )

    try:
        receptor_pdbqt = result_dir / f"{pid}_protein.pdbqt"
        if not receptor_pdbqt.exists():
            subprocess.run(
                # xr flag bc makes protein rigid?
                ["obabel", str(entry_dir / f"{pid}_protein.pdb"), "-O", str(receptor_pdbqt), "-xr", "-h"],
                check=True, capture_output=True,
            )
        
        # 1. ADFR prep protein
        # receptor_pdbqt = result_dir / f"{pid}_protein.pdbqt"
        # if not receptor_pdbqt.exists():
        #     config = ADFRReceptorPrepConfig(output_dir=result_dir.resolve())
        #     ADFRReceptorPrep(config).run(
        #         PDB(path=(entry_dir / f"{pid}_protein.pdb").resolve())
        #     )

        # 2. Prep ligand from mol2 -> pdbqt directly
        ligand_mol2 = entry_dir / f"{pid}_ligand.mol2"
        ligand_pdbqt = result_dir / f"{pid}_ligand.pdbqt"

        if not ligand_pdbqt.exists():
            subprocess.run(
                ["obabel", str(ligand_mol2), "-O", str(ligand_pdbqt), "-h"],
                check=True, capture_output=True,
            )

        # # 2. Prep ligand from mol2 -> pdb -> centered -> pdbqt
        # ligand_mol2 = entry_dir / f"{pid}_ligand.mol2"
        # ligand_pdb = result_dir / f"{pid}_ligand.pdb"
        # ligand_pdbqt = result_dir / "ligand_c0_p0.pdbqt"

        # if not ligand_pdbqt.exists():
        #     # mol2 -> pdb via obabel
        #     if not ligand_pdb.exists():
        #         subprocess.run(
        #             ["obabel", str(ligand_mol2), "-O", str(ligand_pdb), "-h"],
        #             check=True, capture_output=True,
        #         )

        #     inp = LigandPrepInput(
        #         ligand=Ligand(name=pid, pdb=PDB(path=ligand_pdb.resolve())),
        #         pocket_center=center,
        #         conf_index=0,
        #         pose_index=0,
        #     )

        #     # pdb -> pdbqt via meeko
        #     ligand_pdbqt = result_dir / f"{pid}_ligand.pdbqt"
        #     if not ligand_pdbqt.exists():
        #         inp = LigandPrepInput(
        #             ligand=Ligand(name=pid, pdb=PDB(path=ligand_pdb.resolve())),
        #             pocket_center=center,  # unused but required by type
        #             conf_index=0,
        #             pose_index=0,
        #         )
        #         meeko_config = MeekoLigandPrepConfig(output_dir=result_dir.resolve())
        #         MeekoLigandPrep(meeko_config).run(inp)
        #     # cleanup
        #     (result_dir / "ligand_c0_p0.pdb").unlink(missing_ok=True)
        #     (result_dir / "ligand_c0_p0_protonated.pdb").unlink(missing_ok=True)

        # 3. Vina dock
        config = VinaDockingConfig(
            output_dir=result_dir.resolve(), box_size=(30, 30, 30)
        )

        VinaDocking(config).run(
            DockingInput(
                receptor=PDBQT(path=receptor_pdbqt.resolve()),
                ligand=PDBQT(path=ligand_pdbqt.resolve()),
                pocket_center=center,
                conf_index=0,
                pose_index=0,
            )
        )


        return {"pdb_id": pid, "status": "done"}

    except Exception as e:
        return {"pdb_id": pid, "status": f"error: {e}"}


def main():
    OUT_DIR.mkdir(exist_ok=True, parents=True)

    rows = []
    with open(BENCHMARK) as f:
        for row in csv.DictReader(f):
            rows.append(row)

    print(f"Running Vina on {len(rows)} crystal structures")

    results = []
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(prep_and_dock, row): row for row in rows}
        for fut in tqdm(
            as_completed(futures), total=len(futures), desc="Vina validate"
        ):
            results.append(fut.result())

    done = sum(1 for r in results if r["status"] == "done")
    cached = sum(1 for r in results if r["status"] == "cached")
    failed = [r for r in results if r["status"] not in ("done", "cached")]
    print(f"Done: {done} | Cached: {cached} | Failed: {len(failed)}")
    for f in failed[:10]:
        print(f"  {f['pdb_id']}: {f['status']}")


if __name__ == "__main__":
    main()
