# export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
from utils import load_proteins, DATA_DIR

from varidock.stages.adfr_protein_receptor_prep import (
    ADFRReceptorPrepConfig,
    ADFRReceptorPrep,
)
from varidock.types import PDB


MAX_WORKERS = 6
NUM_CONFS = 11



def gather_conformations(proteins: list[str]) -> list[PDB]:
    conformations = []

    for pid in tqdm(proteins, desc="Checking conformations"):
        protein_dir = DATA_DIR / pid

        for i in range(NUM_CONFS):
            conf_pdb = protein_dir / "conformation" / f"protein_conf{i}.pdb"
            if not conf_pdb.exists():
                continue

            centers = protein_dir / "pockets" / f"protein_conf{i}" / "centers.txt"
            if not centers.exists() or centers.stat().st_size == 0:
                continue

            pdbqt = protein_dir / "protein_receptors" / f"protein_conf{i}.pdbqt"
            if pdbqt.exists() and pdbqt.stat().st_size > 0:
                continue

            conformations.append(PDB(path=conf_pdb.resolve()))

    return conformations


def process_one(protein: PDB) -> None:
    out_dir = protein.path.parent.parent / "protein_receptors"
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        config = ADFRReceptorPrepConfig(output_dir=out_dir.resolve())
        ADFRReceptorPrep(config).run(protein)
    except Exception as e:
        print(f"ERROR: {protein.path}: {e}")
        raise


def main():
    proteins = load_proteins()
    conformations = gather_conformations(proteins)

    print(f"{len(conformations)} conformations to process")

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = [pool.submit(process_one, conf) for conf in conformations]

        for fut in tqdm(as_completed(futures), total=len(futures), desc="ADFR prep"):
            try:
                fut.result()
            except Exception as e:
                print(f"ERROR: {e}")


if __name__ == "__main__":
    main()
