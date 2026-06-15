from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from utils import load_proteins, DATA_DIR

from varidock.types import CIF
from varidock.stages.cif_to_pdb import CIFToPDB, CIFToPDBConfig
from varidock.stages.insert_plddt_to_pdb import InsertPLDDT, InsertPLDDTConfig
from varidock.pipeline import Pipeline

MAX_WORKERS = 16


def process_protein(protein_id: str) -> None:
    pdb_out = DATA_DIR / protein_id / f"{protein_id}.pdb"
    if pdb_out.exists():
        print (f"{protein_id} already has PDB, skipping")
        return

    cif_path = (
        DATA_DIR
        / protein_id
        / "inference_af3"
        / "af_output"
        / protein_id
        / f"{protein_id}_model.cif"
    )
    if not cif_path.exists():
        return

    out_dir = DATA_DIR / protein_id
    config_cif = CIFToPDBConfig(output_dir=out_dir)
    config_plddt = InsertPLDDTConfig(output_dir=out_dir)

    Pipeline(CIFToPDB(config_cif), InsertPLDDT(config_plddt)).run(CIF(path=cif_path))


def main():
    proteins = load_proteins(DATA_DIR)

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = [ex.submit(process_protein, pid) for pid in proteins]

        for fut in tqdm(as_completed(futures), total=len(futures), desc="CIF → PDB"):
            fut.result()


if __name__ == "__main__":
    main()
