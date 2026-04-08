from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from utils import load_proteins, DATA_DIR

from varidock.types import PDB
from varidock.stages.vmd_equil_prep import VMDEquilPrep, VMDEquilPrepConfig
from varidock.pipeline import Pipeline


MAX_WORKERS = 8
TOPPAR_DIR = Path("./smba/src/toppar/")
TEMPLATE_DIR = Path("./smba/src/equil/")




def xsc_is_zero(xsc_path: Path) -> bool:
    if not xsc_path.exists():
        return True

    with open(xsc_path) as f:
        lines = f.readlines()

    data = [line for line in lines if not line.startswith("#")]
    if not data:
        return True

    nums = list(map(float, data[-1].split()[1:]))
    return all(abs(v) < 1e-6 for v in nums[:9])


def process_protein(protein_id: str) -> None:
    protein_dir = DATA_DIR / protein_id
    pdb_path = protein_dir / f"{protein_id}.pdb"
    eq_dir = protein_dir / "equilibration"
    psf = eq_dir / "system.psf"
    xsc = eq_dir / "system.xsc"

    if not pdb_path.exists():
        return

    if not xsc_is_zero(xsc) and psf.exists():
        return

    eq_dir.mkdir(parents=True, exist_ok=True)
    config = VMDEquilPrepConfig(
        toppar_dir=TOPPAR_DIR,
        template_dir=TEMPLATE_DIR,
        output_dir=eq_dir,
    )

    Pipeline(VMDEquilPrep(config)).run(PDB(path=pdb_path))


def main():
    proteins = load_proteins()

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = [ex.submit(process_protein, pid) for pid in proteins]

        for fut in tqdm(
            as_completed(futures), total=len(futures), desc="Prepping equilibration"
        ):
            fut.result()


if __name__ == "__main__":
    main()
