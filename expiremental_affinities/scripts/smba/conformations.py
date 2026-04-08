from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from utils import load_proteins, DATA_DIR

from varidock.stages.vmd_frame_extract import (
    VMDFrameExtraction,
    VMDFrameExtractionConfig,
)
from varidock.types import PDB, PSF, Trajectory

MAX_WORKERS = 8
NUM_RUNS = 11



def process_protein(protein_id: str) -> None:
    protein_dir = DATA_DIR / protein_id
    eq_dir = protein_dir / "equilibration"
    conf_dir = protein_dir / "conformation"

    if (conf_dir / "protein_conf11.pdb").exists():
        return

    if not (eq_dir / f"system_run{NUM_RUNS - 1:03d}.coor").exists():
        return

    traj = Trajectory(
        psf=PSF(path=(eq_dir / "system.psf").resolve()),
        pdb=PDB(path=(eq_dir / "system.pdb").resolve()),
        coor_files=[
            (eq_dir / f"system_run{i:03d}.coor").resolve() for i in range(NUM_RUNS)
        ],
    )

    config = VMDFrameExtractionConfig(output_dir=conf_dir)
    VMDFrameExtraction(config).run(traj)


def main():
    proteins = load_proteins(DATA_DIR)

    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futures = [ex.submit(process_protein, pid) for pid in proteins]

        for fut in tqdm(
            as_completed(futures), total=len(futures), desc="Extracting frames"
        ):
            fut.result()


if __name__ == "__main__":
    main()
