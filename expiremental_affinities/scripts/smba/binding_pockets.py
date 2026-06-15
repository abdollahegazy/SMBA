# had to export LD_LIBRARY_PATH=$(python3 -c "import os,nvidia,glob; b=os.path.dirname(nvidia.__file__); print(':'.join(glob.glob(b+'/*/lib')))"):$LD_LIBRARY_PATH
# see https://claude.ai/share/366b06e4-d054-45e7-abdc-2c2cd8ee1934

from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from utils import load_proteins, DATA_DIR, _worker_init
import signal
import os
import traceback


os.environ["TF_NUM_INTRAOP_THREADS"] = "2"
os.environ["TF_NUM_INTEROP_THREADS"] = "2"
os.environ["OMP_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"


MAX_WORKERS = 1
NUM_CONFS = 11


def is_deepsurf_done(protein_id: str) -> bool:
    pocket_dir = DATA_DIR / protein_id / "pockets"
    return all(
        (pocket_dir / f"protein_conf{i}" / "centers.txt").exists()
        for i in range(NUM_CONFS)
    )


def process_protein(protein_id: str) -> None:
    from varidock.types import PDB
    from varidock.stages.deepsurf_pockets import DeepSurfPocketConfig, DeepSurfPockets

    protein_dir = DATA_DIR / protein_id
    conf_dir = protein_dir / "conformation"
    pocket_dir = protein_dir / "pockets"
    print(conf_dir)
    if not conf_dir.exists():
        return

    pocket_dir.mkdir(parents=True, exist_ok=True)

    for conf in conf_dir.iterdir():
        if not conf.name.endswith(".pdb"):
            continue

        if (pocket_dir / conf.stem / "centers.txt").exists():
            continue


        config = DeepSurfPocketConfig(output_dir=pocket_dir.resolve())
        DeepSurfPockets(config).run(PDB(path=conf.resolve()))


def main():
    proteins = load_proteins(DATA_DIR)

    remaining = [pid for pid in proteins if not is_deepsurf_done(pid)]
    print(f"Processing {len(remaining)}/{len(proteins)} proteins")

    pool = ProcessPoolExecutor(max_workers=MAX_WORKERS, initializer=_worker_init,max_tasks_per_child=1)
    futures = {pool.submit(process_protein, pid): pid for pid in remaining}
    try:
        for fut in tqdm(as_completed(futures), total=len(futures), desc="DeepSurf"):
            try:
                fut.result()
            except Exception as e:
                print(f"ERROR: {e}")
                traceback.print_exc()
    except KeyboardInterrupt:
        for fut in futures:
            fut.cancel()
        for pid in list(pool._processes):
            try:
                os.killpg(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass


if __name__ == "__main__":
    main()
