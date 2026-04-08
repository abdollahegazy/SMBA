import time
from tqdm import tqdm
from varidock.utils.slurm import _sbatch as sbatch, get_slurm_queue_count, get_running_job_names, get_job_name
from varidock.utils.namd import is_namd_done
from utils import load_proteins, DATA_DIR

MAX_QUEUED = 990
POLL_INTERVAL = 60
TARGET_NS = 1.002


def main():
    proteins = load_proteins()

    all_dirs = []
    for pid in proteins:
        eq_dir = DATA_DIR / pid / "equilibration"
        if eq_dir.exists():
            all_dirs.append(eq_dir)

    running = get_running_job_names()

    remaining = []
    for eq_dir in tqdm(all_dirs, desc="Checking"):
        script = eq_dir / "eq.sh"
        name = get_job_name(script)
        if name and name in running:
            continue

        if is_namd_done(eq_dir / "system_eq000.log", TARGET_NS):
            continue

        remaining.append(eq_dir)

    print(f"{len(remaining)}/{len(all_dirs)} need submission")

    with tqdm(total=len(remaining), desc="Submitting eq.sh") as pbar:
        while remaining:
            current = get_slurm_queue_count()
            slots = MAX_QUEUED - current

            if slots <= 0:
                time.sleep(POLL_INTERVAL)
                continue

            batch = remaining[:slots]
            remaining = remaining[slots:]

            for eq_dir in batch:
                sbatch(eq_dir / "eq.sh")
                pbar.update(1)


if __name__ == "__main__":
    main()