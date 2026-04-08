import time
import subprocess
from pathlib import Path
from tqdm import tqdm
from utils import load_proteins, DATA_DIR

from varidock.utils.slurm import (
    _sbatch as sbatch,
    get_slurm_queue_count,
    get_running_job_names,
    get_job_name,
)


MAX_QUEUED = 990
POLL_INTERVAL = 60
NUM_RUNS = 11




def next_incomplete_run_idx(eq_dir: Path) -> int | None:
    for i in range(NUM_RUNS):
        if not (eq_dir / f"system_run{i:03d}.coor").exists():
            return i
    return None


def main():
    proteins = load_proteins()

    all_dirs = []
    for pid in proteins:
        eq_dir = DATA_DIR / pid / "equilibration"
        if eq_dir.exists():
            all_dirs.append(eq_dir)

    remaining = []
    remaining_calls: dict[Path, int] = {}

    for eq_dir in tqdm(all_dirs, desc="Checking"):
        if not (eq_dir / "system_eq001.coor").exists():
            continue

        if (eq_dir / f"system_run{NUM_RUNS - 1:03d}.coor").exists():
            continue

        nxt = next_incomplete_run_idx(eq_dir)
        if nxt is None:
            continue

        remaining_calls[eq_dir] = NUM_RUNS - nxt
        remaining.append(eq_dir)

    total_submissions = sum(remaining_calls.values())
    print(
        f"{len(remaining)}/{len(all_dirs)} dirs need work; {total_submissions} submissions remaining"
    )

    attempted = 0

    with tqdm(total=total_submissions, desc="Submitting run.sh") as pbar:
        while remaining and attempted < total_submissions:
            current = get_slurm_queue_count()
            slots = max(0, MAX_QUEUED - current)

            if slots <= 0:
                time.sleep(POLL_INTERVAL)
                continue

            running = get_running_job_names()
            submitted = 0
            i = 0

            while i < len(remaining) and submitted < slots:
                eq_dir = remaining[i]

                nxt = next_incomplete_run_idx(eq_dir)
                if nxt is None:
                    remaining.pop(i)
                    continue

                name = get_job_name(eq_dir / "run.sh")
                if name and name in running:
                    i += 1
                    continue

                try:
                    sbatch(eq_dir / "run.sh")
                    if name:
                        running.add(name)
                except subprocess.CalledProcessError as e:
                    print(f"sbatch failed for {eq_dir}: {e.stderr or e}")
                    i += 1
                    continue

                pbar.update(1)
                submitted += 1
                attempted += 1
                i += 1


if __name__ == "__main__":
    main()
