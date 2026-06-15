"""Submit eq STAGE 1 (the 2 ns run) to SLURM, throttled.

Stage 0 (system_eq000) is done locally first (run_local_eq.py). This queues eq.sh
for every sim that has system_eq000.coor but not yet system_eq001.coor -- one job
per sim (eq.sh runs system_eq.namd once, which continues from 000 -> 001). After
this, submit_runs.py picks up the dirs that now have system_eq001.coor.
"""

import subprocess
import time
from pathlib import Path

from tqdm import tqdm

from varidock.utils.slurm import _sbatch as sbatch, get_slurm_queue_count, get_job_name

DATA = Path("../../data/predictions")
SIMS = ["boltz", "smba", "smba_af3_pocket"]
STAGE0_DONE = "system_eq000.coor"
STAGE1_DONE = "system_eq001.coor"
MAX_QUEUED = 990
POLL_INTERVAL = 60


def sim_dirs() -> list[Path]:
    dirs = []
    for p in sorted(DATA.iterdir()):
        if not (p / "docking_af3_pocket").is_dir():
            continue
        for s in SIMS:
            d = p / "simulations" / s
            if (d / "eq.sh").exists():
                dirs.append(d)
    return dirs


def queued_names() -> set[str]:
    res = subprocess.run(["squeue", "--me", "-h", "--format=%j"],
                         capture_output=True, text=True, check=True)
    return set(res.stdout.split())


def main():
    running = queued_names()
    remaining = []
    no_stage0 = 0
    for d in sim_dirs():
        if not (d / STAGE0_DONE).exists():
            no_stage0 += 1            # run_local_eq.py hasn't done stage 0 here yet
            continue
        if (d / STAGE1_DONE).exists():
            continue                  # already equilibrated
        name = get_job_name(d / "eq.sh")
        if name and name in running:
            continue                  # already queued/running
        remaining.append(d)

    print(f"{len(remaining)} eq stage-1 jobs to submit "
          f"({no_stage0} still missing local stage 0 -- run run_local_eq.py)")

    with tqdm(total=len(remaining), desc="Submitting eq.sh") as pbar:
        while remaining:
            slots = MAX_QUEUED - get_slurm_queue_count()
            if slots <= 0:
                time.sleep(POLL_INTERVAL)
                continue
            batch, remaining = remaining[:slots], remaining[slots:]
            for d in batch:
                try:
                    sbatch(d / "eq.sh")
                    tqdm.write(f"  -> eq.sh @ {d.relative_to(DATA)}")
                except subprocess.CalledProcessError as e:
                    tqdm.write(f"  x {d.relative_to(DATA)}: {e.stderr or e}")
                pbar.update(1)


if __name__ == "__main__":
    main()
