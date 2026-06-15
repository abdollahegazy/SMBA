"""Throttled production submitter for the af3_pocket sims.

For every sim whose eq is done (system_eq001.coor) and that hasn't reached
TARGET_NS yet, submit BOTH run.sh (normal nodes) and run_gh.sh (grace nodes) --
can't request both node types in one job, so we queue both and whichever starts
first wins (the scripts' safety check kills the loser). Keeps the queue under
MAX_QUEUED and loops until every sim hits TARGET_NS.

ns are read from the latest system_run*.log via varidock's namd parser (2 fs).
"""

import subprocess
import time
from pathlib import Path

from tqdm import tqdm

from varidock.utils.slurm import _sbatch as sbatch, get_slurm_queue_count, get_job_name
from varidock.utils.namd import get_namd_ns

DATA = Path("../../data/predictions")
SIMS = ["boltz", "smba", "smba_af3_pocket"]
EQ_DONE = "system_eq001.coor"
TARGET_NS = 250.0            # 125e6 steps * 2 fs in system_run.namd
TIMESTEP_FS = 2.0
MAX_QUEUED = 990
POLL_INTERVAL = 20
MAX_ATTEMPTS = 100          # per-dir safety so a wedged sim can't loop forever


def sim_dirs() -> list[Path]:
    dirs = []
    for p in sorted(DATA.iterdir()):
        if not (p / "docking_af3_pocket").is_dir():
            continue
        for s in SIMS:
            d = p / "simulations" / s
            if (d / "run.sh").exists():
                dirs.append(d)
    return dirs


def total_ns(sim_dir: Path) -> float:
    logs = sorted(sim_dir.glob("system_run[0-9][0-9][0-9].log"))
    if not logs:
        return 0.0
    res = get_namd_ns(logs[-1], timestep_fs=TIMESTEP_FS)
    return res[0] if res else 0.0


def needs_more(sim_dir: Path) -> bool:
    return total_ns(sim_dir) < TARGET_NS


def pending_or_running_names() -> set[str]:
    res = subprocess.run(["squeue", "--me", "-h", "--format=%j"],
                         capture_output=True, text=True, check=True)
    return set(res.stdout.split())


def main():
    all_dirs = sim_dirs()
    remaining = []
    for d in tqdm(all_dirs, desc="Checking"):
        if not (d / EQ_DONE).exists():
            continue          # eq not finished -> skip (run_local_eq.py first)
        if needs_more(d):
            remaining.append(d)

    print(f"{len(remaining)}/{len(all_dirs)} sims need more sampling (target {TARGET_NS:.0f} ns)")
    attempts = dict.fromkeys(remaining, 0)

    with tqdm(total=len(remaining), desc="Submitting") as pbar:
        while remaining:
            ns_now = sum(total_ns(d) for d in remaining)
            pbar.set_postfix_str(f"{ns_now:.0f}/{len(remaining) * TARGET_NS:.0f} ns")

            slots = max(0, MAX_QUEUED - get_slurm_queue_count())
            if slots <= 0:
                time.sleep(POLL_INTERVAL)
                continue

            running = pending_or_running_names()
            submitted = 0
            i = 0
            while i < len(remaining) and submitted < slots:
                d = remaining[i]

                if not needs_more(d):
                    remaining.pop(i)
                    pbar.update(1)
                    continue

                name = get_job_name(d / "run.sh")
                name_gh = get_job_name(d / "run_gh.sh")
                if (name and name in running) or (name_gh and name_gh in running):
                    i += 1
                    continue

                for script, jobname in [("run.sh", name), ("run_gh.sh", name_gh)]:
                    sp = d / script
                    if not sp.exists() or (jobname and jobname in running):
                        continue
                    try:
                        sbatch(sp)
                        submitted += 1
                        if jobname:
                            running.add(jobname)
                        tqdm.write(f"  -> {script} @ {d.relative_to(DATA)}")
                    except subprocess.CalledProcessError as e:
                        tqdm.write(f"  x {sp}: {e.stderr or e}")

                attempts[d] += 1
                if attempts[d] >= MAX_ATTEMPTS:
                    tqdm.write(f"max attempts for {d.relative_to(DATA)}, dropping")
                    remaining.pop(i)
                    pbar.update(1)
                    continue
                i += 1

            time.sleep(POLL_INTERVAL)


if __name__ == "__main__":
    main()
