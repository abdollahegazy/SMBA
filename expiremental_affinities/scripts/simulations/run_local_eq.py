"""Run eq STAGE 0 locally for every built sim.

eq is two namd stages (system_eq.namd, indexed by existing .coor count):
  stage 0 -> system_eq000  (minimize + brief warmup; quick) <-- THIS script, local
  stage 1 -> system_eq001  (2 ns; up to a 4h job)           <-- submit_eq.py, SLURM

So this only runs system_eq.namd once per dir, producing system_eq000.*. The 2 ns
stage 1 is queued separately (submit_eq.py), then production reads system_eq001.

Tune DEVICE / WORKERS to your local GPUs.
"""

import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

DATA = Path("../../data/predictions")
SIMS = ["boltz", "smba", "smba_af3_pocket"]
STAGE0_DONE = "system_eq000.coor"   # what one local run produces
DEVICE = "1"                        # local GPU id for namd +devices
WORKERS = 3                        # parallel sim dirs at once (one GPU each)


def sim_dirs() -> list[Path]:
    dirs = []
    for p in sorted(DATA.iterdir()):
        if not (p / "docking_af3_pocket").is_dir():
            continue
        for s in SIMS:
            d = p / "simulations" / s
            if (d / "system_eq.namd").exists():
                dirs.append(d)
    return dirs


def run_stage0(sim_dir: Path) -> str:
    rel = sim_dir.relative_to(DATA)
    if (sim_dir / STAGE0_DONE).exists():
        return f"skip  {rel} (stage 0 already done)"
    with open(sim_dir / "system_eq000.log", "w") as out, \
         open(sim_dir / "system_eq000.err", "w") as err:
        rc = subprocess.run(
            ["namd3", "+p8", "+setcpuaffinity", "+devices", DEVICE, "system_eq.namd"],
            cwd=sim_dir, stdout=out, stderr=err,
        ).returncode
    if rc != 0:
        return f"FAIL  {rel} (see system_eq000.err)"
    return (f"ok    {rel}" if (sim_dir / STAGE0_DONE).exists()
            else f"FAIL  {rel} (no {STAGE0_DONE} produced)")


def main():
    dirs = sim_dirs()
    print(f"{len(dirs)} sim dirs; running eq stage 0 locally "
          f"(DEVICE={DEVICE}, WORKERS={WORKERS})")
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        for msg in pool.map(run_stage0, dirs):
            print(msg, flush=True)


if __name__ == "__main__":
    main()
