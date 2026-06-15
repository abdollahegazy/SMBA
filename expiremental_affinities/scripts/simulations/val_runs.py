"""Report eq + production progress across all af3_pocket sims.

eq stage0 = system_eq000.coor present (local minimize+warmup)
eq stage1 = system_eq001.coor present (SLURM 2 ns; what production reads)
prod ns   = ns of the latest system_run*.log (2 fs timestep)
finished  = prod ns >= TARGET_NS
"""

from pathlib import Path

from varidock.utils.namd import get_namd_ns

DATA = Path("../../data/predictions")
SIMS = ["boltz", "smba", "smba_af3_pocket"]
STAGE0_DONE = "system_eq000.coor"
EQ_DONE = "system_eq001.coor"
TARGET_NS = 250.0
TIMESTEP_FS = 2.0


def total_ns(sim_dir: Path) -> float:
    logs = sorted(sim_dir.glob("system_run[0-9][0-9][0-9].log"))
    if not logs:
        return 0.0
    res = get_namd_ns(logs[-1], timestep_fs=TIMESTEP_FS)
    return res[0] if res else 0.0


def main():
    total = stage0 = eq_done = finished = 0
    ns_sum = 0.0
    no_stage0, no_eq, in_progress = [], [], []

    for p in sorted(DATA.iterdir()):
        if not (p / "docking_af3_pocket").is_dir():
            continue
        for s in SIMS:
            d = p / "simulations" / s
            if not d.exists():
                continue
            total += 1
            rel = d.relative_to(DATA)
            if (d / STAGE0_DONE).exists():
                stage0 += 1
            if not (d / EQ_DONE).exists():
                (no_eq if (d / STAGE0_DONE).exists() else no_stage0).append(str(rel))
                continue
            eq_done += 1
            ns = total_ns(d)
            ns_sum += min(ns, TARGET_NS)
            if ns >= TARGET_NS:
                finished += 1
            else:
                in_progress.append((str(rel), ns))

    print(f"sims:               {total}")
    print(f"eq stage0 (local):  {stage0}/{total}")
    print(f"eq stage1 (slurm):  {eq_done}/{total}")
    print(f"production done:    {finished}/{total}  (>= {TARGET_NS:.0f} ns)")
    print(f"aggregate ns:       {ns_sum:.0f} / {total * TARGET_NS:.0f}")

    if no_stage0:
        print(f"\nno local stage 0 yet ({len(no_stage0)}) -- run run_local_eq.py:")
        for r in no_stage0[:20]:
            print(f"  {r}")
        if len(no_stage0) > 20:
            print(f"  ... +{len(no_stage0) - 20} more")

    if no_eq:
        print(f"\nstage 0 done, awaiting slurm stage 1 ({len(no_eq)}) -- run submit_eq.py:")
        for r in no_eq[:20]:
            print(f"  {r}")
        if len(no_eq) > 20:
            print(f"  ... +{len(no_eq) - 20} more")

    if in_progress:
        print(f"\nin progress ({len(in_progress)}):")
        for r, ns in sorted(in_progress, key=lambda x: x[1])[:20]:
            print(f"  {ns:6.1f} ns  {r}")
        if len(in_progress) > 20:
            print(f"  ... +{len(in_progress) - 20} more")


if __name__ == "__main__":
    main()
