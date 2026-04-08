"""
Compare Vina scores on crystal structures vs experimental dG.
"""

import csv
from pathlib import Path
import numpy as np
from scipy import stats

BENCHMARK = Path("../pdbbind/monomer_benchmark_small.csv")
VINA_DIR = Path("../data/pdbbind/vina_results")


def parse_best_affinity(log_file: Path) -> float | None:
    if not log_file.exists() or log_file.stat().st_size == 0:
        return None
    best = None
    in_table = False
    with open(log_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("-----+"):
                in_table = True
                continue
            if in_table:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        affinity = float(parts[1])
                        if best is None or affinity < best:
                            best = affinity
                    except ValueError:
                        in_table = False
                else:
                    in_table = False
    return best


exp_dG = []
vina_scores = []
failed = []

with open(BENCHMARK) as f:
    for row in csv.DictReader(f):
        pid = row["pdb_id"]
        log = VINA_DIR / pid / "out_c0_p0.log"
        score = parse_best_affinity(log)
        if score is None:
            failed.append(pid)
            continue
        exp_dG.append(float(row["dG_kcal_mol"]))
        vina_scores.append(score)

exp_dG = np.array(exp_dG)
vina_scores = np.array(vina_scores)

r, p = stats.pearsonr(exp_dG, vina_scores)
rho, _ = stats.spearmanr(exp_dG, vina_scores)
rmse = np.sqrt(np.mean((exp_dG - vina_scores) ** 2))

print(f"Entries: {len(exp_dG)} | Failed: {len(failed)}")
print(f"Pearson R:  {r:.3f} (p={p:.2e})")
print(f"Spearman ρ: {rho:.3f}")
print(f"RMSE:       {rmse:.2f} kcal/mol")
print(f"Mean exp dG:   {exp_dG.mean():.2f} ± {exp_dG.std():.2f}")
print(f"Mean Vina:     {vina_scores.mean():.2f} ± {vina_scores.std():.2f}")
