# all claude'ed

"""Compare SMBA (Vina) and Boltz against experimental Ki/Kd. Saves figures."""

import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_CSV = REPO_ROOT / "data" / "docking_results.csv"
FIG_DIR = REPO_ROOT / "data" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

LIGAND_NAMES = {
    "7930": "Ergotamine",
    "4450907": "Morphine",
    "10194105": "Atropine",
    "2424": "Caffeine",
    "84989": "Quinine",
    "917": "Nicotine",
}

SMBA_COLOR = "#888780"   # warm gray
BOLTZ_COLOR = "#534AB7"  # purple
TRUTH_COLOR = "#1D9E75"  # teal


def as_float(value: str) -> float | None:
    value = value.strip()
    if not value or value == "...":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def pearson(xs, ys):
    if len(xs) < 2:
        return None
    xs = np.asarray(xs); ys = np.asarray(ys)
    if xs.std() == 0 or ys.std() == 0:
        return None
    return float(np.corrcoef(xs, ys)[0, 1])


def spearman(xs, ys):
    if len(xs) < 2:
        return None
    return pearson(_ranks(xs), _ranks(ys))


def _ranks(values):
    order = np.argsort(values)
    ranks = np.empty_like(order, dtype=float)
    ranks[order] = np.arange(1, len(values) + 1)
    return ranks


def load_rows():
    """Returns dict: ligand_id -> list of (truth, vina, boltz). Only ki/kd rows."""
    by_ligand = defaultdict(list)
    with open(RESULTS_CSV, newline="") as f:
        for row in csv.DictReader(f):
            if row["known_affinity_type"].strip() not in {"ki_nM", "kd_nM"}:
                continue
            truth = as_float(row["known_affinity_kcal_mol"])
            if truth is None:
                continue
            vina = as_float(row["vina_affinity"])
            boltz = as_float(row["boltz_affinity"])
            by_ligand[row["ligand"].strip()].append((truth, vina, boltz))
    return by_ligand


def compute_metrics(rows, idx):
    """rows is list of (truth, vina, boltz). idx=1 for vina, 2 for boltz."""
    pairs = [(r[0], r[idx]) for r in rows if r[idx] is not None]
    if not pairs:
        return None
    truth = np.array([p[0] for p in pairs])
    pred = np.array([p[1] for p in pairs])
    err = pred - truth
    return {
        "n": len(pairs),
        "mae": float(np.mean(np.abs(err))),
        "rmse": float(np.sqrt(np.mean(err**2))),
        "bias": float(np.mean(err)),
        "pearson": pearson(pred, truth),
        "spearman": spearman(pred, truth),
        "avg_truth": float(np.mean(truth)),
    }


# ---------- Figure 1: per-ligand metrics side-by-side ----------
def figure_metrics(by_ligand):
    # Sort by avg truth (strongest binders first = most negative)
    ligand_order = sorted(
        by_ligand,
        key=lambda lid: np.mean([r[0] for r in by_ligand[lid]]),
    )
    names = [LIGAND_NAMES.get(lid, lid) for lid in ligand_order]
    smba_m = [compute_metrics(by_ligand[lid], 1) for lid in ligand_order]
    boltz_m = [compute_metrics(by_ligand[lid], 2) for lid in ligand_order]

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    fig.suptitle("SMBA vs Boltz — per ligand (Ki/Kd subset)", fontsize=14, fontweight="500")

    metric_keys = [("mae", "MAE (kcal/mol)"), ("rmse", "RMSE (kcal/mol)"),
                   ("bias", "Bias (kcal/mol)"), ("spearman", "Spearman ρ")]

    x = np.arange(len(names))
    width = 0.38

    for ax, (key, label) in zip(axes.flat, metric_keys):
        smba_vals = [m[key] if m and m[key] is not None else np.nan for m in smba_m]
        boltz_vals = [m[key] if m and m[key] is not None else np.nan for m in boltz_m]
        ax.bar(x - width/2, smba_vals, width, label="SMBA", color=SMBA_COLOR)
        ax.bar(x + width/2, boltz_vals, width, label="Boltz", color=BOLTZ_COLOR)
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=25, ha="right", fontsize=9)
        ax.set_title(label, fontsize=11, fontweight="500")
        ax.axhline(0, color="#000", linewidth=0.5, alpha=0.3)
        ax.grid(axis="y", alpha=0.2, linewidth=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", labelsize=9)

    axes[0, 0].legend(loc="upper left", frameon=False, fontsize=10)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIG_DIR / "per_ligand_metrics.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------- Figure 2: predicted vs truth scatter, per ligand ----------
def figure_scatter(by_ligand):
    ligand_order = sorted(
        by_ligand,
        key=lambda lid: np.mean([r[0] for r in by_ligand[lid]]),
    )

    n = len(ligand_order)
    cols = 3
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(12, 4 * rows), sharex=False, sharey=False)
    axes = np.atleast_2d(axes).flatten()
    fig.suptitle("Predicted vs experimental ΔG — per ligand", fontsize=14, fontweight="500")

    for ax, lid in zip(axes, ligand_order):
        data = by_ligand[lid]
        truth = np.array([r[0] for r in data])
        vina_pts = [(r[0], r[1]) for r in data if r[1] is not None]
        boltz_pts = [(r[0], r[2]) for r in data if r[2] is not None]

        if vina_pts:
            v_t, v_p = zip(*vina_pts)
            rho = spearman(list(v_p), list(v_t))
            label = f"SMBA (ρ={rho:.2f})" if rho is not None else "SMBA"
            ax.scatter(v_t, v_p, color=SMBA_COLOR, s=42, alpha=0.8,
                       label=label, edgecolor="white", linewidth=0.6)
        if boltz_pts:
            b_t, b_p = zip(*boltz_pts)
            rho = spearman(list(b_p), list(b_t))
            label = f"Boltz (ρ={rho:.2f})" if rho is not None else "Boltz"
            ax.scatter(b_t, b_p, color=BOLTZ_COLOR, s=42, alpha=0.8,
                       label=label, edgecolor="white", linewidth=0.6)

        # y = x reference line
        all_vals = list(truth)
        for r in data:
            if r[1] is not None: all_vals.append(r[1])
            if r[2] is not None: all_vals.append(r[2])
        lo, hi = min(all_vals) - 0.5, max(all_vals) + 0.5
        ax.plot([lo, hi], [lo, hi], "--", color="#000", alpha=0.25, linewidth=0.8, label="y = x")

        name = LIGAND_NAMES.get(lid, lid)
        ax.set_title(f"{name} (n={len(data)})", fontsize=11, fontweight="500")
        ax.set_xlabel("Experimental ΔG", fontsize=9)
        ax.set_ylabel("Predicted ΔG", fontsize=9)
        ax.legend(fontsize=8, frameon=False, loc="best")
        ax.grid(alpha=0.2, linewidth=0.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", labelsize=8)

    for ax in axes[n:]:
        ax.set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIG_DIR / "scatter_per_ligand.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------- Figure 3: overall summary ----------
def figure_overall(by_ligand):
    all_rows = [r for rows in by_ligand.values() for r in rows]
    smba = compute_metrics(all_rows, 1)
    boltz = compute_metrics(all_rows, 2)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    fig.suptitle(f"Pooled across all Ki/Kd (n={smba['n']})", fontsize=13, fontweight="500")

    # Left: bar chart metrics
    ax = axes[0]
    labels = ["RMSE", "Bias", "Spearman ρ"]
    smba_vals = [smba["rmse"], smba["bias"], smba["spearman"]]
    boltz_vals = [boltz["rmse"], boltz["bias"], boltz["spearman"]]
    x = np.arange(len(labels))
    width = 0.38
    ax.bar(x - width/2, smba_vals, width, label="SMBA", color=SMBA_COLOR)
    ax.bar(x + width/2, boltz_vals, width, label="Boltz", color=BOLTZ_COLOR)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2, linewidth=0.5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.set_title("Metrics", fontsize=11, fontweight="500")

    # Right: scatter pooled
    ax = axes[1]
    pairs_v = [(r[0], r[1]) for r in all_rows if r[1] is not None]
    pairs_b = [(r[0], r[2]) for r in all_rows if r[2] is not None]
    if pairs_v:
        t, p = zip(*pairs_v)
        ax.scatter(t, p, color=SMBA_COLOR, s=30, alpha=0.65, label="SMBA", edgecolor="white", linewidth=0.5)
    if pairs_b:
        t, p = zip(*pairs_b)
        ax.scatter(t, p, color=BOLTZ_COLOR, s=30, alpha=0.65, label="Boltz", edgecolor="white", linewidth=0.5)
    all_v = [r[0] for r in all_rows] + [r[1] for r in all_rows if r[1] is not None] + [r[2] for r in all_rows if r[2] is not None]
    lo, hi = min(all_v) - 0.5, max(all_v) + 0.5
    ax.plot([lo, hi], [lo, hi], "--", color="#000", alpha=0.25, linewidth=0.8, label="y = x")
    ax.set_xlabel("Experimental ΔG (kcal/mol)")
    ax.set_ylabel("Predicted ΔG (kcal/mol)")
    ax.set_title("Pooled scatter", fontsize=11, fontweight="500")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.2, linewidth=0.5)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIG_DIR / "overall_summary.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def print_summary(by_ligand):
    all_rows = [r for rows in by_ligand.values() for r in rows]
    smba = compute_metrics(all_rows, 1)
    boltz = compute_metrics(all_rows, 2)

    for name, m in [("SMBA vs truth", smba), ("Boltz vs truth", boltz)]:
        print(name)
        print(f"  n:        {m['n']}")
        print(f"  MAE:      {m['mae']:.3f} kcal/mol")
        print(f"  RMSE:     {m['rmse']:.3f} kcal/mol")
        print(f"  Bias:     {m['bias']:.3f} kcal/mol")
        print(f"  Pearson:  {m['pearson']:.3f}")
        print(f"  Spearman: {m['spearman']:.3f}")
        print()


def main():
    by_ligand = load_rows()
    print_summary(by_ligand)
    p1 = figure_overall(by_ligand)
    p2 = figure_metrics(by_ligand)
    p3 = figure_scatter(by_ligand)
    print("Saved:")
    print(f"  {p1}")
    print(f"  {p2}")
    print(f"  {p3}")


if __name__ == "__main__":
    main()