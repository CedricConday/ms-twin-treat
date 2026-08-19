"""How much of Kang's measured IFN-beta delta is even predictable?

A score is only interpretable against what a *perfect* model could score. The
harness measures a model against the delta we happened to measure, and that
measurement carries sampling noise — badly so where a cell type is rare. This
script estimates the ceiling directly from the data, with no model involved:

  split the control cells into two halves and the stimulated cells into two
  halves, compute a delta from each half, and correlate the two. That
  split-half r is the reliability of a half-sample delta; Spearman-Brown lifts
  it to the full sample, and sqrt() of that is the highest delta_pearson any
  model could reach against the observed delta.

The headline it produces: Megakaryocytes (63 control / 69 stimulated cells) have
a split-half reliability near zero — that fold is scoring noise, not biology,
and it drags every model's aggregate down by roughly 0.05. Report the aggregate,
but read the per-cell-type table next to this ceiling before believing any of it.

Runs in the main pinned `.venv` (needs scanpy only).

    . .venv/bin/activate && python scripts/noise_ceiling.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import scanpy as sc
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data.kang import fetch_kang  # noqa: E402


def split_half_reliability(X: np.ndarray, ctrl_idx: np.ndarray, stim_idx: np.ndarray,
                           n_rep: int, rng: np.random.Generator) -> float:
    rs = []
    for _ in range(n_rep):
        ca, cb = np.array_split(rng.permutation(ctrl_idx), 2)
        sa, sb = np.array_split(rng.permutation(stim_idx), 2)
        da = X[sa].mean(axis=0) - X[ca].mean(axis=0)
        db = X[sb].mean(axis=0) - X[cb].mean(axis=0)
        rs.append(stats.pearsonr(da, db).statistic)
    return float(np.mean(rs))


def main(min_cells: int = 50, n_rep: int = 40, seed: int = 0) -> None:
    adata = fetch_kang().copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    X = np.asarray(adata.X.todense() if hasattr(adata.X, "todense") else adata.X,
                   dtype=np.float64)
    ct = adata.obs["cell_type"].values
    lab = adata.obs["label"].values
    rng = np.random.default_rng(seed)

    print("Noise ceiling on the Kang IFN-beta delta (no model involved).\n")
    print(f"{'cell type':<20} {'n_ctrl':>7} {'n_stim':>7} {'split-half r':>13} {'ceiling':>9}")
    ceilings = []
    for c in sorted({str(v) for v in ct}):
        ci = np.where((ct == c) & (lab == "ctrl"))[0]
        si = np.where((ct == c) & (lab == "stim"))[0]
        if len(ci) < min_cells or len(si) < min_cells:
            continue
        rh = split_half_reliability(X, ci, si, n_rep, rng)
        r_full = 2 * rh / (1 + rh)                    # Spearman-Brown
        ceil = float(np.sqrt(max(r_full, 0.0)))
        ceilings.append(ceil)
        print(f"{c:<20} {len(ci):>7} {len(si):>7} {rh:>13.3f} {ceil:>9.3f}")

    print(f"\nmean attainable delta_pearson (perfect model) = {np.mean(ceilings):.4f}")
    print("\nRead this before quoting any aggregate: a fold whose split-half r is "
          "near zero\ncontributes noise to the mean, and no model — ours or the "
          "null — can score it.")


if __name__ == "__main__":
    main()
