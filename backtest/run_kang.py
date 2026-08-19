"""Run the validated harness against REAL biology — Kang IFN-beta PBMCs.

This is the milestone: the same ruler that passed selftest.py, now pointed at a
real MS-relevant perturbation. It reports the two nulls on real cell types:
  - identity (no change)      -> should score ~0: IFN-beta changes a lot.
  - global mean shift         -> should score high: the interferon signature is
                                 largely shared across immune cell types. This is
                                 the REAL bar a cell-type-specific model must beat.

No model beats a null here yet — that is honest. This run establishes the true
baseline numbers on real data. The next brick (a cell-state model, e.g. scGPT)
has to clear the mean-shift bar on held-out cell types to earn any claim.

Run:  python -m backtest.run_kang        (repo root, venv active)
"""

from __future__ import annotations

from data.kang import to_benchmark
from bricks.baselines import GlobalMeanShiftNull, IdentityNull


def main() -> int:
    print("Loading Kang 2018 (IFN-beta PBMCs) and building the backtest...")
    bench = to_benchmark()
    print(f"  cell types scored: {len(bench.cell_types)}")
    for ct in bench.cell_types:
        print(f"    - {ct}")

    identity = IdentityNull()
    mean_shift = GlobalMeanShiftNull(bench.global_mean_delta())

    print("\nPer-cell-type delta_pearson (real IFN-beta response):")
    id_tbl = bench.evaluate(identity)[["delta_pearson", "deg_precision_at_k"]]
    ms_tbl = bench.evaluate(mean_shift)[["delta_pearson", "deg_precision_at_k"]]
    joined = id_tbl.join(ms_tbl, lsuffix="_identity", rsuffix="_meanshift")
    print(joined.round(3).to_string())

    print("\nAggregate (mean across cell types):")
    print(f"  identity   delta_pearson = {id_tbl['delta_pearson'].mean(): .4f}")
    print(f"  mean-shift delta_pearson = {ms_tbl['delta_pearson'].mean(): .4f}   <-- the bar to beat")
    print("\nThe ruler works on real biology. Now build a brick that beats the bar.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
