"""Run the validated harness against REAL biology — Kang IFN-beta PBMCs.

This is the milestone: the same ruler that passed selftest.py, now pointed at a
real MS-relevant perturbation. It reports the nulls on real cell types:
  - identity (no change)      -> should score ~0: IFN-beta changes a lot.
  - leave-one-out mean shift  -> THE CANONICAL BAR. Predicts a cell type from the
                                 average of the OTHERS, which is exactly the
                                 information a leave-one-cell-type-out model has.
  - global mean shift         -> the same idea but averaged over ALL cell types
                                 including the held-out one, so it is leaky in
                                 the null's favour. Retained as a HARDER
                                 secondary check (see bricks/baselines.py).

It also reports each fold's CEILING — the best delta_pearson anything could
score against a measurement that carries sampling noise. A score without its
ceiling is not interpretable: on Kang, one fold (Megakaryocytes, 63/69 cells)
has essentially no measurable signal, and it silently drags every aggregate.

bricks/cell_scgpt.py clears both bars (0.8732). Grade new bricks the same way:
canonical bar first, leaky bar as the stronger claim, ceiling always in view.

Run:  python -m backtest.run_kang        (repo root, venv active)
"""

from __future__ import annotations

from data.kang import to_benchmark
from bricks.baselines import (GlobalMeanShiftNull, IdentityNull,
                              LeaveOneOutMeanShiftNull)


def main() -> int:
    print("Loading Kang 2018 (IFN-beta PBMCs) and building the backtest...")
    bench = to_benchmark()
    print(f"  cell types scored: {len(bench.cell_types)}")
    for ct in bench.cell_types:
        print(f"    - {ct}")

    identity = IdentityNull()
    canonical = LeaveOneOutMeanShiftNull(bench)
    leaky = GlobalMeanShiftNull(bench.global_mean_delta())

    id_tbl = bench.evaluate(identity)
    loo_tbl = bench.evaluate(canonical)
    ms_tbl = bench.evaluate(leaky)

    print("\nPer-cell-type delta_pearson (real IFN-beta response):")
    header = f"  {'cell type':<20}{'identity':>10}{'LOO-shift':>11}{'global':>9}"
    if bench.reliability:
        header += f"{'ceiling':>9}"
    print(header)
    for ct in bench.cell_types:
        row = (f"  {ct:<20}{id_tbl.loc[ct, 'delta_pearson']:10.3f}"
               f"{loo_tbl.loc[ct, 'delta_pearson']:11.3f}"
               f"{ms_tbl.loc[ct, 'delta_pearson']:9.3f}")
        if bench.reliability:
            ceil = bench.ceiling(ct)
            flag = "  <-- noise fold" if ct in bench.unreliable_cell_types else ""
            row += f"{ceil:9.3f}{flag}"
        print(row)

    print("\nAggregate (mean across cell types):")
    print(f"  identity                       {id_tbl['delta_pearson'].mean(): .4f}")
    print(f"  leave-one-out mean shift       {loo_tbl['delta_pearson'].mean(): .4f}"
          f"   <-- THE CANONICAL BAR")
    print(f"  global mean shift (leaky)      {ms_tbl['delta_pearson'].mean(): .4f}"
          f"   <-- harder secondary; it sees the held-out delta")
    if bench.reliability:
        ceilings = [bench.ceiling(ct) for ct in bench.cell_types]
        print(f"  ceiling (a PERFECT model)      {sum(ceilings)/len(ceilings): .4f}")
        noisy = bench.unreliable_cell_types
        if noisy:
            keep = [ct for ct in bench.cell_types if ct not in noisy]
            print(f"\n  {len(noisy)} fold(s) cannot be scored at all: {noisy}")
            print(f"  Restricted to the {len(keep)} measurable cell types:")
            print(f"    LOO-shift bar {loo_tbl.loc[keep, 'delta_pearson'].mean():.4f}"
                  f"   global {ms_tbl.loc[keep, 'delta_pearson'].mean():.4f}"
                  f"   ceiling {sum(bench.ceiling(c) for c in keep)/len(keep):.4f}")
    print("\nThe ruler works on real biology, and now reports what is reachable.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
