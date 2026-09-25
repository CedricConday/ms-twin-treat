"""Rebuild the response tables at 1024 seeds, into their own files.

Gap G2 in docs/GATE_GAP_ANALYSIS.md: a perfect model has about 4.5pp to win
on the point-ARR exam and one prediction carries about 6pp of stochastic
noise at 128 seeds (13.7% CV of the median, scripts/measure_qsp_variance.py).
Bootstrapped CV fell 21.2% -> 13.7% from n=48 to n=128; at n=1024 it should
sit near 5%. This builds the same tables with the same cells and the same
median statistic, eight times the cohort, and writes them beside the
existing ones:

    results/mechanism_curve_n1024.json     the five quantified patterns
    results/exam_v2_curve_n1024.json       the three patterns exam v2 adds

The 128-seed files are never touched. Both scorers read the 128-seed tables
by default; pass the n1024 paths explicitly to compare.

Cost: ~8x the 128-seed build (~8 min), so about an hour per file on this box.
Run detached:

    PYTHONPATH=. nohup python scripts/build_curve_n1024.py > results/build_n1024.log 2>&1 &
"""

from __future__ import annotations

import time
from pathlib import Path

import backtest.exam_v2 as exam_v2
import backtest.lomo as lomo

N = 1024
RESULTS = Path(__file__).resolve().parent.parent / "results"


def main() -> int:
    lomo.SEEDS = tuple(range(N))
    lomo.CACHE = RESULTS / f"mechanism_curve_n{N}.json"
    exam_v2.SEEDS = lomo.SEEDS
    exam_v2.CURVE = RESULTS / f"exam_v2_curve_n{N}.json"

    t0 = time.time()
    print(f"building {lomo.CACHE.name} at {N} seeds ...", flush=True)
    lomo.build_table()
    print(f"  done in {(time.time() - t0) / 60:.1f} min", flush=True)

    t1 = time.time()
    print(f"building {exam_v2.CURVE.name} at {N} seeds ...", flush=True)
    # load_curve reads lomo.load(), which now points at the n1024 base file,
    # and builds only the patterns that file lacks.
    exam_v2.load_curve()
    print(f"  done in {(time.time() - t1) / 60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
