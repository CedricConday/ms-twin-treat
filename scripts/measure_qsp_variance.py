"""Re-derive the QSP model's run-to-run variance, and the cohort size it forces.

Why this exists. `backtest/lomo.py` uses 128 seeds and scores the MEDIAN of
paired per-seed ratios. Those are not taste; they are forced by how noisy this
model is, and an earlier version of that file used 4 seeds and a mean and
reported a headline that was entirely its own noise. This repo's rule is that a
number in a comment comes with a way to re-derive it (see
`scripts/derive_suppressive_strength.py`), so here it is.

    PYTHONPATH=. python scripts/measure_qsp_variance.py

What it prints, and what each part is for:

1. **The spread.** Untreated damage across independent infection histories. The
   claim in `lomo.py` is a ~266-fold range with the standard error at n=4 larger
   than the mean. If that no longer holds, the cohort size in `lomo.py` is stale.

2. **Why the median.** The distribution is heavily right-tailed, so the mean is
   set by a handful of runaway histories and moves wildly with the sample. The
   mean/median ratio quantifies that.

3. **The noise floor.** Bootstrapped CV of the median at each candidate cohort
   size. This is the resolution limit on every arm-vs-arm number the screening
   stack produces — a 5pp difference between two arms means nothing if the floor
   is 10%.

4. **Paired vs unpaired.** The one design choice that buys real precision.
   Comparing two arms on the SAME seed cancels the infection history; comparing
   them on different draws does not. The ratio of the two spreads is how much
   pairing is worth, and it is the reason `simulate()` takes a seed at all.
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from bricks.profiles import PROFILES
from bricks.qsp_velez import MechanismProfile, simulate

# Cohort sizes to report a noise floor for. 128 is what lomo.py uses.
COHORT_SIZES = (4, 8, 16, 32, 48, 64, 128)

# A treated arm with a known-large effect, used for the pairing comparison.
# alpha_R is the paper's own health/autoimmunity axis and the one dial that
# moves damage by orders of magnitude, so it is the clearest signal available.
PAIRING_PROBE = MechanismProfile(label="alpha_R x1.5", alpha_R=1.5,
                                 source="variance probe, not a drug")


def _damages(profile, seeds, t_end: float) -> np.ndarray:
    """Final total damage per seed; NaN where the run left the model's regime."""
    out = []
    for seed in seeds:
        traj = simulate(profile, t_end=t_end, seed=int(seed))
        out.append(traj["total_damage"][-1] if traj["in_regime"] else np.nan)
    return np.asarray(out, dtype=float)


def _bootstrap_median_cv(values: np.ndarray, n: int, reps: int = 500,
                         rng: np.random.Generator | None = None) -> float:
    """Coefficient of variation of the median at cohort size `n`."""
    rng = rng or np.random.default_rng(0)
    medians = [np.median(rng.choice(values, size=n, replace=True)) for _ in range(reps)]
    centre = float(np.median(values))
    return float(np.std(medians) / centre) if centre else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=200, help="seeds to draw (default 200)")
    ap.add_argument("--t-end", type=float, default=730.0,
                    help="horizon in days (default 730, matching lomo.py)")
    args = ap.parse_args()

    seeds = np.arange(args.n)
    t0 = time.time()
    untreated = _damages(PROFILES["untreated"], seeds, args.t_end)
    ok = untreated[~np.isnan(untreated)]
    if ok.size < 2:
        print("too few in-regime runs to characterise; nothing to say.")
        return 1

    print(f"QSP run-to-run variance — {args.n} seeds, {args.t_end:.0f} days, "
          f"{time.time() - t0:.0f}s\n")

    print("1. SPREAD of untreated damage across infection histories")
    print(f"   in regime      {ok.size}/{untreated.size}")
    print(f"   median         {np.median(ok):>12.3f}")
    print(f"   mean           {ok.mean():>12.3f}")
    print(f"   sd             {ok.std():>12.3f}")
    print(f"   min / max      {ok.min():>12.3f} / {ok.max():.3f}")
    print(f"   fold range     {ok.max() / ok.min():>12.0f}x")
    se4 = ok.std() / np.sqrt(4)
    print(f"   se at n=4      {se4:>12.3f}   "
          f"({'LARGER' if se4 > ok.mean() else 'smaller'} than the mean)")

    print("\n2. WHY THE MEDIAN — right-tailed, so the mean is set by a few runaways")
    print(f"   mean / median  {ok.mean() / np.median(ok):>12.2f}x")
    print(f"   top decile holds "
          f"{np.sort(ok)[int(0.9 * ok.size):].sum() / ok.sum() * 100:.0f}% of total damage")

    print("\n3. NOISE FLOOR — bootstrapped CV of the median, by cohort size")
    rng = np.random.default_rng(0)
    for n in COHORT_SIZES:
        cv = _bootstrap_median_cv(ok, n, rng=rng)
        mark = "   <- lomo.py" if n == 128 else ""
        print(f"   n={n:<4} CV = {cv * 100:>5.1f}%{mark}")

    print("\n4. PAIRING — the same seed for both arms, versus independent draws")
    treated = _damages(PAIRING_PROBE, seeds, args.t_end)
    both = ~np.isnan(untreated) & ~np.isnan(treated)
    paired = treated[both] / untreated[both]
    shuffled = treated[both] / rng.permutation(untreated[both])
    print(f"   paired ratio     median {np.median(paired):>8.4f}   "
          f"IQR {np.subtract(*np.percentile(paired, [75, 25])):>8.4f}")
    print(f"   unpaired ratio   median {np.median(shuffled):>8.4f}   "
          f"IQR {np.subtract(*np.percentile(shuffled, [75, 25])):>8.4f}")
    gain = np.subtract(*np.percentile(shuffled, [75, 25])) / max(
        np.subtract(*np.percentile(paired, [75, 25])), 1e-12)
    print(f"   pairing narrows the IQR by {gain:.1f}x")

    print("\n   Read the n=128 CV as the resolution limit on every arm-vs-arm number")
    print("   the screening stack produces. A difference smaller than that is not a")
    print("   result. (bricks/qsp_velez.py is validated=False; so is this.)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
