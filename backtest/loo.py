"""Leave-one-arm-out: can the rule call an arm it was not fitted on?

BUILD_PLAN §8 blocker (3) asks for this, and it is the only test in the repo that
asks the question that matters -- not "does the setup reproduce the arms that
informed it", but "fit the strength on N-1 arms, predict the Nth."

How it works. The pipeline is deterministic given a cohort and an intervention's
parameters, so `backtest.response_curve` tabulates simulated relapse change as a
function of `treat` once. Fitting is then a search over that curve rather than
thousands of pipeline runs:

  for each held-out arm i:
      k_i  = the single class strength that best explains the OTHER arms
      pred = what the simulation predicts for arm i at k_i
      err  = |pred - arm i's real trial number|

Two properties of this test, stated up front because they are the result:

1. **The null is predict-the-mean.** A model that ignores which drug it is and
   answers with the mean of the training arms' outcomes. `docs/QUALITY.md`
   requires a null beside every benchmark; a rule that cannot beat this one
   carries no drug-specific information, whatever its mechanism story.

2. **Active-comparator arms are unpredictable in principle here.** Ocrelizumab,
   alemtuzumab and ponesimod are scored against another simulated drug arm, and
   under a class-level strength both arms take the same `treat`, so the predicted
   difference is exactly 0 for every possible k. They are reported separately
   rather than folded into the headline, because their error measures the rule's
   shape, not its calibration.

Run:  PYTHONPATH=. python -m backtest.loo
"""

from __future__ import annotations

import numpy as np

from backtest.clinical import KNOWN_OUTCOMES
from backtest.response_curve import TREAT_GRID, load, predict
from bricks.grounding import SUPPRESSIVE
from bricks.intervention import LIBRARY


def _quantified_arms() -> list:
    """Arms with a real number attached, which is what can be fitted or scored."""
    return [o for o in KNOWN_OUTCOMES
            if o.relapse_change_pct is not None and o.arm != "untreated"]


def _predicted_change(outcome, treat: float, table: dict) -> float:
    """What the simulation says this arm's trial would have measured, at `treat`.

    Placebo-controlled arms read straight off the curve. Active-comparator arms
    are the ratio of two points on it, because that is what their trial measured.
    """
    drug = predict(treat, 0.0, table)
    if outcome.comparator == "untreated":
        return drug
    comp_arm = LIBRARY[outcome.comparator]
    comp_treat = treat if comp_arm.mechanism == SUPPRESSIVE else comp_arm.treat
    comp = predict(comp_treat, 0.0, table)
    # At treat=1.0 the simulation removes relapses entirely, so a ratio against the
    # comparator arm divides by zero. That is a real boundary of the model, not a
    # numerical accident: with no relapses left in the comparator there is no
    # percentage difference to report. Undefined, and excluded rather than clamped.
    if abs(1 + comp / 100) < 1e-9:
        return float("nan")
    return ((1 + drug / 100) / (1 + comp / 100) - 1) * 100


def _fit_strength(training: list, table: dict) -> float:
    """The one class strength that best explains the training arms (least squares)."""
    best_k, best_sse = TREAT_GRID[0], float("inf")
    for k in np.arange(0.0, 1.0001, 0.01):
        preds = [_predicted_change(o, float(k), table) for o in training]
        if any(np.isnan(p) for p in preds):
            continue  # a strength where some arm's comparison is undefined is not comparable
        sse = sum((p - o.relapse_change_pct) ** 2 for p, o in zip(preds, training, strict=True))
        if sse < best_sse:
            best_k, best_sse = float(k), sse
    return best_k


def run_loo(table: dict | None = None) -> dict:
    table = table or load()
    arms = _quantified_arms()

    rows = []
    for held_out in arms:
        training = [o for o in arms if o.arm != held_out.arm]
        k = _fit_strength(training, table)
        pred = _predicted_change(held_out, k, table)
        null = float(np.mean([o.relapse_change_pct for o in training]))
        rows.append({
            "arm": held_out.arm,
            "comparator": held_out.comparator,
            "known": held_out.relapse_change_pct,
            "fitted_treat": k,
            "predicted": pred,
            "error": abs(pred - held_out.relapse_change_pct),
            "null_predicted": null,
            "null_error": abs(null - held_out.relapse_change_pct),
        })

    placebo = [r for r in rows if r["comparator"] == "untreated"]
    active = [r for r in rows if r["comparator"] != "untreated"]
    return {
        "rows": rows,
        "mae": float(np.mean([r["error"] for r in rows])),
        "null_mae": float(np.mean([r["null_error"] for r in rows])),
        "mae_placebo_only": float(np.mean([r["error"] for r in placebo])) if placebo else float("nan"),
        "null_mae_placebo_only": float(np.mean([r["null_error"] for r in placebo])) if placebo else float("nan"),
        "n_active_comparator": len(active),
    }


def main() -> int:
    r = run_loo()
    print("LEAVE-ONE-ARM-OUT — fit the class strength on N-1 arms, predict the Nth\n")
    print(f"{'held-out arm':<20} {'vs':<15} {'known':>8} {'fitted k':>9} {'predicted':>10} "
          f"{'error':>7} {'null err':>9}")
    print("-" * 84)
    for row in r["rows"]:
        print(f"{row['arm']:<20} {row['comparator']:<15} {row['known']:>+7.1f}% "
              f"{row['fitted_treat']:>9.2f} {row['predicted']:>+9.1f}% {row['error']:>6.1f} "
              f"{row['null_error']:>8.1f}")
    print("-" * 84)
    print(f"out-of-sample MAE: {r['mae']:.1f}pp   predict-the-mean null: {r['null_mae']:.1f}pp")
    print(f"  placebo-controlled arms only: {r['mae_placebo_only']:.1f}pp "
          f"vs null {r['null_mae_placebo_only']:.1f}pp")
    print(f"  ({r['n_active_comparator']} active-comparator arms predict exactly 0% at every k,")
    print("   because one strength per class gives both arms of those trials the same value)")
    verdict = "BEATS" if r["mae"] < r["null_mae"] else "DOES NOT BEAT"
    print(f"\nThe rule {verdict} the null.")
    if r["mae"] >= r["null_mae"]:
        print("  A single class strength carries no drug-specific information: fitted on other")
        print("  arms, it answers with roughly their average whatever drug it is asked about.")
        print("  Per-drug potency from independent data (BUILD_PLAN §8 blocker 4) is what this")
        print("  test exists to grade, and it has not landed yet.")
    print("\n  Trial numbers are real and cited (docs/TRIAL_ANCHORS.md). Simulation numbers are")
    print("  proxies from toy models. Nothing here is evidence about multiple sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
