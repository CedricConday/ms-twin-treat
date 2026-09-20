"""The best score ANY model using this representation could get. No simulation.

A drug here is a set of multipliers on the Vélez model's named dials
(`bricks/profiles.py`), and every arm sharing a dial pattern gets the same
prediction at a given potency. So there is a ceiling on the whole approach, and
it can be computed from the trial numbers alone — predict each arm from its own
dial group and see what error survives.

THREE VARIANTS, BECAUSE THE FIRST ONE FLATTERS
-----------------------------------------------
The obvious version fits each group its own mean INCLUDING the arm being
predicted, and two arms (glatiramer, daclizumab) are alone on their dials, so
they are reproduced exactly by construction and contribute zero error. Both
effects make the ceiling look better than any model could achieve:

  1. in sample, all arms          every group predicts its own mean
  2. singleton groups dropped     removes the two arms fitted exactly for free
  3. out of sample, within group  predict each arm from the OTHER arms in its
                                  group — what a model actually has to do

THE NULL PROTOCOL, STATED BECAUSE IT MOVED THE ANSWER
-------------------------------------------------------
The parallel session computed these independently and matched the model side
exactly (7.9pp and 10.9pp) while differing on the null by 0.2-0.4pp. The cause
is whether the null is also scored leave-one-out. It must be: comparing an
out-of-sample model against an in-sample null charges the model for information
the null is handed free. Both are printed so the difference stays visible
instead of being split.

Run:  PYTHONPATH=. python scripts/dial_ceiling.py
"""

from __future__ import annotations

import numpy as np

from backtest.clinical import KNOWN_OUTCOMES
from bricks.profiles import PROFILES, touched_points


def _groups() -> dict[tuple[str, ...], list[tuple[str, float]]]:
    known = {o.arm: o.relapse_change_pct for o in KNOWN_OUTCOMES
             if o.relapse_change_pct is not None and o.arm != "untreated"}
    groups: dict[tuple[str, ...], list[tuple[str, float]]] = {}
    for arm, pct in known.items():
        groups.setdefault(touched_points(PROFILES[arm]), []).append((arm, pct))
    return groups


def variant(keep_singletons: bool, out_of_sample: bool) -> dict:
    """One row of the table. `out_of_sample` scores BOTH model and null that way."""
    groups = _groups()
    arms = [(a, p, g) for g, members in groups.items() for a, p in members
            if keep_singletons or len(members) > 1]
    errs, nulls_loo, nulls_in = [], [], []
    overall = float(np.mean([p for _, p, _ in arms]))

    for arm, pct, key in arms:
        others = [q for b, q in groups[key] if b != arm]
        if out_of_sample and not others:
            continue
        pred = float(np.mean(others if out_of_sample
                             else [q for _, q in groups[key]]))
        errs.append(abs(pred - pct))
        rest = [q for b, q, _ in arms if b != arm]
        nulls_loo.append(abs(float(np.mean(rest)) - pct))
        nulls_in.append(abs(overall - pct))

    null = float(np.mean(nulls_loo if out_of_sample else nulls_in))
    mae = float(np.mean(errs))
    return {"n": len(errs), "mae": mae, "null": null, "headroom": null - mae,
            "null_in_sample": float(np.mean(nulls_in)),
            "null_leave_one_out": float(np.mean(nulls_loo))}


def dial_ceiling() -> dict:
    return {
        "in_sample": variant(True, False),
        "singletons_dropped": variant(False, False),
        "out_of_sample": variant(False, True),
    }


def main() -> int:
    r = dial_ceiling()
    print("THE DIAL-LEVEL CEILING — the best any model of this shape can score\n")
    print(f"{'variant':<30} {'MAE':>8} {'null':>8} {'headroom':>10} {'n':>4}")
    print("-" * 64)
    for label, key in (("in sample, all arms", "in_sample"),
                       ("singleton groups dropped", "singletons_dropped"),
                       ("out of sample, within group", "out_of_sample")):
        v = r[key]
        print(f"{label:<30} {v['mae']:>6.1f}pp {v['null']:>6.1f}pp "
              f"{v['headroom']:>8.1f}pp {v['n']:>4}")
    print("-" * 64)

    oos = r["out_of_sample"]
    print(f"\n  THE REAL PRIZE IS {oos['headroom']:.1f}pp.")
    print("  Read it in both directions. The representation is NOT incapable — a")
    print("  perfect dial-level model still beats the null. But out of sample,")
    print(f"  scored the way every other scorer here is scored, it wins by "
          f"{oos['headroom']:.1f}pp,")
    print("  with no simulation or fitting error charged against it. The measured")
    print("  LOMO is 45.9pp.")
    print("\n  So a replacement model does not need to be better. It needs to land")
    print(f"  within about {oos['headroom']:.0f}pp of PERFECT to clear a "
          "predict-the-mean null here.")
    print("\n  AND THE HEADROOM IS A PROPERTY OF THE ARMS, NOT THE MODEL. Twelve")
    print("  quantified arms in three multi-member dial groups is a thin exam.")
    print("  More arms PER DIAL is what the out-of-sample variant is starved of,")
    print("  and it is far cheaper than porting a richer model.")
    print(f"\n  Null protocol: out-of-sample rows score the null leave-one-out too "
          f"({oos['null_leave_one_out']:.1f}pp);")
    print(f"  the in-sample null on the same arms is {oos['null_in_sample']:.1f}pp. "
          "Comparing an")
    print("  out-of-sample model against an in-sample null would charge the model")
    print("  for information the null gets free.")
    print("\n  Trial numbers are real and cited (docs/TRIAL_ANCHORS.md). Nothing here")
    print("  is evidence about multiple sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
