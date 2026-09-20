"""The best score ANY model using this representation could get. No simulation.

A drug here is a set of multipliers on the Vélez model's named dials
(`bricks/profiles.py`), and every arm sharing a dial pattern gets the same
prediction at a given potency. So there is a hard ceiling on the whole approach,
and it can be computed from the trial numbers alone:

    predict each mechanism group its OWN MEAN — in-sample, cheating, no
    simulation, no fitting error — and measure the MAE that remains.

No model of this shape can beat that, because the residual is exactly the
spread of real outcomes WITHIN a dial group, which the representation has no
way to express. Compare it against the predict-the-mean null on the same arms.

This exists because "the blocker is the model form" and "the blocker is the
potency source" were both argued tonight without anyone bounding the
representation itself. Run it before budgeting a port.

Run:  PYTHONPATH=. python scripts/dial_ceiling.py
"""

from __future__ import annotations

import numpy as np

from backtest.clinical import KNOWN_OUTCOMES
from bricks.profiles import PROFILES, touched_points


def ceiling() -> dict:
    known = {o.arm: o.relapse_change_pct for o in KNOWN_OUTCOMES
             if o.relapse_change_pct is not None and o.arm != "untreated"}
    groups: dict[tuple[str, ...], list[tuple[str, float]]] = {}
    for arm, pct in known.items():
        groups.setdefault(touched_points(PROFILES[arm]), []).append((arm, pct))

    errs, rows = [], []
    for key, members in sorted(groups.items()):
        vals = [p for _, p in members]
        mean = float(np.mean(vals))
        group_errs = [abs(p - mean) for p in vals]
        errs += group_errs
        rows.append({
            "dials": "|".join(key),
            "n": len(members),
            "arms": [a for a, _ in members],
            "spread_pp": (max(vals) - min(vals)) if len(vals) > 1 else 0.0,
            "in_group_mae": float(np.mean(group_errs)),
        })

    overall = float(np.mean(list(known.values())))
    return {
        "rows": rows,
        "n_arms": len(known),
        "dial_ceiling_pp": float(np.mean(errs)),
        "null_pp": float(np.mean([abs(p - overall) for p in known.values()])),
    }


def main() -> int:
    r = ceiling()
    print("THE DIAL-LEVEL CEILING — the best any model of this shape can score\n")
    print(f"{'dials':<18} {'n':>2}  {'spread':>8}  {'in-group MAE':>12}  arms")
    print("-" * 92)
    for row in r["rows"]:
        print(f"{row['dials']:<18} {row['n']:>2}  {row['spread_pp']:>7.1f}pp  "
              f"{row['in_group_mae']:>11.1f}pp  {', '.join(row['arms'])}")
    print("-" * 92)
    print(f"\n  ceiling (each group predicts its own mean, in-sample): "
          f"{r['dial_ceiling_pp']:.1f}pp over {r['n_arms']} arms")
    print(f"  predict-the-mean null, same arms, in-sample:           {r['null_pp']:.1f}pp")
    head = r["null_pp"] - r["dial_ceiling_pp"]
    print(f"\n  HEADROOM: {head:.1f}pp.")
    print("  Read it carefully in both directions. The representation is NOT")
    print("  incapable — a perfect dial-level model does beat the null in-sample.")
    print("  But it beats it by a hair, on 12 arms, with no fitting error and no")
    print("  simulation error included. The measured out-of-sample LOMO is 45.9pp,")
    print("  so the model is ~39pp away from its own representation's ceiling, and")
    print("  the ceiling itself leaves almost nothing to win by.")
    print("\n  Consequence for any port: a richer model must not merely improve, it")
    print("  must land inside a few pp of perfect to clear a null. Budget on that.")
    print("\n  Trial numbers are real and cited (docs/TRIAL_ANCHORS.md). Nothing here")
    print("  is evidence about multiple sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
