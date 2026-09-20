"""LOMO under the capacity extension: does the sign fix survive being held out?

`backtest/lomo.py` is the gate the screen is blocked on. On the pure Vélez
transcription it reports 45.9pp out-of-sample MAE against a 12.3pp
predict-the-mean null — the model is worse than answering with the training
arms' average, so `screen.rank_candidates` refuses to rank.

BUILD_PLAN §8.4 records why, and it is structural: damage goes as `(E/a)^n`, so
it is set by peak excursions, and effectors recruit their own regulators. Remove
effectors and the brake releases into a larger excursion. Five of the nine
quantified arms sit on that dial. The one structural change that flips the sign
is a carrying capacity on effector proliferation that BINDS AT THE OPERATING
POINT (`simulate(carrying_capacity=K)`, off by default, `bricks/qsp_velez.py`).

That flip was measured on the gamma_E dial in isolation. This file asks the only
question that matters for a screen: **does it move the out-of-sample number?**

THE FIT IS HONEST ABOUT K, WHICH IS THE WHOLE POINT
-----------------------------------------------------
K is a second free parameter, and choosing it by looking at the headline would
be fitting the test set — the exact failure the LOMO exists to prevent. So K is
fitted the same way the potency is: **inside each fold, on the training arms
only**. The held-out mechanism sees neither. Reported alongside is the per-K
table with K fixed globally, which is NOT out-of-sample in K and is labelled as
a diagnostic, never as the headline.

WHAT A PASS WOULD AND WOULD NOT LICENCE
-----------------------------------------
Inverting this unblocks ranking within the mechanisms the model can represent.
It does not repair the extension's magnitude: at K=2000, three-fold depletion
buys 14% less damage against trials reporting 55-68% relapse reductions. And the
capacity is an extension, not transcription — anything scored through it is a
different model from the published one, and must say so.

THE ANSWER, MEASURED 2026-09-20
---------------------------------
**It does not lift the gate.** 45.6pp against the same 12.3pp null, versus
45.9pp uncapped — inside the noise of no change. Every fixed-K diagnostic loses
too; the best, K=1000, reaches 43.4pp, still 3.5x the null.

What the capacity actually bought is visible in the gamma_E fold: the sign is no
longer wrong, and the prediction is -3.6% where the five trials report -30% to
-68%. So the defect was never only the sign. Removing effectors in this model
cannot produce a large damage reduction, because damage is set by peak
excursions and the capacity only damps them. Right direction, wrong order of
magnitude, and the fold error barely moves (53.2pp -> 49.6pp).

The daclizumab fold is worse under the extension (108.8pp): `alpha_R` down with
a binding cap predicts +63.8% where DECIDE reports -45%. One arm, so read it
beside the others, not alone.

**Conclusion: the screen stays kill-only.** This closes the cheapest remaining
route to ranking, and it closes it with a measurement rather than an argument.

Run:  PYTHONPATH=. python -m backtest.lomo_capacity        (builds caches on first run)
"""

from __future__ import annotations

import json

import numpy as np

from backtest.clinical import KNOWN_OUTCOMES
from backtest.lomo import RESULTS, _predicted, load, mechanism_groups, run_lomo

# Capacities probed. None is the published model (no cap). The rest bracket the
# sweep in bricks/qsp_velez.py, where the sign flips between K=3000 and K=2000
# against a baseline effector level of ~1000: above that the cap only clips
# excursions, below it the population is genuinely resource-limited.
K_GRID: tuple[float | None, ...] = (None, 10000.0, 3000.0, 2000.0, 1500.0, 1000.0)

# Same potency resolution as backtest/lomo.py's own fit, so the two headlines
# differ only by the model, not by the search.
S_GRID = [round(float(x), 2) for x in np.arange(0.0, 0.951, 0.01)]

OUT = RESULTS / "lomo_capacity.json"


def _tables(rebuild: bool = False) -> dict[float | None, dict]:
    """One response table per capacity. Cached per K by `lomo.cache_path`."""
    return {k: load(rebuild=rebuild, carrying_capacity=k) for k in K_GRID}


def _known() -> dict[str, float]:
    return {o.arm: o.relapse_change_pct for o in KNOWN_OUTCOMES
            if o.relapse_change_pct is not None and o.arm != "untreated"}


def _fit(training_arms: list[str], known: dict[str, float],
         tables: dict[float | None, dict]) -> tuple[float | None, float, float]:
    """Best (capacity, potency) on the TRAINING arms only. Returns (K, s, sse)."""
    best: tuple[float | None, float, float] = (None, 0.0, float("inf"))
    for k, table in tables.items():
        for s in S_GRID:
            sse, ok = 0.0, True
            for arm in training_arms:
                p = _predicted(arm, s, table)
                if np.isnan(p):
                    ok = False
                    break
                sse += (p - known[arm]) ** 2
            if ok and sse < best[2]:
                best = (k, s, sse)
    return best


def run(tables: dict[float | None, dict] | None = None) -> dict:
    """Fold over mechanism groups, fitting (K, s) blind to the held-out group."""
    tables = tables or _tables()
    groups = mechanism_groups()
    known = _known()

    rows, folds = [], []
    for held_pattern, held_arms in groups.items():
        training_arms = [a for pat, arms in groups.items() if pat != held_pattern
                         for a in arms]
        k, s, _ = _fit(training_arms, known, tables)
        null = float(np.mean([known[a] for a in training_arms]))

        for arm in held_arms:
            pred = _predicted(arm, s, tables[k])
            rows.append({
                "held_mechanism": "|".join(held_pattern),
                "arm": arm,
                "known": known[arm],
                "fitted_capacity": k,
                "fitted_potency": s,
                "predicted": pred,
                "error": abs(pred - known[arm]),
                "null_predicted": null,
                "null_error": abs(null - known[arm]),
            })
        folds.append({
            "mechanism": "|".join(held_pattern),
            "n_arms": len(held_arms),
            "fitted_capacity": k,
            "fitted_potency": s,
            "mae": float(np.mean([r["error"] for r in rows if
                                  r["held_mechanism"] == "|".join(held_pattern)])),
            "null_mae": float(np.mean([abs(null - known[a]) for a in held_arms])),
        })

    fixed = {}
    for k, table in tables.items():
        r = run_lomo(table=table)
        fixed[str(k)] = {"mae": r["mae"], "null_mae": r["null_mae"]}

    return {
        "rows": rows,
        "folds": folds,
        "fixed_k_diagnostic": fixed,
        "mae": float(np.mean([r["error"] for r in rows])),
        "null_mae": float(np.mean([r["null_error"] for r in rows])),
        "note": ("Capacity is an EXTENSION to Velez de Mendizabal 2011, not "
                 "transcription. K and potency are fitted per fold on training "
                 "arms only; the held-out mechanism sees neither. validated=False."),
    }


def main() -> int:
    r = run()
    print("LOMO UNDER THE CAPACITY EXTENSION — K fitted per fold, on training arms only\n")
    print(f"{'held mechanism':<20} {'arm':<20} {'known':>8} {'K':>8} {'potency':>8} "
          f"{'predicted':>10} {'error':>7} {'null err':>9}")
    print("-" * 97)
    for row in r["rows"]:
        k = row["fitted_capacity"]
        print(f"{row['held_mechanism']:<20} {row['arm']:<20} {row['known']:>+7.1f}% "
              f"{('none' if k is None else f'{k:g}'):>8} {row['fitted_potency']:>8.2f} "
              f"{row['predicted']:>+9.1f}% {row['error']:>6.1f} {row['null_error']:>8.1f}")
    print("-" * 97)
    for f in r["folds"]:
        k = f["fitted_capacity"]
        print(f"  {f['mechanism']:<20} {f['n_arms']} arm(s)   "
              f"K {('none' if k is None else f'{k:g}'):>6}   "
              f"MAE {f['mae']:>6.1f}pp   null {f['null_mae']:>6.1f}pp")

    print(f"\nout-of-sample MAE: {r['mae']:.1f}pp   predict-the-mean null: {r['null_mae']:.1f}pp")
    beats = r["mae"] < r["null_mae"]
    print(f"\nThe extended model {'BEATS' if beats else 'DOES NOT BEAT'} "
          "the null on unseen mechanisms.")

    print("\nDIAGNOSTIC — K fixed globally (NOT out-of-sample in K, do not quote as the gate):")
    for k, v in r["fixed_k_diagnostic"].items():
        print(f"  K={k:<8} MAE {v['mae']:>6.1f}pp   null {v['null_mae']:>6.1f}pp")

    OUT.write_text(json.dumps(r, indent=2))
    print(f"\nwritten: {OUT}")
    print("  Capacity is an extension, not transcription. Trial numbers are real and")
    print("  cited (docs/TRIAL_ANCHORS.md); simulation numbers are proxies and the")
    print("  Sormani map is applied outside what it was fitted on. Nothing here is")
    print("  evidence about multiple sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
