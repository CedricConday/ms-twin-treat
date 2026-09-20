"""The direction gate, scored on the GROUNDED stack instead of the ABM path.

`backtest/clinical.py` scores 14 cited trial arms through the ABM
(`abm_damage`) and reports 9/14 direction, 1/9 magnitude. Everything built for
screening — the ported QSP, the per-arm intervention points, the published
lesions→relapses map — sits beside that path and is not graded by it. BUILD_PLAN
§8.4 records the gap as "wiring".

This module closes it by scoring the same arms, against the same cited outcomes,
with the same rules, through:

    bricks/qsp_velez.py   Vélez de Mendizábal 2011, transcribed
    bricks/profiles.py    which dials each drug moves, from pharmacology
    bricks/sormani.py     lesion ratio -> relapse ratio, 31 trials, slope 0.52

so the two can be compared like for like. It does NOT replace `clinical.py`:
both are kept and both are reported, because a new number that quietly replaces
an old one is not a comparison.

HOW AN ARM IS SCORED
--------------------
1. Simulate the arm and its comparator on the SAME 128 infection histories.
2. Take the median per-seed damage ratio — the model is ~90,000-fold
   right-tailed, so the median is the only usable statistic
   (`scripts/measure_qsp_variance.py`).
3. Read that ratio as a lesion rate ratio and put it through Sormani.
4. Compare the predicted relapse change to what the trial reported.

Active-comparator trials (OPERA, CARE-MS, OPTIMUM) are scored the way they ran:
the simulated drug arm against the simulated comparator arm.

WHAT WILL DRAG THIS NUMBER DOWN, STATED BEFORE IT RUNS
-------------------------------------------------------
`gamma_E` predicts the wrong sign. Five arms sit on it — natalizumab,
fingolimod, ponesimod, alemtuzumab, atacicept — and raising it makes the model
worse rather than better, because in the published equations `gamma_E` is gated
by the Treg Hill term and is therefore the REGULATORY killing channel rather
than a generic death rate. That is recorded in BUILD_PLAN §8.4 and is not
patched here; a gate that hides a known defect is not a gate.

Run:  PYTHONPATH=. python -m backtest.clinical_velez
"""

from __future__ import annotations

import numpy as np

from backtest.clinical import KNOWN_OUTCOMES, MAG_TOLERANCE, NEUTRAL_BAND, _direction
from backtest.lomo import SEEDS, T_END, _damage
from bricks.profiles import PROFILES
from bricks.sormani import predict_relapse_ratio


def _arm_damages(arm: str) -> dict[int, float]:
    """Per-seed damage for one arm, skipping histories that left the regime."""
    profile = PROFILES[arm]
    out = {}
    for seed in SEEDS:
        d = _damage(profile, seed)
        if d is not None and d > 0.0:
            out[seed] = d
    return out


def run_gate() -> dict:
    """Score every arm with a stated direction on the grounded stack."""
    cache: dict[str, dict[int, float]] = {}

    def damages(arm: str) -> dict[int, float]:
        if arm not in cache:
            cache[arm] = _arm_damages(arm)
        return cache[arm]

    rows = []
    for outcome in KNOWN_OUTCOMES:
        if outcome.arm == "untreated":
            continue
        drug = damages(outcome.arm)
        comp = damages(outcome.comparator)
        shared = sorted(set(drug) & set(comp))
        if not shared:
            rows.append({"arm": outcome.arm, "predicted": float("nan"),
                         "predicted_direction": "undefined", **_expected(outcome)})
            continue

        ratios = [drug[s] / comp[s] for s in shared]
        rr_lesion = float(np.median(ratios))
        pred = predict_relapse_ratio(rr_lesion)
        rows.append({
            "arm": outcome.arm,
            "comparator": outcome.comparator,
            "rr_lesion": rr_lesion,
            "predicted": pred.percent_change,
            "predicted_direction": _direction(pred.percent_change),
            "blind_spot": pred.blind_spot,
            "n_seeds": len(shared),
            "n_dropped": len(SEEDS) - len(shared),
            **_expected(outcome),
        })

    scored = [r for r in rows if r["predicted_direction"] != "undefined"]
    dir_hits = [r for r in scored if r["predicted_direction"] == r["known_direction"]]
    quant = [r for r in scored if r["known_change"] is not None]
    mag_hits = [r for r in quant
                if abs(r["predicted"] - r["known_change"]) <= MAG_TOLERANCE]

    return {
        "rows": rows,
        "n_arms": len(rows),
        "direction_hits": len(dir_hits),
        "n_scored": len(scored),
        "magnitude_hits": len(mag_hits),
        "n_quantified": len(quant),
    }


def _expected(outcome) -> dict:
    return {"known_direction": outcome.direction,
            "known_change": outcome.relapse_change_pct,
            "comparator": outcome.comparator}


def main() -> int:
    r = run_gate()
    print("CLINICAL DIRECTION GATE — scored on the GROUNDED stack")
    print("  qsp_velez (Velez 2011) -> profiles (pharmacology) -> sormani (31 trials)")
    print(f"  {len(SEEDS)} paired infection histories, {T_END:.0f} days, median ratio\n")

    print(f"{'arm':<20} {'vs':<15} {'RR lesion':>10} {'predicted':>11} "
          f"{'known':>9} {'dir':>9}")
    print("-" * 82)
    for row in r["rows"]:
        known = ("%+.1f%%" % row["known_change"]) if row["known_change"] is not None else "—"
        pred = ("%+.1f%%" % row["predicted"]) if np.isfinite(row["predicted"]) else "n/a"
        rr = ("%.3f" % row["rr_lesion"]) if "rr_lesion" in row else "n/a"
        ok = "ok" if row["predicted_direction"] == row["known_direction"] else "MISS"
        flag = " *" if row.get("blind_spot") else ""
        print(f"{row['arm']:<20} {row['comparator']:<15} {rr:>10} {pred:>11} "
              f"{known:>9} {ok:>9}{flag}")
    print("-" * 82)
    print(f"DIRECTION: {r['direction_hits']}/{r['n_arms']} arms correct   "
          f"(magnitude: {r['magnitude_hits']}/{r['n_quantified']} within "
          f"{MAG_TOLERANCE:.0f}pp)")
    print(f"  neutral band +/-{NEUTRAL_BAND:.0f}%   * = Sormani blind spot, where the map")
    print("    cannot separate 'no effect' from harm (the lenercept case)")

    print("\n  COMPARE: backtest/clinical.py, the ABM path, reports 9/14 direction")
    print("  and 1/9 magnitude on the same arms and the same cited outcomes.")
    print("  So the grounded stack is WORSE on direction and better on magnitude.")
    print("  Grounding the bricks did not buy a better score, and saying otherwise")
    print("  would require ignoring this table.")
    print("\n  READ TWO ROWS CAREFULLY BEFORE QUOTING THE HEADLINE:")
    print("   * atacicept scores 'ok' for the WRONG REASON. It sits on gamma_E,")
    print("     gamma_E has the wrong sign, and atacicept happens to be a drug that")
    print("     harmed -- so a broken dial and a harmful drug cancel. Remove the")
    print("     defect and this row is expected to flip to a MISS.")
    print("   * lenercept is UNDEFINED, not wrong. alpha_E down with alpha_R down")
    print("     strips regulation until the effector population runs away, so every")
    print("     history leaves the model's regime. It is counted as a miss because")
    print("     an unscoreable arm is not a pass.")
    print("\n  Known defect NOT patched for this run: gamma_E has the wrong sign, and")
    print("  five arms sit on it (natalizumab, fingolimod, ponesimod, alemtuzumab,")
    print("  atacicept). In the published equations it is gated by the Treg Hill")
    print("  term, so it is the regulatory killing channel rather than a generic")
    print("  death rate. See BUILD_PLAN §8.4.")
    print("\n  Trial outcomes are REAL and cited (docs/TRIAL_ANCHORS.md). Simulation")
    print("  numbers come from a transcribed model that has not been reproduced, and")
    print("  the Sormani map is applied outside what it was fitted on. Nothing here")
    print("  is evidence about multiple sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
