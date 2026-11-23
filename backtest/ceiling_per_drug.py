"""How well could the Pernice port score if every DRUG had its own magnitude?

Not an exam: the magnitudes here are chosen with the outcome in view, so the result is a
ceiling on what any independent magnitude source could buy, and a list of the arms no
magnitude can rescue. Two things differ from exam v2's table:

D1  The comparator keeps its OWN potency (the one that reproduces its placebo-controlled
    arm) instead of sharing the tested drug's. Exam v2 scores "drug at s against
    comparator at s", which caps how far an active-comparator arm can fall.
D2  Dials that are raised (the two death rates) are swept as a steady-state reduction of
    the cell count, a = 1 - 1/multiplier, up to 0.95 (multiplier 20). Exam v2 stops at
    multiplier 1.95, i.e. at most a 49 % reduction, which no depleting antibody is.

Usage: python -m backtest.ceiling_per_drug      writes results/ceiling_per_drug.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from backtest import exam_v2 as E
from bricks import profiles_pernice, qsp_pernice
from bricks.qsp_pernice import PerniceProfile
from bricks.sormani import predict_relapse_ratio

OUT = Path(__file__).resolve().parent.parent / "results" / "ceiling_per_drug.json"
REDUCTIONS = (0.0, 0.2, 0.33, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
FINE = np.linspace(0.0, 0.95, 381)


def _run(profile: PerniceProfile) -> float | None:
    tr = qsp_pernice.simulate(profile, sample_hours=24.0, readings=qsp_pernice.READING_VARIANTS["s2"])
    return tr["lesion_load"] if tr["in_regime"] and tr["lesion_load"] > 0 else None


def depletion_curve(arm: str, untreated_load: float) -> tuple[list[float], list[float]]:
    """Percent relapse change against untreated, by steady-state count reduction."""
    base = profiles_pernice.PROFILES[arm]
    xs, ys = [], []
    for a in REDUCTIONS:
        mult = {k: 1.0 / (1.0 - a) for k in base.multipliers}
        load = untreated_load if a == 0 else _run(PerniceProfile(label=f"{arm}@a{a:g}", multipliers=mult))
        if load is not None:
            xs.append(a)
            ys.append(predict_relapse_ratio(load / untreated_load).percent_change)
    return xs, ys


def main() -> int:
    E.select_model("pernice", "s2")
    curve = E.load_curve(verbose=False)
    untreated_load = _run(qsp_pernice.UNTREATED_PROFILE)
    raised = {a for a, p in profiles_pernice.PROFILES.items()
              if p.multipliers and all(v > 1.0 for v in p.multipliers.values()) and p.dac_dose == 0
              and set(p.multipliers) <= {"Teff_death", "Treg_death"}}
    dep = {a: depletion_curve(a, untreated_load) for a in sorted(raised)}

    def absolute(arm: str, x: float) -> float:
        """The arm against untreated at magnitude x (potency, or count reduction for D2 arms)."""
        if arm in dep:
            return float(np.interp(x, *dep[arm]))
        return E._curve_value(E._pattern(arm), x, curve)

    outcomes = {o.arm: o for o in E.treated()}
    own, rows = {}, []
    # pass 1: placebo-controlled arms fix each drug's own magnitude
    for arm, o in outcomes.items():
        if o.comparator != "untreated":
            continue
        iv = E.interval(o)
        d = [E.dist(absolute(arm, x), iv) for x in FINE]
        i = int(np.nanargmin(d))
        own[arm] = float(FINE[i])
        rows.append({"arm": arm, "comparator": "untreated", "interval": list(iv), "magnitude": own[arm],
                     "predicted": absolute(arm, own[arm]), "error": float(d[i])})
    # pass 2: active-comparator arms, comparator at its own magnitude
    for arm, o in outcomes.items():
        if o.comparator == "untreated":
            continue
        iv = E.interval(o)
        comp = absolute(o.comparator, own[o.comparator]) if o.comparator in own else float("nan")
        preds = [((1 + absolute(arm, x) / 100) / (1 + comp / 100) - 1) * 100 for x in FINE]
        d = [E.dist(p, iv) for p in preds]
        i = int(np.nanargmin(d))
        rows.append({"arm": arm, "comparator": o.comparator, "comparator_predicted": comp, "interval": list(iv),
                     "magnitude": float(FINE[i]), "predicted": float(preds[i]), "error": float(d[i])})
    v2 = E.run_interval_lomo(curve)
    null = {r["arm"]: r["null_error"] for r in v2["rows"]}
    for r in rows:
        r["null_error"] = null[r["arm"]]
        r["count_reduction_dial"] = r["arm"] in dep
    mae = float(np.mean([r["error"] for r in rows]))
    out = {"ceiling_mae": mae, "exam_v2_s1": v2["mae"], "null_mae": v2["null_mae"],
           "green_needs": 0.8 * v2["null_mae"], "n_arms": len(rows),
           "unreachable": [r["arm"] for r in rows if r["error"] > 5.0],
           "depletion_curves": {a: {"reduction": xs, "percent_change": ys} for a, (xs, ys) in dep.items()},
           "rows": rows}
    OUT.write_text(json.dumps(out, indent=2, default=lambda v: None if not np.isfinite(v) else v))
    fmt = lambda v: "inf" if not np.isfinite(v) else f"{v:.0f}"  # noqa: E731
    for r in sorted(rows, key=lambda r: -r["error"]):
        print(f"  {r['arm']:<18} vs {r['comparator'][:14]:<14} truth ({fmt(r['interval'][0])},{fmt(r['interval'][1])})"
              f"  best {r['predicted']:6.1f} at {r['magnitude']:.2f}  error {r['error']:5.1f}  null {r['null_error']:5.1f}")
    print(f"\nceiling MAE {mae:.1f}pp | exam v2 S1 {v2['mae']:.1f} | null {v2['null_mae']:.1f} | green needs {0.8 * v2['null_mae']:.1f}")
    print("no magnitude rescues:", ", ".join(out["unreachable"]) or "none")
    for a, (xs, ys) in dep.items():
        print(f"  {a} by count reduction: " + ", ".join(f"{x:g}:{y:.0f}" for x, y in zip(xs, ys)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
