"""Exam v3: per-drug magnitude from pharmacodynamic anchors (docs/EXAM_V3_PREREG.md).

    python -m backtest.exam_v3 [--anchors data/pd_anchors.tsv]

V3-P takes each drug's anchor as its magnitude, unchanged; nothing is fitted. V3-K scales
every anchor by one kappa fitted on the arms outside the held fold. Arms, intervals, folds
and null are exam v2's. The model side is the Pernice port with the comparator at its own
magnitude (D1) and raised dials swept as a count reduction (D2).
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from backtest import exam_v2 as E
from backtest.ceiling_per_drug import _run, depletion_curve
from bricks import profiles_pernice, qsp_pernice

ROOT = Path(__file__).resolve().parent.parent
ANCHORS = ROOT / "data" / "pd_anchors.tsv"
OUT = ROOT / "results" / "exam_v3.json"
KAPPA_GRID = np.round(np.arange(0.0, 2.0001, 0.02), 2)
CAP = 0.95


def load_anchors(path: Path) -> dict[str, float]:
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("status", "").strip() == "sourced" and row.get("anchor", "").strip():
                out[row["drug"].strip()] = min(CAP, max(0.0, float(row["anchor"])))
    return out


class Model:
    def __init__(self) -> None:
        E.select_model("pernice", "s2")
        self.curve = E.load_curve(verbose=False)
        load = _run(qsp_pernice.UNTREATED_PROFILE)
        self.dep = {a: depletion_curve(a, load) for a, p in profiles_pernice.PROFILES.items()
                    if p.multipliers and p.dac_dose == 0 and all(v > 1.0 for v in p.multipliers.values())
                    and set(p.multipliers) <= {"Teff_death", "Treg_death"}}

    def absolute(self, arm: str, a: float) -> float:
        """Percent relapse change of the arm against untreated at magnitude a."""
        if arm in self.dep:
            return float(np.interp(a, *self.dep[arm]))
        return E._curve_value(E._pattern(arm), a, self.curve)

    def predict(self, o, anchors: dict[str, float], kappa: float = 1.0) -> float | None:
        if o.arm not in anchors or (o.comparator != "untreated" and o.comparator not in anchors):
            return None
        drug = self.absolute(o.arm, min(CAP, kappa * anchors[o.arm]))
        if o.comparator == "untreated":
            return drug
        comp = self.absolute(o.comparator, min(CAP, kappa * anchors[o.comparator]))
        return ((1.0 + drug / 100.0) / (1.0 + comp / 100.0) - 1.0) * 100.0


def _fit_kappa(model: Model, training: list, anchors: dict[str, float]) -> float:
    best, best_sse = 1.0, float("inf")
    for k in KAPPA_GRID:
        preds = [(model.predict(o, anchors, float(k)), o) for o in training]
        preds = [(p, o) for p, o in preds if p is not None and np.isfinite(p)]
        if not preds:
            continue
        sse = sum(E.dist(p, E.interval(o)) ** 2 for p, o in preds)
        if sse < best_sse:
            best, best_sse = float(k), sse
    return best


def run(model: Model, anchors: dict[str, float], fitted: bool) -> dict:
    groups = E.all_groups()
    rows, folds = [], []
    for held, arms in groups.items():
        training = [o for p, os_ in groups.items() if p != held for o in os_]
        null = float(np.mean([E.anchor(o) for o in training]))
        kappa = _fit_kappa(model, training, anchors) if fitted else 1.0
        errs, nerrs = [], []
        for o in arms:
            iv = E.interval(o)
            pred = model.predict(o, anchors, kappa)
            anchored = bool(pred is not None and np.isfinite(pred))
            ne = E.dist(null, iv)
            e = E.dist(pred, iv) if anchored else ne
            errs.append(e)
            nerrs.append(ne)
            rows.append({"held": "|".join(held), "arm": o.arm, "comparator": o.comparator, "interval": list(iv),
                         "anchored": anchored, "anchor": anchors.get(o.arm), "kappa": kappa,
                         "predicted": pred if anchored else null, "error": e, "null_predicted": null, "null_error": ne})
        folds.append({"held": "|".join(held), "n": len(arms), "kappa": kappa,
                      "mae": float(np.mean(errs)), "null_mae": float(np.mean(nerrs))})

    def summary(sel: list) -> dict:
        if not sel:
            return {"n": 0}
        mae, nmae = float(np.mean([r["error"] for r in sel])), float(np.mean([r["null_error"] for r in sel]))
        return {"n": len(sel), "mae": mae, "null_mae": nmae, "ratio": mae / nmae if nmae else float("nan")}

    lo, hi = E._paired_ci([f["mae"] - f["null_mae"] for f in folds])
    allr, anch = summary(rows), summary([r for r in rows if r["anchored"]])
    plac = summary([r for r in rows if r["comparator"] == "untreated"])
    green = bool(allr["ratio"] <= 0.8 and hi < 0 and plac["ratio"] <= 0.8)
    return {"fitted_kappa": fitted, "all": allr, "anchored_only": anch, "placebo_controlled": plac,
            "ci_fold_gap": [lo, hi], "green": green, "folds": folds, "rows": rows}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--anchors", type=Path, default=ANCHORS)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    anchors = load_anchors(args.anchors)
    model = Model()
    res = {"anchors_file": str(args.anchors), "anchors": anchors,
           "V3-P": run(model, anchors, fitted=False), "V3-K": run(model, anchors, fitted=True)}
    args.out.write_text(json.dumps(res, indent=2, default=float))
    fmt = lambda v: "inf" if not np.isfinite(v) else f"{v:.0f}"  # noqa: E731
    for name in ("V3-P", "V3-K"):
        r = res[name]
        print(f"\n{name}: S1 {r['all']['mae']:.1f}pp vs null {r['all']['null_mae']:.1f}pp "
              f"({100 * r['all']['ratio']:.0f} % of null; green needs 80 %), fold-gap CI "
              f"[{r['ci_fold_gap'][0]:+.1f}, {r['ci_fold_gap'][1]:+.1f}], placebo subset "
              f"{r['placebo_controlled']['mae']:.1f} vs {r['placebo_controlled']['null_mae']:.1f}, anchored arms only "
              f"(n={r['anchored_only']['n']}) {r['anchored_only'].get('mae', float('nan')):.1f} vs "
              f"{r['anchored_only'].get('null_mae', float('nan')):.1f}  ->  {'GREEN' if r['green'] else 'red'}")
        for row in r["rows"]:
            a = "  -  " if row["anchor"] is None else f"{row['anchor']:.2f}"
            print(f"    {row['arm']:<18} vs {row['comparator'][:13]:<13} truth ({fmt(row['interval'][0])},{fmt(row['interval'][1])})"
                  f"  anchor {a} k {row['kappa']:.2f}  pred {row['predicted']:6.1f}  err {row['error']:5.1f}  null {row['null_error']:5.1f}"
                  f"{'' if row['anchored'] else '  [unanchored: null]'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
