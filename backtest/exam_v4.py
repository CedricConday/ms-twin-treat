"""Exam v4: exam v3 with three registered changes (docs/EXAM_V4_PREREG.md).

    python -m backtest.exam_v4

C1  an active comparator's effect is its observed arm when that arm is outside the held
    fold, the model's otherwise;
C2  the interferons are unanchored;
C3  the S1P modulators carry one anchor per barrier-passage dial (effector memory T cells
    on pPass_BBB_teff, regulatory T cells on pPass_BBB_treg).
Nothing is fitted.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from backtest import exam_v2 as E
from backtest.ceiling_per_drug import _run
from backtest.exam_v3 import ANCHORS, CAP, Model, load_anchors
from bricks import qsp_pernice
from bricks.qsp_pernice import PerniceProfile
from bricks.sormani import predict_relapse_ratio

ROOT = Path(__file__).resolve().parent.parent
ANCHORS_V4 = ROOT / "data" / "pd_anchors_v4.tsv"
OUT = ROOT / "results" / "exam_v4.json"
INTERFERONS = ("IFN-beta", "IFN-beta-1b")
S1P = ("fingolimod", "ponesimod", "ozanimod")


def load_s1p(path: Path) -> dict[str, tuple[float, float]]:
    """Drug -> (effector anchor, regulatory anchor); only drugs with both sourced."""
    got: dict[str, dict[str, float]] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("status", "").strip() == "sourced" and row.get("anchor", "").strip():
                got.setdefault(row["drug"].strip(), {})[row["cell_type"].strip()] = min(CAP, max(0.0, float(row["anchor"])))
    return {d: (v["Teff"], v["Treg"]) for d, v in got.items() if d in S1P and {"Teff", "Treg"} <= set(v)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--anchors-v4", type=Path, default=ANCHORS_V4)
    ap.add_argument("--extra-anchors", type=Path, default=None,
                    help="further single-anchor drugs in pd_anchors.tsv format (exam v5)")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    anchors = {d: a for d, a in load_anchors(ANCHORS).items() if d not in INTERFERONS and d not in S1P}  # C2
    if args.extra_anchors is not None:
        anchors.update({d: a for d, a in load_anchors(args.extra_anchors).items() if d not in INTERFERONS and d not in S1P})
    s1p = load_s1p(args.anchors_v4)
    model = Model()
    untreated = _run(qsp_pernice.UNTREATED_PROFILE)

    def absolute(arm: str) -> float | None:
        """Model prediction of the arm against untreated at its anchor; None if unanchored."""
        if arm in s1p:  # C3
            ae, ar = s1p[arm]
            load = _run(PerniceProfile(label=arm, multipliers={"pPass_BBB_teff": 1.0 - ae, "pPass_BBB_treg": 1.0 - ar}))
            return None if load is None else predict_relapse_ratio(load / untreated).percent_change
        if arm in anchors:
            return model.absolute(arm, anchors[arm])
        return None

    outcomes = {o.arm: o for o in E.treated()}
    groups = E.all_groups()
    fold_of = {o.arm: held for held, os_ in groups.items() for o in os_}
    rows, folds = [], []
    for held, arms in groups.items():
        training = [o for p, os_ in groups.items() if p != held for o in os_]
        null = float(np.mean([E.anchor(o) for o in training]))
        errs, nerrs = [], []
        for o in arms:
            iv = E.interval(o)
            pred, how = absolute(o.arm), "model"
            if pred is not None and o.comparator != "untreated":
                c = outcomes.get(o.comparator)
                if c is not None and c.comparator == "untreated" and fold_of[c.arm] != held:  # C1
                    comp, how = E.anchor(c), "observed comparator"
                else:
                    comp, how = absolute(o.comparator), "modelled comparator"
                pred = None if comp is None else ((1.0 + pred / 100.0) / (1.0 + comp / 100.0) - 1.0) * 100.0
            anchored = bool(pred is not None and np.isfinite(pred))
            ne = E.dist(null, iv)
            e = E.dist(pred, iv) if anchored else ne
            errs.append(e)
            nerrs.append(ne)
            rows.append({"held": "|".join(held), "arm": o.arm, "comparator": o.comparator, "interval": list(iv),
                         "anchored": anchored, "anchor": list(s1p[o.arm]) if o.arm in s1p else anchors.get(o.arm),
                         "conversion": how if anchored else None, "predicted": pred if anchored else null,
                         "error": e, "null_predicted": null, "null_error": ne})
        folds.append({"held": "|".join(held), "n": len(arms), "mae": float(np.mean(errs)), "null_mae": float(np.mean(nerrs))})

    def summary(sel: list) -> dict:
        mae, nmae = float(np.mean([r["error"] for r in sel])), float(np.mean([r["null_error"] for r in sel]))
        return {"n": len(sel), "mae": mae, "null_mae": nmae, "ratio": mae / nmae if nmae else float("nan")}

    lo, hi = E._paired_ci([f["mae"] - f["null_mae"] for f in folds])
    allr, anch = summary(rows), summary([r for r in rows if r["anchored"]])
    plac = summary([r for r in rows if r["comparator"] == "untreated"])
    green = bool(allr["ratio"] <= 0.8 and hi < 0 and plac["ratio"] <= 0.8)
    res = {"anchors": anchors, "s1p_anchors": {k: list(v) for k, v in s1p.items()}, "all": allr, "anchored_only": anch,
           "placebo_controlled": plac, "ci_fold_gap": [lo, hi], "green": green, "folds": folds, "rows": rows}
    args.out.write_text(json.dumps(res, indent=2, default=float))
    fmt = lambda v: "inf" if not np.isfinite(v) else f"{v:.0f}"  # noqa: E731
    print(f"V4: S1 {allr['mae']:.1f}pp vs null {allr['null_mae']:.1f}pp ({100 * allr['ratio']:.0f} % of null; green needs 80 %), "
          f"fold-gap CI [{lo:+.1f}, {hi:+.1f}], placebo subset {plac['mae']:.1f} vs {plac['null_mae']:.1f} "
          f"({100 * plac['ratio']:.0f} %), anchored arms only (n={anch['n']}) {anch['mae']:.1f} vs {anch['null_mae']:.1f}"
          f"  ->  {'GREEN (candidate)' if green else 'red'}")
    for r in rows:
        a = "   -   " if r["anchor"] is None else ("/".join(f"{x:.2f}" for x in r["anchor"]) if isinstance(r["anchor"], list) else f"{r['anchor']:.2f}")
        print(f"    {r['arm']:<18} vs {r['comparator'][:13]:<13} truth ({fmt(r['interval'][0])},{fmt(r['interval'][1])})  anchor {a:<9}"
              f"  pred {r['predicted']:6.1f}  err {r['error']:5.1f}  null {r['null_error']:5.1f}  {r['conversion'] or '[unanchored: null]'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
