"""Exam H, the confirmation exam (docs/EXAM_H_PREREG.md).

    python -m backtest.exam_h predict     # writes results/exam_h_predictions.json; reads no outcome
    python -m backtest.exam_h score       # reads data/confirmation_outcomes.tsv and the predictions

The pipeline is exam v5's, frozen: the anchor is the magnitude, nothing is fitted, an
active comparator's effect is its observed arm among the 23 exam arms, interferons are
unanchored, S1P modulators carry two anchors. `predict` must be run and committed before
the outcome file is opened.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

from backtest import exam_v2 as E
from backtest.ceiling_per_drug import _run
from backtest.exam_v3 import ANCHORS, CAP, Model, load_anchors
from backtest.exam_v4 import ANCHORS_V4, INTERFERONS
from bricks import qsp_pernice
from bricks.qsp_pernice import PerniceProfile
from bricks.sormani import predict_relapse_ratio

ROOT = Path(__file__).resolve().parent.parent
ANCHORS_H = ROOT / "data" / "pd_anchors_h.tsv"
OUTCOMES = ROOT / "data" / "confirmation_outcomes.tsv"
PRED = ROOT / "results" / "exam_h_predictions.json"
OUT = ROOT / "results" / "exam_h.json"

# (stratum, drug, trial, comparator, profile the drug is run on). The profile is the exam
# arm whose dials the registration assigns: "s1p" = both passage dials with two anchors.
ARMS = [
    ("A", "siponimod", "BOLD", "untreated", "s1p"),
    ("A", "mitoxantrone", "Millefiorini 1997", "untreated", "cladribine"),
    ("A", "evobrutinib", "evolutionRMS 1", "teriflunomide", "abatacept"),
    ("A", "evobrutinib", "evolutionRMS 2", "teriflunomide", "abatacept"),
    ("A", "tolebrutinib", "GEMINI 1", "teriflunomide", "abatacept"),
    ("A", "tolebrutinib", "GEMINI 2", "teriflunomide", "abatacept"),
    ("A", "fenebrutinib", "FENhance 1", "teriflunomide", "abatacept"),
    ("A", "fenebrutinib", "FENhance 2", "teriflunomide", "abatacept"),
    ("A", "peginterferon beta-1a", "ADVANCE", "untreated", None),
    ("B", "fingolimod", "FREEDOMS II", "untreated", "s1p"),
    ("B", "fingolimod", "TRANSFORMS", "IFN-beta", "s1p"),
    ("B", "teriflunomide", "TOWER", "untreated", "teriflunomide"),
    ("B", "dimethyl fumarate", "CONFIRM", "untreated", "dimethyl fumarate"),
    ("B", "ocrelizumab", "OPERA II", "IFN-beta", "ocrelizumab"),
    ("B", "ofatumumab", "ASCLEPIOS II", "teriflunomide", "ofatumumab"),
    ("B", "ublituximab", "ULTIMATE II", "teriflunomide", "ublituximab"),
    ("B", "alemtuzumab", "CARE-MS II", "IFN-beta", "alemtuzumab"),
    ("B", "ozanimod", "SUNBEAM", "IFN-beta", "s1p"),
    ("B", "daclizumab", "SELECT", "untreated", None),
    ("B", "rituximab", "HERMES", "untreated", "rituximab"),
]
# Firategrast is dropped under the registration's own rule: its source reports a relapse
# rate per treatment period, not an annualised rate (decided 2026-11-23 before the outcome
# file was opened, from the collector's description of the source alone).
DROPPED = {"firategrast": "relapse rate per treatment period, not annualised"}


def _two_anchor(path: Path, key: str) -> dict[str, tuple[float, float]]:
    """Drug -> (effector, regulatory) from a file whose `key` column names the cell type."""
    got: dict[str, dict[str, float]] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cell = row.get(key, "").strip()[:4]
            if row.get("status", "").strip() == "sourced" and row.get("anchor", "").strip() and cell in ("Teff", "Treg"):
                got.setdefault(row["drug"].strip(), {})[cell] = min(CAP, max(0.0, float(row["anchor"])))
    return {d: (v["Teff"], v["Treg"]) for d, v in got.items() if {"Teff", "Treg"} <= set(v)}


def predict() -> int:
    single = {d: a for d, a in load_anchors(ANCHORS).items() if d not in INTERFERONS}
    with open(ANCHORS_H, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("status", "").strip() == "sourced" and row.get("anchor", "").strip()
                    and not row.get("quantity", "").startswith(("Teff", "Treg"))):
                single[row["drug"].strip()] = min(CAP, max(0.0, float(row["anchor"])))
    double = {**_two_anchor(ANCHORS_V4, "cell_type"), **_two_anchor(ANCHORS_H, "quantity")}
    model = Model()
    untreated = _run(qsp_pernice.UNTREATED_PROFILE)
    observed = {o.arm: E.anchor(o) for o in E.treated() if o.comparator == "untreated"}
    null_abs = float(np.mean([E.anchor(o) for o in E.treated()]))

    def convert(vs_untreated: float, comparator: str) -> float:
        if comparator == "untreated":
            return vs_untreated
        return ((1.0 + vs_untreated / 100.0) / (1.0 + observed[comparator] / 100.0) - 1.0) * 100.0

    rows = []
    for stratum, drug, trial, comparator, profile in ARMS:
        absolute, anchor = None, None
        if profile == "s1p" and drug in double:
            anchor = list(double[drug])
            load = _run(PerniceProfile(label=drug, multipliers={"pPass_BBB_teff": 1.0 - anchor[0], "pPass_BBB_treg": 1.0 - anchor[1]}))
            absolute = None if load is None else predict_relapse_ratio(load / untreated).percent_change
        elif profile not in (None, "s1p") and drug in single:
            anchor = single[drug]
            absolute = model.absolute(profile, anchor)
        anchored = bool(absolute is not None and np.isfinite(absolute))
        null = convert(null_abs, comparator)
        rows.append({"stratum": stratum, "drug": drug, "trial": trial, "comparator": comparator, "profile": profile,
                     "anchor": anchor, "anchored": anchored,
                     "predicted": convert(absolute, comparator) if anchored else null, "null_predicted": null})
    PRED.write_text(json.dumps({"pipeline": "exam v5 at 85a32fc", "null_vs_untreated": null_abs,
                                "observed_comparators": observed, "dropped": DROPPED, "rows": rows}, indent=2, default=float))
    for r in rows:
        a = "-" if r["anchor"] is None else ("/".join(f"{x:.2f}" for x in r["anchor"]) if isinstance(r["anchor"], list) else f"{r['anchor']:.2f}")
        print(f"  {r['stratum']} {r['drug']:<22} {r['trial']:<18} vs {r['comparator']:<13} anchor {a:<9} "
              f"predicted {r['predicted']:6.1f}  null {r['null_predicted']:6.1f}{'' if r['anchored'] else '  [unanchored: null]'}")
    print(f"wrote {PRED}")
    return 0


def score() -> int:
    pred = json.loads(PRED.read_text())
    truth = {}
    with open(OUTCOMES, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("status", "").strip() == "sourced" and row.get("relapse_change_pct", "").strip():
                truth[(row["drug"].strip().lower(), row["trial"].strip().lower())] = float(row["relapse_change_pct"])

    def find(drug: str, trial: str) -> float | None:
        for (d, t), v in truth.items():
            if drug.lower().split()[0] in d and trial.lower().replace(" ", "") in t.replace(" ", ""):
                return v
        return None

    rows, missing = [], []
    for r in pred["rows"]:
        t = find(r["drug"], r["trial"])
        if t is None:
            missing.append(f"{r['drug']} {r['trial']}")
            continue
        rows.append({**r, "truth": t, "error": abs(r["predicted"] - t), "null_error": abs(r["null_predicted"] - t)})

    def summary(sel: list) -> dict:
        if not sel:
            return {"n": 0}
        mae, nmae = float(np.mean([x["error"] for x in sel])), float(np.mean([x["null_error"] for x in sel]))
        rng = np.random.default_rng(E.RNG_SEED)
        d = np.array([x["error"] - x["null_error"] for x in sel])
        boots = [float(np.mean(rng.choice(d, size=len(d), replace=True))) for _ in range(E.BOOT)]
        return {"n": len(sel), "mae": mae, "null_mae": nmae, "ratio": mae / nmae,
                "ci_gap": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]}

    res = {"all": summary(rows), "stratum_A": summary([x for x in rows if x["stratum"] == "A"]),
           "stratum_B": summary([x for x in rows if x["stratum"] == "B"]),
           "anchored_only": summary([x for x in rows if x["anchored"]]), "missing": missing, "rows": rows}
    OUT.write_text(json.dumps(res, indent=2, default=float))
    for name in ("all", "stratum_A", "stratum_B", "anchored_only"):
        s = res[name]
        if s["n"]:
            print(f"{name:<14} n={s['n']:<3} S1 {s['mae']:5.1f}pp vs null {s['null_mae']:5.1f}pp ({100 * s['ratio']:.0f} % of null), "
                  f"gap CI [{s['ci_gap'][0]:+.1f}, {s['ci_gap'][1]:+.1f}]")
    for x in rows:
        print(f"  {x['stratum']} {x['drug']:<22} {x['trial']:<18} vs {x['comparator']:<13} truth {x['truth']:6.1f}  pred {x['predicted']:6.1f}"
              f"  err {x['error']:5.1f}  null {x['null_error']:5.1f}{'' if x['anchored'] else '  [unanchored]'}")
    if missing:
        print("no outcome matched:", "; ".join(missing))
    return 0


if __name__ == "__main__":
    raise SystemExit({"predict": predict, "score": score}[sys.argv[1]]())
