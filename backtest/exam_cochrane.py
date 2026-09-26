"""Exam C: the Cochrane 2024 network as a second exam of the cached curves.

Rules C1-C7 in docs/EXAM_COCHRANE_PREREG.md, committed before this module was
first run. This docstring only says where each rule lands.

    C1  COCHRANE_24M / NOT_IN_ARM_SET      which arms, from the review's SoF 2
    C2  interval() / anchor()              (RR - 1) x 100, CI as the interval
    C3  run_lomo()                         leave-one-dial-out, each model's dials
    C4  exam_v2's cached curves, FIT_GRID, no comparator adjustment
    C5  run_lomo() CS1, run_oracle() CS2, run_tau() CS3
    C6  read in main()'s printout order

Run:  PYTHONPATH=. python -m backtest.exam_cochrane
"""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import numpy as np

from backtest import exam_v2
from backtest.clinical import KNOWN_OUTCOMES
from backtest.exam_v2 import FIT_GRID, INF, _curve_value, _paired_ci, _tau, dist

RESULTS = Path(__file__).resolve().parent.parent / "results"
OUT = RESULTS / "exam_cochrane.json"

SOURCE = ("Gonzalez-Lorenzo et al. 2024, Cochrane Database Syst Rev 1:CD011381 "
          "(PMC10765473), Summary of findings 2, relapses at 24 months, network RR "
          "vs placebo with 95% CrI; open-access XML via Europe PMC, read 2026-09-26")

# C1. Arm -> (RR, lo, hi, GRADE certainty, direct evidence as the table states it).
# Names are this repo's arm labels; the review's row label is in the comment.
COCHRANE_24M: dict[str, tuple[float, float, float, str, str]] = {
    "alemtuzumab": (0.57, 0.47, 0.68, "high", "no direct evidence"),
    "cladribine": (0.53, 0.44, 0.64, "high", "1 RCT, 1326"),
    "dimethyl fumarate": (0.62, 0.55, 0.70, "moderate", "2 RCTs, 2307"),      # Dimethylfumarate
    "fingolimod": (0.54, 0.48, 0.60, "moderate", "2 RCTs, 2355"),
    "glatiramer acetate": (0.84, 0.76, 0.93, "moderate", "3 RCTs, 1014"),
    "IFN-beta-1b": (0.85, 0.76, 0.94, "low", "1 RCT, 372"),                   # Interferon beta-1b (Betaferon)
    "IFN-beta": (0.84, 0.78, 0.91, "moderate", "3 RCTs, 1629"),               # Interferon beta-1a (Avonex, Rebif)
    "natalizumab": (0.56, 0.48, 0.65, "high", "1 RCT, 942"),
    "ponesimod": (0.58, 0.48, 0.70, "moderate", "no direct evidence"),
    "teriflunomide": (0.82, 0.71, 0.94, "very low", "1 RCT, 1088"),
}

# Rows of the same table this arm set has no arm for. Listed, not scored (C1).
NOT_IN_ARM_SET: dict[str, tuple[float, float, float]] = {
    "azathioprine": (0.77, 0.51, 1.18),
    "immunoglobulins": (0.73, 0.59, 0.90),
    "interferon beta 1a-1b": (1.21, 0.66, 2.19),
    "laquinimod": (0.83, 0.76, 0.91),
    "mitoxantrone": (0.47, 0.27, 0.80),
}

MODELS: tuple[tuple[str, str], ...] = (
    ("velez", "s2"), ("minimal", "s2"), ("pernice", "s2"), ("pernice", "memread"),
)


def pct(rr: float) -> float:
    return (rr - 1.0) * 100.0


def interval(arm: str) -> tuple[float, float]:
    _, lo, hi, _, _ = COCHRANE_24M[arm]
    return (pct(lo), pct(hi))


def anchor(arm: str) -> float:
    return pct(COCHRANE_24M[arm][0])


def arms() -> list[str]:
    known = {o.arm for o in KNOWN_OUTCOMES}
    missing = set(COCHRANE_24M) - known
    if missing:
        raise RuntimeError(f"exam C names arms the arm set lacks: {sorted(missing)}")
    return list(COCHRANE_24M)


def groups() -> dict[tuple[str, ...], list[str]]:
    g: dict[tuple[str, ...], list[str]] = {}
    for a in arms():
        g.setdefault(exam_v2._pattern(a), []).append(a)
    return g


def _predict(arm: str, s: float, curve: dict) -> float:
    return _curve_value(exam_v2._pattern(arm), s, curve)


def _fit(training: list[str], curve: dict) -> float:
    best_s, best = 0.0, INF
    for s in FIT_GRID:
        sse = 0.0
        for a in training:
            p = _predict(a, s, curve)
            if np.isnan(p):
                sse = INF
                break
            sse += dist(p, interval(a)) ** 2
        if sse < best:
            best_s, best = s, sse
    return best_s


# CS1
def run_lomo(curve: dict) -> dict:
    g = groups()
    rows, folds = [], []
    for held, members in g.items():
        training = [a for p, ms in g.items() if p != held for a in ms]
        s = _fit(training, curve)
        # C3: a single fold has no training arms; the null is then undefined and
        # the exam has no resolving power for that model (C6), stated as nan.
        null = float(np.mean([anchor(a) for a in training])) if training else float("nan")
        errs, nerrs = [], []
        for a in members:
            pred = _predict(a, s, curve)
            iv = interval(a)
            e, ne = dist(pred, iv), dist(null, iv)
            errs.append(e)
            nerrs.append(ne)
            rows.append({"held": "|".join(held), "arm": a, "interval": list(iv),
                         "observed": anchor(a), "fitted_potency": s, "predicted": pred,
                         "error": e, "null_predicted": null, "null_error": ne})
        folds.append({"held": "|".join(held), "n": len(members), "fitted_potency": s,
                      "mae": float(np.mean(errs)), "null_mae": float(np.mean(nerrs))})
    lo, hi = _paired_ci([f["mae"] - f["null_mae"] for f in folds])
    return {"mae": float(np.mean([r["error"] for r in rows])),
            "null_mae": float(np.mean([r["null_error"] for r in rows])),
            "n_arms": len(rows), "n_folds": len(folds), "ci_fold_gap": [lo, hi],
            "folds": folds, "rows": rows}


# CS2
def run_oracle() -> dict:
    g = groups()
    grid = np.arange(-100.0, 100.0, 0.5)
    rows = []
    for pattern, members in g.items():
        ivs = [interval(a) for a in members]
        cost = [sum(dist(float(c), iv) for iv in ivs) for c in grid]
        centre = float(grid[int(np.argmin(cost))])
        others = [a for p, ms in g.items() if p != pattern for a in ms]
        null = float(np.mean([anchor(a) for a in others])) if others else float("nan")
        for a, iv in zip(members, ivs, strict=True):
            rows.append({"held": "|".join(pattern), "arm": a, "oracle": centre,
                         "error": dist(centre, iv), "null_predicted": null,
                         "null_error": dist(null, iv)})
    mae = float(np.mean([r["error"] for r in rows]))
    nmae = float(np.mean([r["null_error"] for r in rows]))
    return {"mae": mae, "null_mae": nmae, "headroom": nmae - mae,
            "relative": (nmae - mae) / nmae if nmae > 0 else float("nan"), "rows": rows}


# CS3
def run_tau(curve: dict) -> dict:
    s = _fit(arms(), curve)
    pred = {a: _predict(a, s, curve) for a in arms()}
    obs = {a: anchor(a) for a in arms()}
    pairs = [(pred[a], obs[a], pred[b], obs[b]) for a, b in combinations(arms(), 2)]
    return {"fitted_potency_all_arms": s, "tau": _tau(pairs), "n_pairs": len(pairs),
            "predicted": pred, "in_sample": True}


def main() -> int:
    print("exam C -- rules in docs/EXAM_COCHRANE_PREREG.md")
    print(f"  {len(arms())} arms scored; {len(NOT_IN_ARM_SET)} review rows without an arm: "
          f"{', '.join(NOT_IN_ARM_SET)}")
    out: dict = {"source": SOURCE, "arms": COCHRANE_24M, "not_in_arm_set": NOT_IN_ARM_SET,
                 "models": {}}
    for name, variant in MODELS:
        exam_v2.select_model(name, variant)
        curve = exam_v2.load_curve(verbose=False)
        cs1 = run_lomo(curve)
        cs2 = run_oracle()
        cs3 = run_tau(curve)
        key = name if variant == "s2" else f"{name}_{variant}"
        out["models"][key] = {"CS1": cs1, "CS2": cs2, "CS3": cs3}
        print(f"\n{key}: {cs1['n_folds']} dial folds over {cs1['n_arms']} arms")
        print(f"  CS1 interval-LOMO  MAE {cs1['mae']:.1f}pp vs null {cs1['null_mae']:.1f}pp  "
              f"[fold-gap 95% CI {cs1['ci_fold_gap'][0]:+.1f}, {cs1['ci_fold_gap'][1]:+.1f}]")
        for f in cs1["folds"]:
            print(f"      {f['held'] or '(none)':<40} n={f['n']}  s={f['fitted_potency']:.2f}  "
                  f"{f['mae']:6.1f} vs {f['null_mae']:5.1f}")
        print(f"  CS2 oracle         MAE {cs2['mae']:.1f}pp vs null {cs2['null_mae']:.1f}pp  "
              f"headroom {cs2['headroom']:.1f}pp ({cs2['relative']:.0%})")
        tau = cs3["tau"]
        tau_s = f"{tau:+.2f}" if np.isfinite(tau) else "undefined (all predictions tied)"
        print(f"  CS3 tau (in sample, s={cs3['fitted_potency_all_arms']:.2f})  "
              f"{tau_s} over {cs3['n_pairs']} pairs")
    OUT.write_text(json.dumps(out, indent=2, default=float))
    print(f"\nwritten {OUT.relative_to(RESULTS.parent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
