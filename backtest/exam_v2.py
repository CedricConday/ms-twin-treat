"""Exam v2: every treated arm counts, scored as an interval.

Pre-registered in docs/EXAM_V2_PREREG.md, which fixes every rule below and was
committed before this module was first run. Read that file first; this
docstring only says where each rule lands in the code.

    R1  interval()        outcome -> [lo, hi]
    R2  anchor()          what the null averages per arm
    R3  all_groups()      folds over all 23 treated arms; missing response
                          columns built into results/exam_v2_curve.json
    R4  _predict()        comparator-adjusted (primary) or not (secondary)
    S1  run_interval_lomo()
    S2  run_oracle()
    S3  run_direction()
    S4  run_mri_rank()

Nothing here touches gate/criterion.py or backtest/lomo.py. Both exams run;
both are reported.

Run:  PYTHONPATH=. python -m backtest.exam_v2                 (Velez port, default)
      PYTHONPATH=. python -m backtest.exam_v2 --model minimal  (Jenner 2026 probe)

`--model` swaps the response-table builder and the arm-to-dial map
(docs/FOUR_DAY_PLAN.md, day 1, item 4) and nothing else: the rules R1-R4 and
the statistics S1-S4 are the registered ones and are shared. The minimal model
is deterministic, so its table is one run per cell, cached to
results/exam_v2_curve_minimal.json, and its result goes to
results/exam_v2_minimal.json. The Velez files are untouched by that switch.
"""

from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

import numpy as np

from backtest.clinical import KNOWN_OUTCOMES, MAG_TOLERANCE, NEUTRAL_BAND
from backtest.lomo import POTENCY_GRID, SEEDS, T_END, _cell, _damage, _median_ratio
from backtest.lomo import load as load_base
from backtest.potency import OBSERVED_LESION_RATIOS
from bricks import profiles_minimal, profiles_pernice, qsp_minimal, qsp_pernice
from bricks.profiles import PROFILES, touched_points
from bricks.sormani import predict_relapse_ratio

RESULTS = Path(__file__).resolve().parent.parent / "results"
CURVE = RESULTS / "exam_v2_curve.json"
OUT = RESULTS / "exam_v2.json"
CURVE_MINIMAL = RESULTS / "exam_v2_curve_minimal.json"
OUT_MINIMAL = RESULTS / "exam_v2_minimal.json"
CURVE_PERNICE = RESULTS / "exam_v2_curve_pernice.json"
OUT_PERNICE = RESULTS / "exam_v2_pernice.json"

# --------------------------------------------------------------------------- #
# The model under examination. "velez" is the transcription the exam was
# registered on; "minimal" is probe A; "pernice" is the port. Only the dial
# map and the table builder read this.
# --------------------------------------------------------------------------- #
MODELS = ("velez", "minimal", "pernice")
MODEL = {"name": "velez", "variant": "s2"}
OUT_BY_MODEL = {"velez": OUT, "minimal": OUT_MINIMAL, "pernice": OUT_PERNICE}


def select_model(name: str, variant: str = "s2") -> None:
    """`variant` is the Pernice port's arc-reading set (qsp_pernice.READING_VARIANTS);
    "s2" is the Figure S2 reproduction, the registered target."""
    if name not in MODELS:
        raise ValueError(f"unknown model {name!r}; one of {MODELS}")
    if variant not in qsp_pernice.READING_VARIANTS:
        raise ValueError(f"unknown variant {variant!r}")
    MODEL["name"] = name
    MODEL["variant"] = variant


def _pernice_paths() -> tuple[Path, Path]:
    v = MODEL["variant"]
    suffix = "" if v == "s2" else f"_{v}"
    return (RESULTS / f"exam_v2_curve_pernice{suffix}.json", RESULTS / f"exam_v2_pernice{suffix}.json")


def _pattern(arm: str) -> tuple[str, ...]:
    """The arm's dial pattern under the selected model. The fold key of S1/S2."""
    if MODEL["name"] == "minimal":
        return profiles_minimal.touched_points(profiles_minimal.PROFILES[arm])
    if MODEL["name"] == "pernice":
        return profiles_pernice.touched_points(profiles_pernice.PROFILES[arm])
    return touched_points(PROFILES[arm])

INF = float("inf")
FIT_GRID = [round(float(s), 2) for s in np.arange(0.0, 0.951, 0.01)]
BOOT = 10_000
RNG_SEED = 0


# --------------------------------------------------------------------------- #
# R1 / R2
# --------------------------------------------------------------------------- #
def interval(o) -> tuple[float, float]:
    if o.relapse_change_pct is not None:
        return (o.relapse_change_pct, o.relapse_change_pct)
    if o.direction == "neutral":
        return (-MAG_TOLERANCE, MAG_TOLERANCE)
    if o.direction == "harms":
        return (NEUTRAL_BAND, INF)
    if o.direction == "improves":
        return (-INF, -NEUTRAL_BAND)
    raise ValueError(o)


def anchor(o) -> float:
    if o.relapse_change_pct is not None:
        return o.relapse_change_pct
    return {"neutral": 0.0, "harms": NEUTRAL_BAND, "improves": -NEUTRAL_BAND}[o.direction]


def dist(p: float, iv: tuple[float, float]) -> float:
    lo, hi = iv
    if lo <= p <= hi:
        return 0.0
    return min(abs(p - lo), abs(p - hi))


def treated() -> list:
    return [o for o in KNOWN_OUTCOMES if o.arm != "untreated"]


def all_groups() -> dict[tuple[str, ...], list]:
    groups: dict[tuple[str, ...], list] = {}
    for o in treated():
        groups.setdefault(_pattern(o.arm), []).append(o)
    return groups


# --------------------------------------------------------------------------- #
# R3: the response table, extended to every pattern the arm set uses
# --------------------------------------------------------------------------- #
def _build_columns(patterns: list[tuple[str, ...]], groups: dict, verbose: bool = True) -> dict:
    """Same cells, same seeds, same statistic as backtest/lomo.build_table."""
    workers = min(os.cpu_count() or 1, 8)
    untreated = {}
    for seed in SEEDS:
        d = _damage(PROFILES["untreated"], seed)
        if d is None or d <= 0.0:
            raise RuntimeError(f"untreated arm produced no damage at seed {seed}")
        untreated[seed] = d
    jobs = [(groups[p][0].arm, s, untreated, None) for p in patterns for s in POTENCY_GRID]
    if verbose:
        print(f"  building {len(patterns)} pattern column(s): {len(jobs)} cells x "
              f"{len(SEEDS)} seeds on {workers} workers ...")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(_cell, jobs, chunksize=1))
    by_cell = {(arm, s): ch for arm, s, ch in results}
    cols = {}
    for p in patterns:
        arm = groups[p][0].arm
        col = {}
        for s in POTENCY_GRID:
            ch = by_cell[(arm, s)]
            col[f"{s:g}"] = {
                "mean": _median_ratio(ch) if ch else None,
                "sd": float(np.std(ch)) if len(ch) > 1 else None,
                "n_ok": len(ch),
                "n_out_of_regime": len(SEEDS) - len(ch),
            }
        cols["|".join(p)] = col
    return cols


def _build_columns_minimal(patterns: list[tuple[str, ...]], groups: dict,
                           verbose: bool = True) -> dict:
    """Probe A's table: deterministic model, one run per (pattern, potency) cell.

    Same cell semantics as `_build_columns`: the pattern's first arm is swept
    over POTENCY_GRID with `profiles_minimal.profile_at`, the lesion ratio goes
    through the unchanged Sormani map, and the column stores the percent
    change. `sd` is None and `n_ok` is 1 because there is no cohort: the model
    has no noise, so there is nothing to take a median over. The lesion ratio
    and the relapse proxy are kept beside the scored number for reading.
    """
    if verbose:
        print(f"  building {len(patterns)} pattern column(s) for the minimal model: "
              f"{len(patterns) * len(POTENCY_GRID)} deterministic runs ...")
    untreated = qsp_minimal.simulate(qsp_minimal.UNTREATED_PROFILE)
    cols = {}
    for p in patterns:
        arm = groups[p][0].arm
        col = {}
        for s in POTENCY_GRID:
            c = qsp_minimal.compare(profiles_minimal.profile_at(arm, s), untreated)
            ok = c["in_regime"] and np.isfinite(c["lesion_ratio"]) and c["lesion_ratio"] > 0
            col[f"{s:g}"] = {
                "mean": predict_relapse_ratio(c["lesion_ratio"]).percent_change if ok else None,
                "sd": None,
                "n_ok": 1 if ok else 0,
                "n_out_of_regime": 0 if ok else 1,
                "lesion_ratio": c["lesion_ratio"] if ok else None,
                "relapses_per_year": c["relapses_per_year"],
                "mean_M": c["mean_M"],
                "in_range": c["in_range"],
            }
        cols["|".join(p)] = col
    return cols


def load_curve_minimal(verbose: bool = True) -> dict:
    """Probe A's whole table lives in exam_v2_curve_minimal.json; no base file."""
    extra = json.loads(CURVE_MINIMAL.read_text()) if CURVE_MINIMAL.exists() else {"table": {}}
    table = dict(extra["table"])
    groups = all_groups()
    missing = [p for p in groups if "|".join(p) not in table]
    if missing:
        cols = _build_columns_minimal(missing, groups, verbose=verbose)
        extra["table"].update(cols)
        extra.update({
            "model": "minimal",
            "potency_grid": POTENCY_GRID,
            "seeds": [0],
            "schedule_months": {"dt": qsp_minimal.DT, "burn_in": qsp_minimal.BURN_IN_MONTHS,
                                "score": qsp_minimal.SCORE_MONTHS},
            "untreated_params": qsp_minimal.UNTREATED_PARAMS,
            "untreated_relapses_per_year": qsp_minimal.compare(qsp_minimal.UNTREATED_PROFILE)
            ["untreated_relapses_per_year"],
            "note": ("Jenner 2026 two-equation model (bricks/qsp_minimal.py), arms from "
                     "bricks/profiles_minimal.py, deterministic: one run per cell. The "
                     "scored number is the Sormani-mapped percent change of the "
                     "time-integrated inflammation ratio, treated over untreated."),
        })
        CURVE_MINIMAL.write_text(json.dumps(extra, indent=2))
        table.update(cols)
    return {"table": table, "potency_grid": POTENCY_GRID}


def _pernice_cell(job: tuple[str, float, float]) -> tuple[str, float, dict]:
    """One (arm, potency) cell of the port: a two-year deterministic run.

    Picklable so the table can fan across processes. The untreated lesion load
    is passed in. The scored number is the Sormani-mapped percent change of the
    lesion-load ratio (time-integrated ODC below Lmax, treated / untreated).
    """
    arm, s, untreated_load = job
    tr = qsp_pernice.simulate(profiles_pernice.profile_at(arm, s), sample_hours=24.0,
                              readings=qsp_pernice.READING_VARIANTS[MODEL["variant"]])
    ok = tr["in_regime"] and untreated_load > 0 and tr["lesion_load"] > 0
    ratio = tr["lesion_load"] / untreated_load if ok else None
    return arm, s, {
        "mean": predict_relapse_ratio(ratio).percent_change if ok else None,
        "sd": None,
        "n_ok": 1 if ok else 0,
        "n_out_of_regime": 0 if ok else 1,
        "lesion_ratio": ratio,
        "irreversibly_damaged": tr["irreversibly_damaged"],
    }


def _build_columns_pernice(patterns: list[tuple[str, ...]], groups: dict,
                           verbose: bool = True) -> dict:
    workers = min(os.cpu_count() or 1, 8)
    untreated = qsp_pernice.simulate(qsp_pernice.UNTREATED_PROFILE, sample_hours=24.0,
                                     readings=qsp_pernice.READING_VARIANTS[MODEL["variant"]])
    if not untreated["in_regime"] or untreated["lesion_load"] <= 0:
        raise RuntimeError("the untreated two-year run left the regime or did no damage")
    jobs = [(groups[p][0].arm, s, untreated["lesion_load"]) for p in patterns for s in POTENCY_GRID]
    if verbose:
        print(f"  building {len(patterns)} pattern column(s) for the Pernice port: "
              f"{len(jobs)} deterministic two-year runs on {workers} workers ...")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(_pernice_cell, jobs, chunksize=1))
    by_cell = {(arm, s): cell for arm, s, cell in results}
    cols = {}
    for p in patterns:
        arm = groups[p][0].arm
        cols["|".join(p)] = {f"{s:g}": by_cell[(arm, s)] for s in POTENCY_GRID}
    return cols, untreated


def load_curve_pernice(verbose: bool = True) -> dict:
    curve_path, _ = _pernice_paths()
    extra = json.loads(curve_path.read_text()) if curve_path.exists() else {"table": {}}
    table = dict(extra["table"])
    groups = all_groups()
    missing = [p for p in groups if "|".join(p) not in table]
    if missing:
        cols, untreated = _build_columns_pernice(missing, groups, verbose=verbose)
        extra["table"].update(cols)
        extra.update({
            "model": "pernice",
            "variant": MODEL["variant"],
            "potency_grid": POTENCY_GRID,
            "seeds": [0],
            "schedule": {"t_end_hours": qsp_pernice.TWO_YEAR_HOURS,
                         "injections_days": list(qsp_pernice.TWO_YEAR_INJECTIONS_DAYS)},
            "readings": qsp_pernice.READING_VARIANTS[MODEL["variant"]],
            "untreated": {"lesion_load": untreated["lesion_load"],
                          "irreversibly_damaged": untreated["irreversibly_damaged"]},
            "note": ("Pernice 2020 port (bricks/qsp_pernice.py), MS parameters, the paper's "
                     "two-year antigen schedule, arms from bricks/profiles_pernice.py, "
                     "deterministic. Scored number: Sormani-mapped percent change of the "
                     "lesion-load ratio (time-integrated ODC below Lmax, treated/untreated)."),
        })
        curve_path.write_text(json.dumps(extra, indent=2))
        table.update(cols)
    return {"table": table, "potency_grid": POTENCY_GRID}


def load_curve(verbose: bool = True) -> dict:
    """The transcription's own table plus the columns the widened exam needs.

    The base file is never written to. Extra columns live in exam_v2_curve.json
    and are rebuilt only if a pattern in the arm set has no column anywhere.
    Under `--model minimal` / `--model pernice` the whole table comes from that
    model's loader.
    """
    if MODEL["name"] == "minimal":
        return load_curve_minimal(verbose=verbose)
    if MODEL["name"] == "pernice":
        return load_curve_pernice(verbose=verbose)
    base = load_base()
    extra = json.loads(CURVE.read_text()) if CURVE.exists() else {"table": {}}
    table = dict(base["table"])
    table.update(extra["table"])
    groups = all_groups()
    missing = [p for p in groups if "|".join(p) not in table]
    if missing:
        cols = _build_columns(missing, groups, verbose=verbose)
        extra["table"].update(cols)
        extra.update({"potency_grid": POTENCY_GRID, "seeds": list(SEEDS), "t_end": T_END,
                      "note": ("Columns for dial patterns the quantified exam never "
                               "needed. Built by backtest/exam_v2 with the same cells, "
                               "seeds and median statistic as backtest/lomo.")})
        CURVE.write_text(json.dumps(extra, indent=2))
        table.update(cols)
    return {"table": table, "potency_grid": base["potency_grid"]}


def _curve_value(pattern: tuple[str, ...], s: float, curve: dict) -> float:
    col = curve["table"]["|".join(pattern)]
    xs = [g for g in curve["potency_grid"] if col[f"{g:g}"]["mean"] is not None]
    ys = [col[f"{g:g}"]["mean"] for g in xs]
    if not xs:
        return float("nan")
    return float(np.interp(s, xs, ys))


# --------------------------------------------------------------------------- #
# R4
# --------------------------------------------------------------------------- #
def _predict(o, s: float, curve: dict, adjust: bool) -> float:
    drug = _curve_value(_pattern(o.arm), s, curve)
    if not adjust or o.comparator == "untreated":
        return drug
    comp = _curve_value(_pattern(o.comparator), s, curve)
    return ((1.0 + drug / 100.0) / (1.0 + comp / 100.0) - 1.0) * 100.0


def _fit(training: list, curve: dict, adjust: bool) -> float:
    best_s, best = 0.0, INF
    for s in FIT_GRID:
        sse = 0.0
        for o in training:
            p = _predict(o, s, curve, adjust)
            if np.isnan(p):
                sse = INF
                break
            sse += dist(p, interval(o)) ** 2
        if sse < best:
            best_s, best = s, sse
    return best_s


def _paired_ci(diffs: list[float]) -> tuple[float, float]:
    if len(diffs) < 2:
        return (-INF, INF)
    rng = np.random.default_rng(RNG_SEED)
    d = np.asarray(diffs)
    means = [float(np.mean(rng.choice(d, size=len(d), replace=True))) for _ in range(BOOT)]
    return (float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5)))


# --------------------------------------------------------------------------- #
# S1
# --------------------------------------------------------------------------- #
def run_interval_lomo(curve: dict, adjust: bool = True) -> dict:
    groups = all_groups()
    rows, folds = [], []
    for held, arms in groups.items():
        training = [o for p, os_ in groups.items() if p != held for o in os_]
        s = _fit(training, curve, adjust)
        null = float(np.mean([anchor(o) for o in training]))
        errs, nerrs = [], []
        for o in arms:
            pred = _predict(o, s, curve, adjust)
            iv = interval(o)
            e, ne = dist(pred, iv), dist(null, iv)
            errs.append(e)
            nerrs.append(ne)
            rows.append({"held": "|".join(held), "arm": o.arm, "interval": list(iv),
                         "fitted_potency": s, "predicted": pred, "error": e,
                         "null_predicted": null, "null_error": ne})
        folds.append({"held": "|".join(held), "n": len(arms), "fitted_potency": s,
                      "mae": float(np.mean(errs)), "null_mae": float(np.mean(nerrs))})
    mae = float(np.mean([r["error"] for r in rows]))
    nmae = float(np.mean([r["null_error"] for r in rows]))
    lo, hi = _paired_ci([f["mae"] - f["null_mae"] for f in folds])
    placebo = [r for r in rows if KNOWN_OUTCOMES_BY[r["arm"]].comparator == "untreated"]
    return {"adjusted": adjust, "mae": mae, "null_mae": nmae, "n_arms": len(rows),
            "n_folds": len(folds), "ci_fold_gap": [lo, hi],
            "placebo_mae": float(np.mean([r["error"] for r in placebo])),
            "placebo_null_mae": float(np.mean([r["null_error"] for r in placebo])),
            "folds": folds, "rows": rows}


# --------------------------------------------------------------------------- #
# S2
# --------------------------------------------------------------------------- #
def run_oracle() -> dict:
    groups = all_groups()
    grid = np.arange(-100.0, 200.0, 0.5)
    rows = []
    for pattern, arms in groups.items():
        ivs = [interval(o) for o in arms]
        cost = [sum(dist(float(c), iv) for iv in ivs) for c in grid]
        centre = float(grid[int(np.argmin(cost))])
        others = [o for p, os_ in groups.items() if p != pattern for o in os_]
        null = float(np.mean([anchor(o) for o in others]))
        for o, iv in zip(arms, ivs, strict=True):
            rows.append({"held": "|".join(pattern), "arm": o.arm, "oracle": centre,
                         "error": dist(centre, iv), "null_predicted": null,
                         "null_error": dist(null, iv)})
    mae = float(np.mean([r["error"] for r in rows]))
    nmae = float(np.mean([r["null_error"] for r in rows]))
    return {"mae": mae, "null_mae": nmae, "headroom": nmae - mae,
            "relative": (nmae - mae) / nmae if nmae > 0 else float("nan"), "rows": rows}


# --------------------------------------------------------------------------- #
# S3
# --------------------------------------------------------------------------- #
def _cls(p: float) -> str:
    if p < -NEUTRAL_BAND:
        return "improves"
    if p > NEUTRAL_BAND:
        return "harms"
    return "neutral"


def _balanced_accuracy(true: list[str], pred: list[str]) -> float:
    recalls = []
    for c in ("improves", "neutral", "harms"):
        idx = [i for i, t in enumerate(true) if t == c]
        if idx:
            recalls.append(np.mean([pred[i] == c for i in idx]))
    return float(np.mean(recalls))


def run_direction(lomo: dict) -> dict:
    true = [KNOWN_OUTCOMES_BY[r["arm"]].direction for r in lomo["rows"]]
    pred = [_cls(r["predicted"]) for r in lomo["rows"]]
    acc = float(np.mean([t == p for t, p in zip(true, pred, strict=True)]))
    bacc = _balanced_accuracy(true, pred)
    majority = max(set(true), key=true.count)
    maj_acc = float(np.mean([t == majority for t in true]))
    maj_bacc = _balanced_accuracy(true, [majority] * len(true))
    rng = np.random.default_rng(RNG_SEED)
    t = np.asarray(true)
    perm = [_balanced_accuracy(list(rng.permutation(t)), pred) for _ in range(BOOT)]
    p = float(np.mean([b >= bacc for b in perm]))
    return {"accuracy": acc, "balanced_accuracy": bacc, "majority_class": majority,
            "majority_accuracy": maj_acc, "majority_balanced_accuracy": maj_bacc,
            "permutation_p": p,
            "confusion": {f"{t}->{pr}": sum(1 for a, b in zip(true, pred, strict=True) if a == t and b == pr)
                          for t in ("improves", "neutral", "harms")
                          for pr in ("improves", "neutral", "harms")}}


# --------------------------------------------------------------------------- #
# S4
# --------------------------------------------------------------------------- #
def _metric(source: str) -> str:
    return source[source.rfind("[") + 1:source.rfind("]")] if "[" in source else "?"


def _tau(pairs: list[tuple[float, float, float, float]]) -> float:
    c = d = 0
    for a1, b1, a2, b2 in pairs:
        s = np.sign((a1 - a2) * (b1 - b2))
        if s > 0:
            c += 1
        elif s < 0:
            d += 1
    return (c - d) / (c + d) if (c + d) else float("nan")


def run_mri_rank() -> dict:
    arms = {o.arm: o for o in treated()
            if o.relapse_change_pct is not None and o.arm in OBSERVED_LESION_RATIOS}
    groups: dict[tuple[str, ...], list[str]] = {}
    for a in arms:
        groups.setdefault(_pattern(a), []).append(a)

    def pairs(same_metric: bool):
        out = []
        for members in groups.values():
            for a, b in combinations(members, 2):
                ma, mb = (_metric(OBSERVED_LESION_RATIOS[x][1]) for x in (a, b))
                if same_metric and ma != mb:
                    continue
                out.append((a, b))
        return out

    def score(pairlist, ratio_of):
        return _tau([(np.log(OBSERVED_LESION_RATIOS[a][0]), ratio_of[a],
                      np.log(OBSERVED_LESION_RATIOS[b][0]), ratio_of[b])
                     for a, b in pairlist])

    result = {}
    for label, same in (("all_pairs", False), ("same_metric_pairs", True)):
        pl = pairs(same)
        observed = {a: arms[a].relapse_change_pct for a in arms}
        tau = score(pl, observed)
        rng = np.random.default_rng(RNG_SEED)
        perms = []
        for _ in range(BOOT):
            shuffled = {}
            for members in groups.values():
                vals = [observed[m] for m in members]
                for m, v in zip(members, rng.permutation(vals), strict=True):
                    shuffled[m] = float(v)
            perms.append(score(pl, shuffled))
        perms = [x for x in perms if not np.isnan(x)]
        p = float(np.mean([x >= tau for x in perms])) if perms and not np.isnan(tau) else float("nan")
        result[label] = {"n_pairs": len(pl), "tau": tau, "permutation_p": p,
                         "pairs": [f"{a} vs {b}" for a, b in pl]}
    return result


KNOWN_OUTCOMES_BY = {o.arm: o for o in KNOWN_OUTCOMES}


# --------------------------------------------------------------------------- #
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", choices=MODELS, default="velez",
                    help="which model's dial map and response table to score (default: velez)")
    ap.add_argument("--variant", choices=tuple(qsp_pernice.READING_VARIANTS), default="s2",
                    help="Pernice port only: which arc-reading set (default: s2, the Figure S2 one)")
    args = ap.parse_args(argv)
    select_model(args.model, args.variant)
    out = _pernice_paths()[1] if args.model == "pernice" else OUT_BY_MODEL[args.model]
    tag = f" (variant {args.variant})" if args.model == "pernice" else ""
    print(f"exam v2 -- rules in docs/EXAM_V2_PREREG.md -- model: {args.model}{tag}")
    curve = load_curve()
    s1 = run_interval_lomo(curve, adjust=True)
    s1u = run_interval_lomo(curve, adjust=False)
    s2 = run_oracle()
    s3 = run_direction(s1)
    s4 = run_mri_rank()

    print(f"\nS1 interval-LOMO, comparator-adjusted (primary): "
          f"MAE {s1['mae']:.1f}pp vs null {s1['null_mae']:.1f}pp  "
          f"[fold-gap 95% CI {s1['ci_fold_gap'][0]:+.1f}, {s1['ci_fold_gap'][1]:+.1f}]  "
          f"placebo-only {s1['placebo_mae']:.1f} vs {s1['placebo_null_mae']:.1f}")
    for f in s1["folds"]:
        print(f"    {f['held']:<18} n={f['n']}  s={f['fitted_potency']:.2f}  "
              f"{f['mae']:6.1f} vs {f['null_mae']:5.1f}")
    print(f"S1 unadjusted (lomo.py convention): MAE {s1u['mae']:.1f}pp vs null {s1u['null_mae']:.1f}pp")
    print(f"\nS2 oracle (perfect dial-level model): MAE {s2['mae']:.1f}pp vs null "
          f"{s2['null_mae']:.1f}pp  headroom {s2['headroom']:.1f}pp ({s2['relative']:.0%})")
    print(f"\nS3 direction: accuracy {s3['accuracy']:.0%}, balanced {s3['balanced_accuracy']:.2f}; "
          f"majority null {s3['majority_accuracy']:.0%} / {s3['majority_balanced_accuracy']:.2f}; "
          f"permutation p = {s3['permutation_p']:.3f}")
    for k, v in s3["confusion"].items():
        if v:
            print(f"    {k:<20} {v}")
    print("\nS4 MRI channel, within-dial rank:")
    for k, v in s4.items():
        print(f"    {k:<18} n_pairs={v['n_pairs']:2d}  tau={v['tau']:+.2f}  p={v['permutation_p']:.3f}")

    out.write_text(json.dumps({"model": args.model, "variant": args.variant, "S1": s1,
                               "S1_unadjusted": s1u, "S2": s2, "S3": s3, "S4": s4},
                              indent=2, default=float))
    print(f"\nwritten {out.relative_to(RESULTS.parent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
