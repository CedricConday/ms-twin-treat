"""Leave-one-MECHANISM-out: can the model call a mechanism it has never seen?

`backtest/loo.py` holds out one ARM and fits on the rest. That grades
calibration: the held-out drug's class is still in the training set, so the fit
has seen something that works the same way. It is the right test for ranking
known drugs and the wrong one for a screen.

A screened candidate is usually a NEW mechanism. Nothing in training acts where
it acts. So the question that grades a screen is:

    hold out an entire mechanism GROUP, fit on the others, predict the group.

That is what this does, and it is deliberately harder than the LOO.

THE PIPELINE THIS RUNS ON
-------------------------
Not the ABM path the clinical gate uses. This runs end to end on the three
pieces built for screening, which is also a check that they compose:

    bricks/qsp_velez.py   the grounded MS QSP model; a drug is a set of
                          multipliers on its named rates
    bricks/profiles.py    which rates each arm moves, from pharmacology
    bricks/sormani.py     published trial-level map from a lesion rate ratio
                          to a relapse rate ratio (slope 0.52, R^2 0.71)

A simulated arm produces a damage trajectory. Damage relative to the untreated
arm on the SAME stochastic infection history is the lesion rate ratio. Sormani
turns that into a predicted relapse-rate change, which is what the trials
report. No step invents a scale.

WHAT IS FITTED
--------------
One global scalar, `s`, the potency every drug is assumed to have at its labelled
dose. Each arm's profile is rebuilt at that potency: a dial the pharmacology says
goes down becomes `1 - s`, a dial that goes up becomes `1 + s`, and dials the
drug does not touch stay at 1. Fitting one number across all mechanisms is the
honest analogue of the LOO's single class strength — and holding out a mechanism
asks whether that number transfers to a place it was never measured.

**Expect it not to.** A group held out is a dial with no training data on it at
all. If the headline MAE here is close to the null, that is not a bug and not a
calibration problem; it is the measurement this file exists to take.

THE ARM SET IS SMALL AND THE GROUPS ARE SMALLER
------------------------------------------------
Nine arms carry a quantified relapse number, and they fall into four mechanism
groups (see `mechanism_groups()`). One group is a single arm. With four groups
the fold count is four, so the headline is an average over four numbers and a
single arm can move it. Reported with the per-group detail, never alone.

Run:  PYTHONPATH=. python -m backtest.lomo
"""

from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from backtest.clinical import KNOWN_OUTCOMES
from bricks.profiles import PROFILES, touched_points
from bricks.qsp_velez import MechanismProfile, simulate
from bricks.sormani import predict_relapse_ratio

CACHE = Path(__file__).resolve().parent.parent / "results" / "mechanism_curve.json"

# Potency grid. 0 is untreated, and 1.0 is excluded: a dial multiplied by 0
# removes a rate entirely, which is outside anything the published model was
# characterised over.
POTENCY_GRID = [round(x, 2) for x in np.arange(0.0, 0.96, 0.05)]

# Cohort. Each seed is one stochastic infection history, shared between the
# treated and untreated arms so the ratio is a treatment effect and not a
# difference in luck.
#
# THE COHORT SIZE IS NOT A TASTE DECISION HERE, it is forced by the model's
# variance, and a first attempt at this file got it badly wrong. Measured over
# 40 five-year untreated runs: mean damage 22.7, SD 49.9, range 0.94 to 251 --
# a 266-fold spread, with the standard error at n=4 LARGER than the mean. The
# first version of this table used four seeds and was measuring nothing but its
# own noise. Only the alpha_R effect, which is about 100-fold, survived it.
#
# Two changes fix it:
#   * the statistic is the MEDIAN of per-seed ratios, not the mean. The damage
#     distribution is heavily right-tailed (at two years, median 1.43 against a
#     mean of 17.1), so the mean is dominated by a handful of runaway histories.
#   * 128 seeds. Bootstrapped, the median's coefficient of variation is ~16% at
#     n=48 and ~10% at n=128, which is the residual noise floor under every
#     number this file reports.
SEEDS = tuple(range(128))

# Two years, not the model's five-year FINAL TIME. The trials this is scored
# against report over 1-2 years, and the shorter horizon is both more comparable
# and far cheaper: runs are ~0.04s each here against ~0.43s at five years.
T_END = 730.0


def mechanism_groups() -> dict[tuple[str, ...], list[str]]:
    """Quantified arms, grouped by which intervention points they move.

    The group key IS the mechanism in this representation, which is why holding
    one out is holding out a mechanism rather than a label someone assigned.
    """
    groups: dict[tuple[str, ...], list[str]] = {}
    for outcome in KNOWN_OUTCOMES:
        if outcome.relapse_change_pct is None or outcome.arm == "untreated":
            continue
        key = touched_points(PROFILES[outcome.arm])
        groups.setdefault(key, []).append(outcome.arm)
    return groups


def _profile_at(arm: str, s: float) -> MechanismProfile:
    """Rebuild an arm's profile at potency `s`, keeping its sourced directions.

    Direction comes from bricks/profiles.py and is never changed here. Only the
    size moves, which is the thing being fitted.

    **This deliberately discards a FITTED magnitude where one exists.**
    Ocrelizumab carries ke=0.85 from Martinez-Pasamar 2013, and that value is
    NOT used here: it is replaced by the swept potency like every other arm.
    Keeping it would make its fold a different and easier test -- the held-out
    mechanism would arrive with a number somebody already measured for it, which
    is exactly the information LOMO is supposed to withhold. The fitted value is
    the right default everywhere else; inside this test it is contamination.
    """
    base = PROFILES[arm]
    points = {}
    for point in touched_points(base):
        points[point] = (1.0 - s) if getattr(base, point) < 1.0 else (1.0 + s)
    return MechanismProfile(label=f"{arm}@{s:g}", source="LOMO potency sweep", **points)


def _damage(profile: MechanismProfile, seed: int) -> float | None:
    """Total damage at the horizon, or None if the run left the model's regime."""
    traj = simulate(profile, t_end=T_END, seed=seed)
    if not traj["in_regime"]:
        return None
    return float(traj["total_damage"][-1])


def _median_ratio(changes: list[float]) -> float:
    """Median of PAIRED per-seed ratios. See the SEEDS comment for why median."""
    return float(np.median(changes))


def _cell(job: tuple[str, float, dict[int, float]]) -> tuple[str, float, list[float]]:
    """One (arm, potency) cell: every seed, returned as percent changes.

    Module-level and picklable so `build_table` can fan cells across processes.
    Each worker re-simulates the treated arm only; the untreated baselines are
    passed in, because they are shared by every cell and re-running them per
    cell would quadruple the work.
    """
    arm, s, untreated = job
    changes = []
    for seed, base in untreated.items():
        treated = _damage(_profile_at(arm, s), seed)
        if treated is None or treated <= 0.0:
            continue
        changes.append(predict_relapse_ratio(treated / base).percent_change)
    return arm, s, changes


def build_table(verbose: bool = True, workers: int | None = None) -> dict:
    """Tabulate predicted relapse change per (mechanism pattern, potency).

    Cached, because arms sharing a pattern share a curve: the five gamma_E arms
    are one column, not five. That is the same degeneracy `profiles.py` reports,
    showing up here as saved compute rather than as information.

    Fanned across processes. The grid is ~10k independent simulations and a
    serial build takes over ten minutes on four cores, which is long enough that
    people stop rebuilding it after changing a profile — and a stale response
    table is a silently wrong answer, not a slow one. Each cell is independent,
    so this is a pure fan-out with no shared state.
    """
    workers = workers or min(os.cpu_count() or 1, 8)

    untreated = {}
    for seed in SEEDS:
        d = _damage(PROFILES["untreated"], seed)
        if d is None or d <= 0.0:
            raise RuntimeError(
                f"untreated arm produced no damage at seed {seed}; the lesion ratio "
                "would be undefined. Widen the horizon or the cohort."
            )
        untreated[seed] = d

    patterns = sorted(mechanism_groups())
    groups = mechanism_groups()
    jobs = [(groups[pattern][0], s, untreated) for pattern in patterns for s in POTENCY_GRID]

    if verbose:
        print(f"  {len(jobs)} cells x {len(SEEDS)} seeds on {workers} workers ...")
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(_cell, jobs, chunksize=1))
    by_cell = {(arm, s): changes for arm, s, changes in results}

    table: dict[str, dict[str, dict]] = {}
    for pattern in patterns:
        arm = groups[pattern][0]   # any arm of the group; same dials
        col: dict[str, dict] = {}
        for s in POTENCY_GRID:
            changes = by_cell[(arm, s)]
            col[f"{s:g}"] = {
                # `mean` is the key name the rest of this module reads; the VALUE
                # is a median. Renaming it would be tidier, but the number that
                # matters is which statistic was computed, and that is stated
                # here and in the SEEDS comment.
                "mean": _median_ratio(changes) if changes else None,
                "sd": float(np.std(changes)) if len(changes) > 1 else None,
                "n_ok": len(changes),
                "n_out_of_regime": len(SEEDS) - len(changes),
            }
        table["|".join(pattern)] = col
        if verbose:
            print(f"  tabulated {'|'.join(pattern):<18} ({len(POTENCY_GRID)} potencies "
                  f"x {len(SEEDS)} seeds)")

    out = {
        "table": table,
        "potency_grid": POTENCY_GRID,
        "seeds": list(SEEDS),
        "t_end": T_END,
        "note": ("Predicted relapse-rate change per mechanism pattern and potency. "
                 "Damage from bricks/qsp_velez (Velez de Mendizabal 2011), lesion "
                 "ratio vs the untreated arm on the SAME seed, converted by the "
                 "Sormani & Bruzzi 2013 trial-level map. validated=False."),
    }
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(out, indent=2))
    return out


def load(rebuild: bool = False) -> dict:
    if rebuild or not CACHE.exists():
        print(f"building {CACHE.name} (a few minutes) ...")
        return build_table()
    return json.loads(CACHE.read_text())


def _predicted(arm: str, s: float, table: dict) -> float:
    """Interpolate this arm's predicted relapse change at potency `s`."""
    key = "|".join(touched_points(PROFILES[arm]))
    col = table["table"][key]
    grid = table["potency_grid"]
    xs = [g for g in grid if col[f"{g:g}"]["mean"] is not None]
    ys = [col[f"{g:g}"]["mean"] for g in xs]
    if not xs:
        return float("nan")
    return float(np.interp(s, xs, ys))


def _fit_potency(training_arms: list[str], known: dict[str, float], table: dict) -> float:
    """The single potency that best explains every training arm (least squares)."""
    best_s, best_sse = 0.0, float("inf")
    for s in np.arange(0.0, 0.951, 0.01):
        sse = 0.0
        ok = True
        for arm in training_arms:
            p = _predicted(arm, float(s), table)
            if np.isnan(p):
                ok = False
                break
            sse += (p - known[arm]) ** 2
        if ok and sse < best_sse:
            best_s, best_sse = float(s), sse
    return best_s


def run_lomo(table: dict | None = None) -> dict:
    table = table or load()
    groups = mechanism_groups()
    known = {o.arm: o.relapse_change_pct for o in KNOWN_OUTCOMES
             if o.relapse_change_pct is not None and o.arm != "untreated"}

    rows, folds = [], []
    for held_pattern, held_arms in groups.items():
        training_arms = [a for pat, arms in groups.items() if pat != held_pattern
                         for a in arms]
        s = _fit_potency(training_arms, known, table)
        null = float(np.mean([known[a] for a in training_arms]))

        fold_errors = []
        for arm in held_arms:
            pred = _predicted(arm, s, table)
            err = abs(pred - known[arm])
            rows.append({
                "held_mechanism": "|".join(held_pattern),
                "arm": arm,
                "known": known[arm],
                "fitted_potency": s,
                "predicted": pred,
                "error": err,
                "null_predicted": null,
                "null_error": abs(null - known[arm]),
            })
            fold_errors.append(err)
        folds.append({
            "mechanism": "|".join(held_pattern),
            "n_arms": len(held_arms),
            "mae": float(np.mean(fold_errors)),
            "null_mae": float(np.mean([abs(null - known[a]) for a in held_arms])),
        })

    return {
        "rows": rows,
        "folds": folds,
        "n_groups": len(groups),
        "mae": float(np.mean([r["error"] for r in rows])),
        "null_mae": float(np.mean([r["null_error"] for r in rows])),
    }


def main() -> int:
    r = run_lomo()
    print("LEAVE-ONE-MECHANISM-OUT — hold out a whole mechanism, predict it blind\n")
    print(f"{'held mechanism':<20} {'arm':<20} {'known':>8} {'potency':>8} "
          f"{'predicted':>10} {'error':>7} {'null err':>9}")
    print("-" * 88)
    for row in r["rows"]:
        print(f"{row['held_mechanism']:<20} {row['arm']:<20} {row['known']:>+7.1f}% "
              f"{row['fitted_potency']:>8.2f} {row['predicted']:>+9.1f}% "
              f"{row['error']:>6.1f} {row['null_error']:>8.1f}")
    print("-" * 88)
    print(f"{r['n_groups']} mechanism groups, {len(r['rows'])} quantified arms\n")
    for f in r["folds"]:
        print(f"  {f['mechanism']:<20} {f['n_arms']} arm(s)   "
              f"MAE {f['mae']:>6.1f}pp   null {f['null_mae']:>6.1f}pp")
    print(f"\nout-of-sample MAE: {r['mae']:.1f}pp   predict-the-mean null: {r['null_mae']:.1f}pp")

    beats = r["mae"] < r["null_mae"]
    print(f"\nThe model {'BEATS' if beats else 'DOES NOT BEAT'} the null on unseen mechanisms.")
    if not beats:
        print("  Expected, and this is the number the screen is gated on. A held-out")
        print("  mechanism is a dial with no training data on it: one global potency")
        print("  fitted elsewhere says nothing about a place it was never measured.")
        print("  Until this inverts, a generated candidate's score is not information.")
    print("\n  Four groups means four folds, and one group is a single arm — read the")
    print("  per-group rows, not the headline alone.")
    print(f"  Cohort is {len(SEEDS)} paired seeds at {T_END:.0f} days, scored on the MEDIAN")
    print("  per-seed ratio: untreated damage is heavily right-tailed (see SEEDS) and")
    print("  the bootstrapped noise floor on any number above is roughly 10%.")
    print("\n  Trial numbers are real and cited (docs/TRIAL_ANCHORS.md). Simulation")
    print("  numbers are proxies from a ported toy model, and the Sormani map is")
    print("  applied outside what it was fitted on. Nothing here is evidence about")
    print("  multiple sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
