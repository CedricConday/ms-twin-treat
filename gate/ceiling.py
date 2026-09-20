"""The oracle ceiling: the best this model form could EVER do on these arms.

`backtest/lomo.py` fits ONE potency across all arms and asks whether it
transfers to a mechanism it never saw. It does not (45.9pp vs a 12.3pp null).
The obvious next move is to get per-drug potency from independent data, which is
what `backtest/potency.py` does through the MRI channel and what BUILD_PLAN
blockers (4) and (6) are about.

This module asks whether that move can possibly work, before anyone spends
another week on it:

    give every arm its OWN potency, fitted IN SAMPLE to its own trial number --
    deliberately cheating, the best a perfect potency oracle could ever supply --
    and measure the error that remains.

That number is a CEILING. No independent potency source, however good, can beat
a potency fitted directly to the answer. So:

    oracle beats the null  ->  calibration is the blocker; better potency data
                               is worth pursuing, and the gap between the oracle
                               and the MRI-fitted potencies is the size of the
                               prize
    oracle loses anyway    ->  the model FORM is the blocker; no amount of
                               potency data lifts the gate, and blockers (4) and
                               (6) are closed routes rather than deferred ones

REACHABLE AND UNREACHABLE ARMS ARE REPORTED SEPARATELY, AND THAT IS NOT OPTIONAL
---------------------------------------------------------------------------------
For five of the twelve quantified arms the model cannot produce a benefit at ANY
potency -- `gamma_E` and a lowered `alpha_R` both RAISE damage here, so the whole
depleting/sequestering class is structurally unable to come out beneficial
(`bricks/qsp_velez.py`, and `backtest/potency.py` reports them OUT OF RANGE).
On those arms the oracle's "best fit" is a boundary value, not a fit, and its
error is a sign error.

Averaging both kinds together produces a single number that reads like a
calibration result while actually being dominated by direction failures. So the
headline here is the REACHABLE subset, the unreachable arms are reported beside
it, and the combined figure is shown last and labelled as what it is.

THE CEILING ALONE IS NEARLY TAUTOLOGICAL, AND THE SECOND HALF IS THE REAL TEST
-------------------------------------------------------------------------------
One free parameter fitted to one target reproduces that target. So "the oracle
beats the null on reachable arms" is close to arithmetic, and read alone it
would be a fine way to talk yourself into a calibration project. What it does
establish is narrower and still worth having: on those arms the model's response
curve PASSES THROUGH the trial's number, so nothing about the curve's shape
forbids a correct answer there. On the unreachable arms it does not pass through
at any potency, and that is a fact about the model, not the fit.

The question that actually decides blockers (4) and (6) is RECOVERABILITY:

    does an INDEPENDENT potency source land near the oracle's potency?

`backtest/potency.py` supplies one -- each drug's magnitude fitted to its trial's
MRI lesion ratio, never to its relapse number. Scoring the arms at those
MRI-fitted potencies is a real out-of-fit prediction, and `recoverability()`
below reports it beside the oracle. If the MRI potencies reproduce the relapse
numbers, the channel works and the gate is a calibration problem. If they are
systematically off, the independent source cannot deliver what the model needs,
whatever its own quality.

COST
----
None to speak of. The response tables are already cached
(`results/mechanism_curve*.json`, ~8 minutes each to build, built by the
mechanism-curve runs), and an in-sample per-arm oracle is an argmin over a
20-point grid with interpolation. Nothing is simulated here.

Run:  PYTHONPATH=. python3 -m gate.ceiling
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from backtest.clinical import KNOWN_OUTCOMES
from backtest.lomo import POTENCY_GRID, _predicted, load

# A fit that lands on the first or last grid point is not a fit: the optimum is
# outside the model's reachable range and the grid edge is merely the closest it
# can get. `backtest/potency.py` calls the same condition OUT OF RANGE.
GRID_LO, GRID_HI = min(POTENCY_GRID), max(POTENCY_GRID)
EDGE_TOL = 1e-6


@dataclass(frozen=True)
class OracleFit:
    """One arm's best possible outcome under a perfect potency oracle."""

    arm: str
    known: float
    potency: float
    predicted: float
    error: float
    null_error: float
    reachable: bool

    def line(self) -> str:
        flag = "" if self.reachable else "   AT GRID EDGE — unreachable, not fitted"
        return (f"{self.arm:<22} {self.known:>+7.1f}% {self.potency:>8.2f} "
                f"{self.predicted:>+10.1f}% {self.error:>7.1f} {self.null_error:>9.1f}{flag}")


def _oracle_potency(arm: str, target: float, table: dict) -> tuple[float, float]:
    """The in-sample best potency for one arm, and what it predicts there.

    Searched on a finer grid than POTENCY_GRID because interpolation between
    tabulated points is free -- a coarse argmin would understate the oracle and
    so overstate the case against calibration.
    """
    best_s, best_err, best_pred = GRID_LO, float("inf"), float("nan")
    for s in np.arange(GRID_LO, GRID_HI + 1e-9, 0.005):
        pred = _predicted(arm, float(s), table)
        if np.isnan(pred):
            continue
        err = abs(pred - target)
        if err < best_err:
            best_s, best_err, best_pred = float(s), err, pred
    return best_s, best_pred


def run_ceiling(table: dict | None = None) -> dict:
    table = table or load()
    # The same arm set backtest/lomo.py scores: quantified, and not the control.
    known = {o.arm: o.relapse_change_pct for o in KNOWN_OUTCOMES
             if o.relapse_change_pct is not None and o.arm != "untreated"}

    # The null every scorer in this repo is measured against: answer with the
    # mean of the OTHER arms' outcomes. Computed leave-one-out so the oracle and
    # the null are asked the same question about the same arm.
    fits = []
    for arm, target in known.items():
        s, pred = _oracle_potency(arm, target, table)
        others = [v for a, v in known.items() if a != arm]
        null_pred = float(np.mean(others))
        reachable = (GRID_LO + EDGE_TOL) < s < (GRID_HI - EDGE_TOL)
        fits.append(OracleFit(
            arm=arm, known=target, potency=s, predicted=pred,
            error=abs(pred - target), null_error=abs(null_pred - target),
            reachable=reachable))

    reach = [f for f in fits if f.reachable]
    unreach = [f for f in fits if not f.reachable]

    def _mae(rows, attr):
        return float(np.mean([getattr(r, attr) for r in rows])) if rows else float("nan")

    return {
        "fits": fits,
        "n_reachable": len(reach),
        "n_unreachable": len(unreach),
        "mae_reachable": _mae(reach, "error"),
        "null_mae_reachable": _mae(reach, "null_error"),
        "mae_unreachable": _mae(unreach, "error"),
        "null_mae_unreachable": _mae(unreach, "null_error"),
        "mae_all": _mae(fits, "error"),
        "null_mae_all": _mae(fits, "null_error"),
    }


def recoverability(table: dict | None = None) -> dict:
    """Do the INDEPENDENT (MRI-fitted) potencies land where the oracle needs them?

    The oracle says what potency reproduces each arm's relapse number.
    `backtest/potency.py` says what potency reproduces each arm's MRI lesion
    ratio, having never seen the relapse number. Agreement between the two is
    the whole case for the potency channel; disagreement is its size.

    Scored only on arms where BOTH are in range -- an arm the model cannot reach
    has no oracle potency to compare against, and an arm the MRI fitter called
    OUT OF RANGE has no independent estimate.
    """
    from backtest.potency import fit_all

    table = table or load()
    ceil = run_ceiling(table)
    oracle = {f.arm: f for f in ceil["fits"]}
    fitted = fit_all(table)

    rows = []
    for arm, fp in fitted.items():
        o = oracle.get(arm)
        if o is None or not o.reachable or not fp.in_range:
            continue
        pred = _predicted(arm, fp.potency, table)
        rows.append({
            "arm": arm,
            "known": o.known,
            "oracle_potency": o.potency,
            "mri_potency": fp.potency,
            "ratio": (fp.potency / o.potency) if o.potency else float("nan"),
            "predicted_at_mri": pred,
            "error": abs(pred - o.known),
            "null_error": o.null_error,
            "blind_spot": fp.blind_spot,
        })

    return {
        "rows": rows,
        "mae": float(np.mean([r["error"] for r in rows])) if rows else float("nan"),
        "null_mae": float(np.mean([r["null_error"] for r in rows])) if rows else float("nan"),
        "median_ratio": float(np.median([r["ratio"] for r in rows])) if rows else float("nan"),
    }


def main() -> int:
    r = run_ceiling()
    print("ORACLE CEILING — every arm gets its own potency, fitted to its own answer\n")
    print(f"{'arm':<22} {'known':>8} {'potency':>8} {'predicted':>11} {'error':>7} {'null err':>9}")
    print("-" * 82)
    for f in r["fits"]:
        print(f.line())
    print("-" * 82)

    print(f"\nREACHABLE ARMS ({r['n_reachable']}) — the honest headline")
    print(f"  oracle MAE {r['mae_reachable']:.1f}pp   vs null {r['null_mae_reachable']:.1f}pp")
    reach_beats = r["mae_reachable"] < r["null_mae_reachable"]
    print(f"  the ceiling {'BEATS' if reach_beats else 'DOES NOT BEAT'} the null")

    if r["n_unreachable"]:
        print(f"\nUNREACHABLE ARMS ({r['n_unreachable']}) — sign errors, not calibration errors")
        print(f"  oracle MAE {r['mae_unreachable']:.1f}pp   vs null "
              f"{r['null_mae_unreachable']:.1f}pp")
        print("  No potency reproduces these: the dial cannot produce a benefit at any")
        print("  strength (bricks/qsp_velez.py). Their error is the model's direction")
        print("  failure, and no potency source of any quality touches it.")

    print(f"\nCOMBINED ({len(r['fits'])}) — shown last because it mixes the two")
    print(f"  oracle MAE {r['mae_all']:.1f}pp   vs null {r['null_mae_all']:.1f}pp")

    print("\n  One free parameter fitted to one number reproduces that number, so the")
    print("  reachable-arm figure above is near-tautological on its own. What it")
    print("  establishes is only that the response curve PASSES THROUGH each of those")
    print("  trials' answers at some potency. Whether an independent source finds that")
    print("  potency is the next table, and it is the one that decides anything.")

    rec = recoverability()
    print("\n\nRECOVERABILITY — do the INDEPENDENT (MRI-fitted) potencies land there?\n")
    if not rec["rows"]:
        print("  No arm has both an oracle potency and an in-range MRI-fitted potency.")
    else:
        print(f"{'arm':<22} {'oracle s':>9} {'MRI s':>8} {'ratio':>7} "
              f"{'pred at MRI s':>14} {'error':>7} {'null err':>9}")
        print("-" * 82)
        for row in rec["rows"]:
            print(f"{row['arm']:<22} {row['oracle_potency']:>9.2f} {row['mri_potency']:>8.2f} "
                  f"{row['ratio']:>7.2f} {row['predicted_at_mri']:>+13.1f}% "
                  f"{row['error']:>7.1f} {row['null_error']:>9.1f}")
        print("-" * 82)
        print(f"\n  MRI-fitted potencies score {rec['mae']:.1f}pp against a "
              f"{rec['null_mae']:.1f}pp null on {len(rec['rows'])} arms.")
        print(f"  Median MRI:oracle potency ratio {rec['median_ratio']:.2f}.")
        rec_beats = rec["mae"] < rec["null_mae"]
        print(f"  The independent channel {'BEATS' if rec_beats else 'DOES NOT BEAT'} "
              "the null where it can be scored at all.")

    print("\nWHAT THIS DECIDES")
    if reach_beats:
        print("  On arms the model can reach, its response curve passes through the")
        print("  trial's answer, so nothing about the curve's SHAPE forbids a correct")
        print("  result there. Half the exam is not reachable at all, and no potency")
        print("  source touches that half.")
        print("  Whether calibration is worth pursuing turns on the recoverability")
        print("  table above, not on the ceiling: a channel that cannot find the")
        print("  oracle's potency cannot deliver the ceiling's error.")
    else:
        print("  The model FORM is the blocker, not the calibration. Even a potency")
        print("  fitted directly to the answer -- which no independent source can beat --")
        print("  loses to answering with the other arms' average. BUILD_PLAN blockers (4)")
        print("  and (6) are closed routes on this model, not deferred ones: more or")
        print("  better potency data cannot lift a ceiling that is already below the null.")
    print("\n  Trial numbers are real and cited (docs/TRIAL_ANCHORS.md). Simulation")
    print("  numbers are proxies from a ported toy model, and the Sormani map is applied")
    print("  outside what it was fitted on. Nothing here is evidence about multiple")
    print("  sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
