"""Which failure binds? Three candidate causes, measured one at a time.

`gate/ceiling.py` showed that per-drug potency from the MRI channel cannot lift
the gate, and named the model's FORM as what remains -- six of twelve arms are
unreachable because depletion, sequestration and transit blockade collapse onto
a dial that cannot turn downward.

Naming it is not measuring it. "The model form is the blocker" is a hypothesis
with an obvious test: if reachability were the binding constraint, then the
scorers should work on the arms that ARE reachable. This module checks that, and
separates three failures that a single headline MAE pools together:

  (1) REACHABILITY     the model cannot produce the observed direction at any
                       potency, on 6 of 12 arms
  (2) RECOVERABILITY   an independent source cannot find the potency the model
                       would need -- measured in gate/ceiling.py: 21.4pp against
                       a 10.6pp null, biased high every time
  (3) EXPRESSIVENESS   ONE shared potency cannot distinguish two drugs on the
                       same dial, whatever its value

(1) and (2) are already measured. This module measures (3), by running the
leave-one-arm-out procedure -- fit one shared potency on N-1 arms, predict the
Nth -- RESTRICTED to the reachable arms, where (1) does not apply.

WHY THE ANSWER MATTERS MORE THAN THE NUMBER
--------------------------------------------
    reachable-only LOO beats the null  ->  reachability is the binding
                                           constraint. Fixing the model's form
                                           is the whole job, and the prize is
                                           the difference between this number
                                           and the full-set one.

    reachable-only LOO loses too       ->  a third failure is present that
                                           neither a better model form nor
                                           better potency data addresses. Fixing
                                           reachability would move the headline
                                           and not the verdict, and anyone
                                           planning that work should know before
                                           starting it.

The null is the same one every scorer here uses: answer with the mean of the
training arms' outcomes, computed leave-one-out over the same restricted set so
that the rule and the null are asked about the same arm from the same evidence.

Nothing is simulated. The response table is cached.

Run:  PYTHONPATH=. python3 -m gate.headroom
"""

from __future__ import annotations

import numpy as np

from backtest.clinical import KNOWN_OUTCOMES
from backtest.lomo import _predicted, load, mechanism_groups
from gate.ceiling import run_ceiling


def _shared_potency(training: list[str], known: dict[str, float], table: dict) -> float:
    """The one potency that best explains every training arm (least squares).

    Same procedure as `backtest/lomo._fit_potency`, restated here rather than
    imported so that restricting the arm set cannot silently change what is
    being fitted.
    """
    best_s, best_sse = 0.0, float("inf")
    for s in np.arange(0.0, 0.951, 0.01):
        sse, ok = 0.0, True
        for arm in training:
            p = _predicted(arm, float(s), table)
            if np.isnan(p):
                ok = False
                break
            sse += (p - known[arm]) ** 2
        if ok and sse < best_sse:
            best_s, best_sse = float(s), sse
    return best_s


def run_headroom(table: dict | None = None) -> dict:
    table = table or load()
    known = {o.arm: o.relapse_change_pct for o in KNOWN_OUTCOMES
             if o.relapse_change_pct is not None and o.arm != "untreated"}

    ceil = run_ceiling(table)
    reachable = [f.arm for f in ceil["fits"] if f.reachable]

    # How many distinct mechanisms survive the restriction? If the reachable
    # arms are all one dial, a leave-one-ARM-out over them is easier than the
    # full test in a way that has nothing to do with reachability, and the
    # comparison would flatter itself. Reported so the reader can discount it.
    groups = mechanism_groups()
    reachable_groups = {k: [a for a in v if a in reachable]
                        for k, v in groups.items()}
    reachable_groups = {k: v for k, v in reachable_groups.items() if v}

    rows = []
    for held in reachable:
        training = [a for a in reachable if a != held]
        s = _shared_potency(training, known, table)
        pred = _predicted(held, s, table)
        null = float(np.mean([known[a] for a in training]))
        rows.append({
            "arm": held,
            "known": known[held],
            "fitted_potency": s,
            "predicted": pred,
            "error": abs(pred - known[held]),
            "null_predicted": null,
            "null_error": abs(null - known[held]),
        })

    # One arm whose mechanism group has no other member is held out with no
    # same-dial training data at all, which is a leave-one-MECHANISM-out in
    # disguise and can dominate a six-arm headline. Reported as a sensitivity so
    # the conclusion does not rest on it.
    singletons = [arms[0] for arms in reachable_groups.values() if len(arms) == 1]
    trimmed = [r for r in rows if r["arm"] not in singletons]

    return {
        "rows": rows,
        "singletons": singletons,
        "mae_excl_singletons": (float(np.mean([r["error"] for r in trimmed]))
                                if trimmed else float("nan")),
        "null_mae_excl_singletons": (float(np.mean([r["null_error"] for r in trimmed]))
                                     if trimmed else float("nan")),
        "n_reachable": len(reachable),
        "n_reachable_groups": len(reachable_groups),
        "reachable_groups": {"|".join(k): v for k, v in reachable_groups.items()},
        "mae": float(np.mean([r["error"] for r in rows])) if rows else float("nan"),
        "null_mae": float(np.mean([r["null_error"] for r in rows])) if rows else float("nan"),
        "oracle_mae": ceil["mae_reachable"],
    }


def main() -> int:
    r = run_headroom()
    print("HEADROOM — leave-one-arm-out RESTRICTED to the arms the model can reach\n")
    print("If reachability is the binding constraint, the scorer should work here.\n")
    print(f"{'held-out arm':<22} {'known':>8} {'potency':>8} {'predicted':>11} "
          f"{'error':>7} {'null err':>9}")
    print("-" * 74)
    for row in r["rows"]:
        print(f"{row['arm']:<22} {row['known']:>+7.1f}% {row['fitted_potency']:>8.2f} "
              f"{row['predicted']:>+10.1f}% {row['error']:>7.1f} {row['null_error']:>9.1f}")
    print("-" * 74)
    print(f"\n  {r['n_reachable']} reachable arms across {r['n_reachable_groups']} "
          f"distinct mechanism(s):")
    for pattern, arms in r["reachable_groups"].items():
        print(f"    {pattern:<22} {', '.join(arms)}")

    print(f"\n  restricted LOO MAE {r['mae']:.1f}pp   vs null {r['null_mae']:.1f}pp")
    print(f"  (the per-arm oracle on this same set is {r['oracle_mae']:.1f}pp — the floor)")

    beats = r["mae"] < r["null_mae"]
    print(f"\n  A shared potency {'BEATS' if beats else 'DOES NOT BEAT'} the null "
          "on reachable arms.")

    if r["singletons"]:
        beats_trim = r["mae_excl_singletons"] < r["null_mae_excl_singletons"]
        one = len(r["singletons"]) == 1
        print(f"\n  SENSITIVITY — {', '.join(r['singletons'])} "
              f"{'is the only arm on its dial' if one else 'are the only arms on their dials'}, so")
        print(f"  holding {'it' if one else 'them'} out removes every same-mechanism "
              "training arm and makes")
        print(f"  {'that' if one else 'those'} fold{'' if one else 's'} a mechanism "
              "holdout wearing an arm holdout's name. Dropping")
        print(f"  {'it' if one else 'them'}: {r['mae_excl_singletons']:.1f}pp vs null "
              f"{r['null_mae_excl_singletons']:.1f}pp — still "
              f"{'BEATS' if beats_trim else 'DOES NOT BEAT'}. The conclusion does not")
        print("  rest on that fold.")

    print("\nWHAT THIS DECIDES")
    if beats:
        print("  Reachability is the binding constraint. On arms the model can reach, the")
        print("  existing fitting procedure carries real information, so fixing the model's")
        print("  form — a CNS compartment, or any change that lets the depleting class")
        print("  produce a benefit — is the whole job. The prize is the gap between this")
        print("  number and the full-set one.")
    else:
        print("  A THIRD failure is present, and it is the one to plan around. Even where")
        print("  the model can reach the answer, one shared potency fitted on other arms")
        print("  does not predict the held-out arm better than their average. Fixing")
        print("  reachability would move the headline MAE and not the verdict, because a")
        print("  single strength cannot express how two drugs on the same dial differ —")
        print("  and gate/ceiling.py has already measured that the per-drug potency which")
        print("  would express it is not recoverable from the independent channel.")
        print("  Both routes out of this are therefore measured shut, which is a stronger")
        print("  statement than either measurement alone.")
    print("\n  Trial numbers are real and cited (docs/TRIAL_ANCHORS.md). Simulation")
    print("  numbers are proxies from a ported toy model. Nothing here is evidence")
    print("  about multiple sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
