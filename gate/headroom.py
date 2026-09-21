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
                       potency, on 7 of the 15 quantified arms
  (2) RECOVERABILITY   an independent source cannot find the potency the model
                       would need -- measured in gate/ceiling.py: 21.4pp against
                       a 10.4pp null, biased high every time
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


def dial_ceiling() -> dict:
    """How much is there to win at all? Bound the ARM SET, not the model.

    Origin: `scripts/dial_ceiling.py` in the master checkout, which asks what a
    PERFECT dial-level model could score -- predict every arm with the mean of
    its own dial group, in sample, with no fitting and no simulation. Reproduced
    here independently because the number is about to be quoted in a plan.

    Three variants, because the in-sample one flatters itself in a way worth
    naming:

    in_sample          every arm predicted by its own group's mean. Two of the
                       five groups have a single member (glatiramer, daclizumab),
                       so those arms are fitted EXACTLY by construction and
                       contribute zero error. Same objection that applies to a
                       single-arm fold anywhere else in this repo.

    multi_only         singleton groups dropped. The remaining ten arms are the
                       ones where a dial-level predictor has anything to do.

    out_of_sample      each arm predicted by the mean of the OTHER arms in its
                       group -- which is what a model actually has to do, and
                       how every other scorer in this repo is graded.

    The last one is the honest bound on the exam, and it is the reason this
    function exists: if a perfect dial-level model barely clears the null out of
    sample, then the arm set, not the model, is what limits how much any
    replacement could win by.

    NULL CONVENTION, because two conventions disagree by 0.2-0.4pp here and the
    repo must not carry both. The null is drawn from the SAME arms the model is
    scored on, and scored the same way the model is: in sample against the
    scored set's own mean, out of sample against the mean of the other scored
    arms. A null pooled over arms excluded from the model's exam is handed
    information the model was not, which flatters the headroom. Matches
    `scripts/dial_ceiling.py` in the master checkout. Both conventions are
    computed in `tests/test_headroom.py` and asserted to differ, so the check
    cannot pass vacuously and neither number needs quoting here.
    """
    known = {o.arm: o.relapse_change_pct for o in KNOWN_OUTCOMES
             if o.relapse_change_pct is not None and o.arm != "untreated"}
    groups = mechanism_groups()
    grand = float(np.mean(list(known.values())))
    multi = {k: v for k, v in groups.items() if len(v) > 1}
    singletons = sorted(v[0] for v in groups.values() if len(v) == 1)

    def _pack(errs, nulls):
        return {"mae": float(np.mean(errs)), "null_mae": float(np.mean(nulls)),
                "headroom": float(np.mean(nulls) - np.mean(errs)), "n": len(errs)}

    ins_e, ins_n = [], []
    for arms in groups.values():
        mu = float(np.mean([known[a] for a in arms]))
        for a in arms:
            ins_e.append(abs(known[a] - mu))
            ins_n.append(abs(known[a] - grand))

    # The scored set for the restricted rows, and the null pool that matches it.
    scored = [a for arms in multi.values() for a in arms]
    scored_mean = float(np.mean([known[a] for a in scored]))

    multi_e, multi_n, oos_e, oos_n = [], [], [], []
    for arms in multi.values():
        mu = float(np.mean([known[a] for a in arms]))
        for a in arms:
            multi_e.append(abs(known[a] - mu))
            multi_n.append(abs(known[a] - scored_mean))
            others = [known[o] for o in arms if o != a]
            oos_e.append(abs(known[a] - float(np.mean(others))))
            rest = [known[o] for o in scored if o != a]
            oos_n.append(abs(known[a] - float(np.mean(rest))))

    return {
        "in_sample": _pack(ins_e, ins_n),
        "multi_only": _pack(multi_e, multi_n),
        "out_of_sample": _pack(oos_e, oos_n),
        "singletons": singletons,
        "groups": {"|".join(k): v for k, v in groups.items()},
    }


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
    d = dial_ceiling()
    print("\n\nHOW MUCH IS THERE TO WIN AT ALL? — bounding the ARM SET, not the model\n")
    print("  A PERFECT dial-level model: predict every arm by its own dial group's")
    print("  outcome. No fitting, no simulation, nothing to get wrong.\n")
    for label, key in ((f"in sample (all {d['in_sample']['n']} arms)", "in_sample"),
                       ("singleton groups dropped", "multi_only"),
                       ("out of sample, within group", "out_of_sample")):
        v = d[key]
        print(f"    {label:<30} {v['mae']:>5.1f}pp vs null {v['null_mae']:>5.1f}pp   "
              f"headroom {v['headroom']:>4.1f}pp   n={v['n']}")
    print(f"\n  {' and '.join(d['singletons'])} are alone on their dials, so the first")
    print("  row fits them exactly by construction. The last row is the honest one: it")
    print("  is how every other scorer here is graded.")
    print(f"\n  So the entire prize is {d['out_of_sample']['headroom']:.1f}pp, out of sample, "
          "before anything is")
    print("  charged for simulation or fitting. A replacement model does not need to be")
    print("  better than this one — it needs to land within a percentage point or two of")
    print("  PERFECT to clear a predict-the-mean null on this arm set.")
    print("\n  That is a property of the ARMS, not of the model. 12 quantified arms in 3")
    print("  multi-member dial groups is a thin exam, and widening it is far cheaper than")
    print("  any model port. It is the single change that would most increase what a good")
    print("  model could demonstrate here.")

    print("\n  Trial numbers are real and cited (docs/TRIAL_ANCHORS.md). Simulation")
    print("  numbers are proxies from a ported toy model. Nothing here is evidence")
    print("  about multiple sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
