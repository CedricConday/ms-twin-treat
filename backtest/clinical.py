"""Clinical backtest gate — does the pipeline reproduce KNOWN trial outcomes?

This gate defines "viable." Before any prediction about the unknown is trusted,
the setup must reproduce known clinical history IN BOTH DIRECTIONS:
  - known SUCCESSES must come out as improvements  (IFN-β, glatiramer)
  - known FAILURES/HARMS must come out as harm      (APL CGP77116)

A pass requires the correct DIRECTION for every arm; magnitude is a secondary,
harder bar. The trial outcomes below are REAL and cited (via PubMed). The
pipeline's numbers are proxies from toy models, so this gate is EXPECTED TO BE
RED today — and that is the point: it turns "viable" into a number to watch turn
green as toy bricks are replaced with data-grounded ones.

The arms are in KNOWN_OUTCOMES below, each with the trial, the comparator that
trial used, and its PMID; docs/TRIAL_ANCHORS.md carries the same table alongside
the mechanism-class assignment for every drug, which is made from pharmacology
and never from the outcome.

Some arms are scored against an ACTIVE COMPARATOR rather than placebo (OPERA vs
IFN beta-1a, CARE-MS vs IFN beta-1a, OPTIMUM vs teriflunomide). Those are scored
the way the trial ran: the simulated drug arm against the simulated comparator
arm, not against untreated.

Run:  python -m backtest.clinical
"""

from __future__ import annotations

import copy
import statistics as stats
from dataclasses import dataclass

from bricks.intervention import LIBRARY
from bricks.vpop import sample_vpop
from spine.pipeline import Pipeline
from spine.run_demo import build_stages


@dataclass(frozen=True)
class ClinicalOutcome:
    arm: str
    direction: str                 # "improves" | "harms" | "neutral"
    relapse_change_pct: float | None  # vs the comparator; negative = fewer relapses. None = unquantified
    source: str
    comparator: str = "untreated"  # the arm the trial actually measured against


# The known ground truth. This list IS the exam the stack must pass.
#
# Every number is the trial's own, against the comparator the trial used --
# active-comparator studies are NOT converted to a "vs placebo" figure, because
# chaining two trials' effects would invent precision that does not exist. Full
# table with PMIDs and the mechanism-class assignments: docs/TRIAL_ANCHORS.md.
#
# Two arms here are deliberate counterexamples: lenercept and atacicept are
# immunosuppressive by mechanism and HARMED patients. The mechanism-class rule in
# bricks/grounding.py predicts benefit for both, so the gate is expected to fail
# them. That is the point of scoring more than four arms.
KNOWN_OUTCOMES = [
    ClinicalOutcome("untreated", "neutral", 0.0, "control arm"),
    ClinicalOutcome("IFN-beta", "improves", -30.0,
                    "PRISMS, Lancet 1998 (PMID 9820297): 27-33% relapse-rate reduction"),
    ClinicalOutcome("glatiramer acetate", "improves", -29.0,
                    "Copolymer 1 / Johnson 1995, Neurology (PMID 11902590); CONFIRM 2012 "
                    "(PMID 22992072) ARR 0.29 vs 0.40 placebo, 29% relative reduction"),
    ClinicalOutcome("natalizumab", "improves", -68.0,
                    "AFFIRM, NEJM 2006 (PMID 16510744): relapse rate at 1 year 68% lower"),
    ClinicalOutcome("fingolimod", "improves", -55.0,
                    "FREEDOMS, NEJM 2010 (PMID 20089952): ARR 0.18 (0.5 mg) vs 0.40 placebo"),
    ClinicalOutcome("teriflunomide", "improves", -31.5,
                    "TEMSO, NEJM 2011 (PMID 21991951): ARR 0.37 (14 mg) vs 0.54 placebo"),
    ClinicalOutcome("dimethyl fumarate", "improves", -53.0,
                    "DEFINE, NEJM 2012 (PMID 22992073): ARR 0.17 (BID) vs 0.36 placebo, "
                    "53% relative reduction"),
    ClinicalOutcome("ocrelizumab", "improves", -46.0,
                    "OPERA I, NEJM 2017 (PMID 28002679): ARR 0.16 vs 0.29 on IFN beta-1a",
                    comparator="IFN-beta"),
    ClinicalOutcome("alemtuzumab", "improves", -55.0,
                    "CARE-MS I, Lancet 2012 (PMID 23122652): rate ratio 0.45 vs IFN beta-1a, "
                    "54.9% improvement",
                    comparator="IFN-beta"),
    ClinicalOutcome("ponesimod", "improves", -30.5,
                    "OPTIMUM, JAMA Neurol 2021 (PMID 33779698): ARR 0.202 vs 0.290 on "
                    "teriflunomide",
                    comparator="teriflunomide"),
    ClinicalOutcome("lenercept", "harms", None,
                    "Lenercept MS Study Group, Neurology 1999 (PMID 10449104): more patients "
                    "with exacerbations, occurring earlier (p=0.006)"),
    ClinicalOutcome("atacicept", "harms", None,
                    "ATAMS, Lancet Neurol 2014 (PMID 24613349): halted early; ARR 0.86 / 0.79 / "
                    "0.98 (25 / 75 / 150 mg) vs 0.38 placebo. Magnitude left unscored: wide CIs, "
                    "non-monotonic in dose, trial stopped early"),
    ClinicalOutcome("IFN-gamma", "harms", None,
                    "Panitch, Lancet 1987 (PMID 2882294): 7 of 18 patients had exacerbations on "
                    "treatment, above both the pre- and post-treatment rate"),
    ClinicalOutcome("APL CGP77116", "harms", None,
                    "Bielekova, Nat Med 2000 (PMID 11017150, doi:10.1038/80516): halted; exacerbations"),
]

NEUTRAL_BAND = 5.0   # |change| < 5% reads as "no effect"
MAG_TOLERANCE = 15.0  # magnitude within 15 percentage points counts as a hit


def _arm_relapse(arm: str, cohort: list[dict]) -> float:
    stages = build_stages(with_data=False, arm=arm)
    pipe = Pipeline(stages, name=f"arm={arm}")
    members = [dict(m, intervention_name=arm) for m in copy.deepcopy(cohort)]
    results = pipe.run_cohort(members, verbose=False)
    return stats.mean(r["readout"]["relapse_proxy"] for r in results)


def _arm_key(arm: str) -> tuple:
    """What the simulation actually depends on: the intervention's parameters.

    Two arms with identical parameters produce identical runs, so they are run
    once. Under today's class-level rule that collapses every suppressive drug
    onto one simulation -- which is not an optimization artefact, it is the
    limitation blocker (4) exists to fix, made visible. When per-drug strengths
    land, the keys separate on their own.
    """
    i = LIBRARY[arm]
    return (i.treat, i.immunogenic, i.dose, i.cns_required)


def _direction(change_pct: float) -> str:
    if change_pct < -NEUTRAL_BAND:
        return "improves"
    if change_pct > NEUTRAL_BAND:
        return "harms"
    return "neutral"


def run_gate(n: int = 12, seed: int = 1) -> dict:
    cohort = sample_vpop(n=n, seed=seed)

    arms = {o.arm for o in KNOWN_OUTCOMES} | {o.comparator for o in KNOWN_OUTCOMES}
    by_params: dict[tuple, float] = {}
    relapse: dict[str, float] = {}
    for arm in sorted(arms):
        key = _arm_key(arm)
        if key not in by_params:
            by_params[key] = _arm_relapse(arm, cohort)
        relapse[arm] = by_params[key]

    rows, dir_pass, mag_pass, mag_total = [], 0, 0, 0
    for o in KNOWN_OUTCOMES:
        base = relapse[o.comparator]
        change = (relapse[o.arm] - base) / base * 100 if base else 0.0
        pred = _direction(change)
        d_ok = pred == o.direction
        dir_pass += d_ok
        m_ok = None
        if o.relapse_change_pct is not None and o.arm != "untreated":
            mag_total += 1
            m_ok = abs(change - o.relapse_change_pct) <= MAG_TOLERANCE
            mag_pass += m_ok
        rows.append((o, change, pred, d_ok, m_ok))
    return {"rows": rows, "dir_pass": dir_pass, "n": len(KNOWN_OUTCOMES),
            "mag_pass": mag_pass, "mag_total": mag_total}


def main() -> int:
    g = run_gate()
    print("CLINICAL BACKTEST GATE — does the stack reproduce known trial history?\n")
    print(f"{'arm':<20} {'vs':<16} {'sim Δ':>7} {'known Δ':>8} {'predicted':>10} "
          f"{'known':>9}  {'DIR':>4}")
    print("-" * 82)
    for o, change, pred, d_ok, _m_ok in g["rows"]:
        mark = "PASS" if d_ok else "FAIL"
        known = f"{o.relapse_change_pct:>+7.1f}%" if o.relapse_change_pct is not None else "      —"
        print(f"{o.arm:<20} {o.comparator:<16} {change:>+6.0f}% {known:>8} {pred:>10} "
              f"{o.direction:>9}  {mark:>4}")
    print("-" * 82)
    print(f"DIRECTION gate: {g['dir_pass']}/{g['n']} arms correct   "
          f"(magnitude: {g['mag_pass']}/{g['mag_total']} within {MAG_TOLERANCE:.0f}pp)")
    passed = g["dir_pass"] == g["n"]
    print(f"\nDIRECTION GATE: {'PASS' if passed else 'INCOMPLETE'} ({g['dir_pass']}/{g['n']} arms)")
    print("    (intervention params come from the 2-PARAMETER MECHANISM-CLASS RULE in")
    print("     bricks/grounding.py, keyed on each drug's independent in-vitro mechanism —")
    print("     never hand-tuned per arm, never fit to a relapse number. The gate scores")
    print("     that rule.)")
    if passed:
        print("  Every arm is directed correctly, including the suppressive-but-harmful ones.")
        print("  That is a CAPABILITY MILESTONE, not validation: these outcomes are all known,")
        print("  and a rule that reproduces them still has to predict an arm it never saw.")
    else:
        fails = [o.arm for o, _, _, d_ok, _ in g["rows"] if not d_ok]
        print(f"  not yet reproduced: {fails}")
        print("  Expected, and the reason the arm set was grown: a rule that maps mechanism")
        print("  class to clinical direction cannot get lenercept or atacicept right — both")
        print("  are immunosuppressive by mechanism and both harmed patients — and one")
        print("  strength per class cannot separate two drugs inside a class, which is what")
        print("  the active-comparator arms ask it to do. See BUILD_PLAN.md §8 blockers (3)")
        print("  and (4); the four-arm version of this gate could not show either failure.")
    print("\n  Trial outcomes are REAL and cited (KNOWN_OUTCOMES above, docs/TRIAL_ANCHORS.md).")
    print("  Pipeline numbers are proxies from toy models. Nothing here is evidence about MS.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
