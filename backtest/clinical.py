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

Trial anchors (via PubMed):
  IFN-β   : PRISMS Study Group, Lancet 1998; PMID 9820297. 27–33% relapse-rate
            reduction vs placebo (mean relapses 2.56 placebo → 1.73–1.82).  IMPROVES.
  Glatiramer: Copolymer 1 / Johnson KP et al., Neurology 1995; PMID 11902590.
            Reduces relapse rate (~29%).  IMPROVES.
  APL CGP77116: Bielekova B et al., Nat Med 2000; PMID 11017150; doi:10.1038/80516.
            Phase II halted; 3 patients had MS exacerbations, 2 drug-linked
            (encephalitogenic MBP 83-99).  HARMS.

Run:  python -m backtest.clinical
"""

from __future__ import annotations

import copy
import statistics as stats
from dataclasses import dataclass

from spine.pipeline import Pipeline
from spine.run_demo import build_stages
from bricks.vpop import sample_vpop


@dataclass(frozen=True)
class ClinicalOutcome:
    arm: str
    direction: str                 # "improves" | "harms" | "neutral"
    relapse_change_pct: float | None  # vs placebo; negative = fewer relapses. None = unquantified
    source: str


# The known ground truth. This list IS the exam the stack must pass.
KNOWN_OUTCOMES = [
    ClinicalOutcome("untreated", "neutral", 0.0, "control arm"),
    ClinicalOutcome("IFN-beta", "improves", -30.0,
                    "PRISMS, Lancet 1998 (PMID 9820297): 27–33% relapse-rate reduction"),
    ClinicalOutcome("glatiramer acetate", "improves", -29.0,
                    "Copolymer 1 / Johnson 1995, Neurology (PMID 11902590)"),
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


def _direction(change_pct: float) -> str:
    if change_pct < -NEUTRAL_BAND:
        return "improves"
    if change_pct > NEUTRAL_BAND:
        return "harms"
    return "neutral"


def run_gate(n: int = 12, seed: int = 1) -> dict:
    cohort = sample_vpop(n=n, seed=seed)
    relapse = {o.arm: _arm_relapse(o.arm, cohort) for o in KNOWN_OUTCOMES}
    base = relapse["untreated"]

    rows, dir_pass, mag_pass, mag_total = [], 0, 0, 0
    for o in KNOWN_OUTCOMES:
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
    print(f"{'arm':<20} {'sim Δrelapse':>12} {'predicted':>10} {'known':>9}  {'DIR':>4}")
    print("-" * 64)
    for o, change, pred, d_ok, m_ok in g["rows"]:
        mark = "PASS" if d_ok else "FAIL"
        print(f"{o.arm:<20} {change:>+11.0f}% {pred:>10} {o.direction:>9}  {mark:>4}")
    print("-" * 64)
    print(f"DIRECTION gate: {g['dir_pass']}/{g['n']} arms correct   "
          f"(magnitude: {g['mag_pass']}/{g['mag_total']} within {MAG_TOLERANCE:.0f}pp)")
    passed = g["dir_pass"] == g["n"]
    print(f"\nDIRECTION GATE: {'PASS' if passed else 'INCOMPLETE'} ({g['dir_pass']}/{g['n']} arms)")
    if passed:
        print("  The stack can now represent BOTH benefit and harm and directs every known")
        print("  arm correctly. This is a CAPABILITY MILESTONE, not validation:")
        print("    - intervention params now come from a 2-PARAMETER MECHANISM-CLASS RULE")
        print("      (bricks/grounding.py), keyed on each drug's independent in-vitro")
        print("      mechanism — NOT hand-tuned per arm, NOT fit to any relapse number.")
        print("      The gate now tests the RULE, not four hand-set values.")
        print("    - but the rule is coarse and reasoned, not data-fit; with only 4 arms a")
        print("      rigorous leave-one-out (fit strengths on N-1, predict the Nth) needs")
        print("      more arms than we have data for. Necessary, not yet sufficient.")
    else:
        fails = [o.arm for o, _, _, d_ok, _ in g["rows"] if not d_ok]
        print(f"  not yet reproduced: {fails}")
    print("\n  Trial outcomes are REAL and cited (see module docstring / RESULTS.md).")
    print("  Pipeline numbers are proxies from toy models. Nothing here is evidence about MS.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
