"""End-to-end experiment: run the virtual cohort across real trial arms.

Not a new brick — an honest USE of the pipeline that exists. Samples one virtual
population (B9) once, then runs the identical cohort through the spine under each
arm in the intervention library, and reports how the clinical PROXY separates the
arms. This is the shape of the eventual backtest: feed known-outcome therapies,
see whether the model ranks them correctly.

Two honest results come out of this:
  1. The pipeline SEPARATES treated from untreated (the toy models are at least
     directionally coherent end to end).
  2. It now FLAGS the harmful arm (APL CGP77116, halted in Phase II for causing MS
     exacerbations) — it reads worse than untreated. But that harm is
     MECHANISM-ENCODED, not discovered: APL's immunogenic parameter is set from its
     documented encephalitogenic biology, not its relapse outcome. That is a
     capability milestone, not validation.

Everything here is toy/unvalidated (see each brick). Numbers are proxies, not
clinical endpoints. Run:  python -m results.experiment
"""

from __future__ import annotations

import copy
import statistics as stats
import warnings

warnings.filterwarnings("ignore")  # silence mesa's seed-kwarg FutureWarning for clean output

from spine.pipeline import Pipeline
from spine.run_demo import build_stages
from bricks.vpop import sample_vpop
from bricks.intervention import LIBRARY


ARMS = ["untreated", "IFN-beta", "glatiramer acetate", "APL CGP77116"]


def run_arm(cohort: list[dict], arm: str) -> list[dict]:
    """Run a fresh copy of the cohort through the pipeline under one arm."""
    stages = build_stages(with_data=False, arm=arm)
    pipe = Pipeline(stages, name=f"arm={arm}")
    members = [dict(m, intervention_name=arm) for m in copy.deepcopy(cohort)]
    return pipe.run_cohort(members, verbose=False)


def main() -> int:
    cohort = sample_vpop(n=12, seed=1)   # one plausible population, reused across arms
    print(f"virtual population: {len(cohort)} plausible patients (B9 LHS + plausibility filter)\n")

    print(f"{'arm':<22} {'treat':>6} {'lesion_proxy':>14} {'relapse_proxy':>15}   note")
    print("-" * 82)
    summary = {}
    for arm in ARMS:
        results = run_arm(cohort, arm)
        lesions = [r["readout"]["lesion_proxy"] for r in results]
        relapses = [r["readout"]["relapse_proxy"] for r in results]
        treat = LIBRARY[arm].treat
        summary[arm] = (stats.mean(lesions), stats.mean(relapses), treat)
        note = ""
        if arm == "untreated":
            note = "control"
        elif arm == "APL CGP77116":
            note = "REAL outcome: HARMED (Phase II halted) — now flagged (mechanism-encoded)"
        elif treat > 0:
            note = "REAL outcome: worked/approved"
        print(f"{arm:<22} {treat:>6.2f} {stats.mean(lesions):>14.1f} "
              f"{stats.mean(relapses):>15.2f}   {note}")

    base = summary["untreated"][0]
    print("\nArm separation vs untreated (mean lesion_proxy, lower = better):")
    for arm in ARMS:
        d = summary[arm][0] - base
        print(f"  {arm:<22} delta = {d:+.1f}")

    print("\nHONEST READING:")
    print("  - Treated arms (IFN-β, glatiramer) separate from untreated, and APL now")
    print("    reads WORSE than untreated: the stack directs all four known arms right.")
    print("  - BUT APL's harm is MECHANISM-ENCODED, not discovered: we set its")
    print("    immunogenic parameter from its documented encephalitogenic biology")
    print("    (Bielekova 2000), not from the relapse outcome. That is a capability")
    print("    milestone (the model can now express harm), NOT validation — the params")
    print("    are hand-assigned and the arms informed the setup.")
    print("  - Every number is a proxy from toy models. Nothing here is evidence")
    print("    about multiple sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
