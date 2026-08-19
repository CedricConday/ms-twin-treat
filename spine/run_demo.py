"""End-to-end demo: one virtual cohort through every scale of the spine.

This is the "built" milestone from BUILD_PLAN §7 — and it is wired BEFORE the
bricks exist, on purpose. Every stage below starts life as a labelled
pass-through. As each brick lands it REPLACES its pass-through, so this script
runs at every commit from now on rather than only at the end. A half-finished
build still demonstrates something; a half-finished build with the wiring last
demonstrates nothing.

Run:
    python spine/run_demo.py            # 3 virtual patients
    python spine/run_demo.py --n 10     # bigger cohort

What this shows: that molecular -> cell -> network -> tissue -> barrier ->
clinical composes as one state, and where each scale plugs in.

What this does NOT show: that any scale is correct. Nearly everything is
STANDIN. See the report at the bottom of the run.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spine.pipeline import (  # noqa: E402
    MultiScaleState,
    PassThroughStage,
    Pipeline,
    is_standin,
)

# --------------------------------------------------------------------------- #
# stage table — the eight bricks, in scale order.
#
# TO REPLACE A STAND-IN: import your brick and swap it in below. Keep `name`
# and the key it produces identical so the contract in PARALLEL_PLAN §4 holds.
# --------------------------------------------------------------------------- #

STAGES = [
    PassThroughStage(
        "B7 intervention", "intervention",
        value={"agent": "IFN-beta", "dose": 1.0, "on": True},
        reason="B7 not built — fixed IFN-beta on/off",
    ),
    PassThroughStage(
        "B2 cell perturbation", "cell_delta",
        value={}, requires=("intervention",),
        reason="B2 not built — no per-cell-type delta predicted yet",
    ),
    PassThroughStage(
        "B3 gene regulatory net", "grn_edges",
        value=[], requires=("cell_delta",),
        reason="B3 not built — no GRN inferred yet",
    ),
    PassThroughStage(
        "B4 QSP / ODE", "qsp_traj",
        value={"t": [], "y": []}, requires=("intervention",),
        reason="B4 not built — no inflammation ODE integrated yet",
    ),
    PassThroughStage(
        "B5 population ABM", "abm_damage",
        value=[], requires=("qsp_traj",),
        reason="B5 not built — no agent-based myelin damage yet",
    ),
    PassThroughStage(
        "B6 barrier / PBPK", "cns_exposure",
        value=0.0, requires=("intervention",),
        reason="B6 not built — no blood/CSF/CNS compartment model yet",
    ),
    PassThroughStage(
        "B8 clinical readout", "readout",
        value={"lesion_proxy": None, "relapse_proxy": None},
        requires=("abm_damage", "cns_exposure"),
        reason="B8 not built — micro-to-clinical map is open research",
    ),
]


def make_cohort(n: int) -> list[MultiScaleState]:
    """Stand-in for B9 VPop.

    A real virtual population samples parameter sets against plausibility
    bounds and weights them to a target prevalence. This just varies a seed so
    the cohort machinery is exercised — it is NOT a virtual population.
    """
    return [{"patient_id": f"vp{i:03d}", "seed": i} for i in range(n)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=3, help="cohort size")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    pipe = Pipeline(STAGES, name="ms-twin-treat v0")
    cohort = make_cohort(args.n)

    print(f"MS-Twin spine — {len(pipe.stages)} stages, cohort of {len(cohort)}")
    print("=" * 72)

    results = []
    for i, member in enumerate(cohort, 1):
        print(f"\npatient {member['patient_id']}  ({i}/{len(cohort)})")
        results.append(pipe.run(member, verbose=not args.quiet))

    print("\n" + "=" * 72)
    print(pipe.report(results[-1]))

    n_stand = sum(1 for k in results[-1] if is_standin(results[-1][k]))
    print("\n" + "=" * 72)
    if n_stand:
        print(f"BUILT, NOT VALIDATED — {n_stand} of {len(results[-1])} state keys are STANDIN.")
        print("The scales compose. Nothing here is a scientific claim.")
    else:
        print("No stand-ins left. Every key came from a real brick.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
