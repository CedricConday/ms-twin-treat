"""End-to-end demo: a virtual cohort through every scale of the spine.

The wiring was built BEFORE the bricks, on purpose (PARALLEL_PLAN §3). Each
brick that lands REPLACES its pass-through, so this script runs at every commit
rather than only at the end. What follows is therefore an honest picture of how
much of the pipeline is real right now — read the summary at the bottom.

Run:
    python spine/run_demo.py                    # 3 patients, no dataset needed
    python spine/run_demo.py --n 10             # bigger cohort
    python spine/run_demo.py --arm untreated    # switch treatment arm
    python spine/run_demo.py --with-data        # also run B2/B3 (loads Kang, slow)

What this shows: that intervention -> cell -> network -> QSP -> population ->
barrier -> clinical composes as one state, and where each scale plugs in.

What this does NOT show: that any scale is correct. Every brick is a toy or a
stand-in. See BUILD_PLAN §5 — built is not validated.
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

# real bricks — no dataset required
from bricks.abm import ABMBrick            # noqa: E402  B5
from bricks.barrier import BarrierStage    # noqa: E402  B6
from bricks.intervention import (          # noqa: E402  B7
    LIBRARY,
    InterventionStage,
)
from bricks.qsp import QSPBrick            # noqa: E402  B4
from bricks.readout import ReadoutStage    # noqa: E402  B8
from bricks.vpop import sample_vpop        # noqa: E402  B9


def build_stages(with_data: bool, arm: str) -> list:
    """Assemble the pipeline. Data-dependent bricks are opt-in."""
    stages: list = [InterventionStage(arm)]                       # B7  REAL

    if with_data:
        # B2/B3 need the Kang expression matrix in memory.
        from backtest.harness import PerturbationBenchmark  # noqa: F401
        from bricks.cell_transfer import CellTransferModel
        from bricks.grn import GRNBrick
        from data.kang import to_benchmark

        bench = to_benchmark()
        stages.append(_CellStage(CellTransferModel(bench)))       # B2  REAL
        stages.append(GRNBrick())                                 # B3  REAL
    else:
        stages.append(PassThroughStage(
            "B2 cell perturbation", "cell_delta", value={},
            requires=("intervention",),
            reason="skipped: needs the Kang matrix, run with --with-data"))
        stages.append(PassThroughStage(
            "B3 gene regulatory net", "grn_edges", value=[],
            requires=("cell_delta",),
            reason="skipped: needs the Kang matrix, run with --with-data"))

    stages += [
        QSPBrick(),                                               # B4  REAL (toy)
        ABMBrick(),                                               # B5  REAL (toy)
        BarrierStage(),                                           # B6  REAL (toy)
        ReadoutStage(),                                           # B8  REAL (toy, barrier-gated)
    ]
    return stages


class _CellStage:
    """Adapts the harness-facing PerturbationModel to a spine Stage.

    B2 implements `predict(control_mean, cell_type)` because that is what the
    backtest harness scores. The spine wants `run(state)`. Rather than give the
    brick two personalities, the adapter lives here.
    """

    name = "B2 cell perturbation"
    requires = ("intervention",)

    def __init__(self, model) -> None:
        self.model = model

    def run(self, state: MultiScaleState) -> MultiScaleState:
        bench = getattr(self.model, "benchmark", None)
        deltas = {}
        if bench is not None:
            # NOTE: cell_types is a property, not a method (harness.py) — no parens.
            for ct in bench.cell_types:
                prof = bench.profiles.get(ct) if hasattr(bench, "profiles") else None
                if prof is None:
                    continue
                # store the delta (perturbed - control), which is what cell_delta means
                deltas[ct] = self.model.predict(prof.control_mean, ct) - prof.control_mean
        state["cell_delta"] = deltas
        state["cell_meta"] = {"validated": False, "model": self.model.name,
                              "note": "beats no null yet; bar is mean-shift 0.85"}
        return state

    __call__ = run


def make_cohort(n: int, arm: str | None = None) -> list[MultiScaleState]:
    """The virtual population — now the REAL B9 (bricks/vpop.py).

    B9 Latin-hypercube samples physiological parameter sets and rejects the
    implausible region via a plausibility filter (Allen-Rieger flavour), then
    hands the accepted cohort to the spine. Method real, filter model a toy
    (validated=False). This replaces the old seed-jitter stand-in.
    """
    return sample_vpop(n=n, seed=0, arm=arm)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=3, help="cohort size")
    ap.add_argument("--arm", default="IFN-beta", choices=sorted(LIBRARY))
    ap.add_argument("--with-data", action="store_true",
                    help="also run B2/B3 (loads the Kang dataset)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    stages = build_stages(args.with_data, args.arm)
    pipe = Pipeline(stages, name="ms-twin-treat v0")
    cohort = make_cohort(args.n, args.arm)

    print(f"MS-Twin spine — {len(pipe.stages)} stages, cohort of {len(cohort)}, arm={args.arm!r}")
    print("=" * 76)

    results = []
    for i, member in enumerate(cohort, 1):
        print(f"\npatient {member['patient_id']}  ({i}/{len(cohort)})  "
              f"bbb_disruption={member['bbb_disruption']}")
        results.append(pipe.run(member, verbose=not args.quiet))

    final = results[-1]
    print("\n" + "=" * 76)
    print(pipe.report(final))

    keys = [k for k in final if not k.endswith("_meta")]
    n_stand = sum(1 for k in keys if pipe._unvalidated(final, k))
    print("\n" + "=" * 76)
    print(f"BUILT, NOT VALIDATED — {n_stand} of {len(keys)} state keys are unvalidated.")
    print("Every brick is a toy or a stand-in. The scales compose; that is the claim.")
    print("Nothing here is evidence about multiple sclerosis.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
