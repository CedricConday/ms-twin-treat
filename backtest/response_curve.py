"""The simulation's dose-response curve: treat -> simulated relapse change.

The clinical gate runs the whole pipeline per arm, which costs ~17s per twelve
patients. Every question worth asking of it -- what would an arm at treat=0.55
predict? what strength best explains a set of arms? -- is a question about the
SAME curve evaluated at different points, so the curve is computed once on a
grid and cached.

This is not an approximation of the model; it is the model, tabulated. The
pipeline is deterministic given (cohort seed, treat, immunogenic), so a grid
point is the exact value the gate would produce for an arm with those
parameters, and anything between grid points is linear interpolation.

Interpolation error, measured rather than assumed: at treat=0.78 the gate itself
produces -77.06% and this table interpolates -76.26%, i.e. **0.8pp** in the
steepest part of the curve (0.75 -> -71.3%, 0.80 -> -79.6%). That is small
against the tens of percentage points the clinical gate and the LOO are
measuring, but it is not zero, and a finer grid is the fix if it ever matters.

Build or refresh the cache:

    PYTHONPATH=. python -m backtest.response_curve

It writes results/response_curve.json. Nothing else in the repo writes that
file; delete it and it is rebuilt from scratch.
"""

from __future__ import annotations

import copy
import json
import statistics as stats
from pathlib import Path

import numpy as np

from bricks.intervention import Intervention
from bricks.vpop import sample_vpop
from spine.pipeline import Pipeline
from spine.run_demo import build_stages

CACHE = Path(__file__).resolve().parent.parent / "results" / "response_curve.json"

# The grid. 0.05 steps over the whole [0,1] range of `treat`, at the two harm
# levels the mechanism rule can produce (0.0, and IMMUNOGENIC_STRENGTH for a
# regulation-disrupting drug).
TREAT_GRID = [round(x, 2) for x in np.arange(0.0, 1.0001, 0.05)]


def _relapse_for(treat: float, immunogenic: float, cohort: list[dict]) -> float:
    """Mean relapse proxy for a hypothetical arm with these parameters."""
    arm = Intervention("grid", treat=treat, immunogenic=immunogenic,
                       mechanism="grid", notes="response-curve grid point")
    stages = build_stages(with_data=False, arm="untreated")
    stages[0].intervention = arm          # B7 is the intervention stage
    pipe = Pipeline(stages, name=f"grid t={treat} i={immunogenic}")
    members = [dict(m) for m in copy.deepcopy(cohort)]
    for m in members:
        m.pop("intervention_name", None)  # do not let a per-patient override win
    results = pipe.run_cohort(members, verbose=False)
    return stats.mean(r["readout"]["relapse_proxy"] for r in results)


def build(n: int = 12, seed: int = 1, immunogenic_levels: tuple[float, ...] = (0.0, 0.4)) -> dict:
    cohort = sample_vpop(n=n, seed=seed)
    base = _relapse_for(0.0, 0.0, cohort)
    curve: dict[str, dict[str, float]] = {}
    for immuno in immunogenic_levels:
        points = {}
        for treat in TREAT_GRID:
            relapse = _relapse_for(treat, immuno, cohort)
            points[f"{treat:.2f}"] = (relapse - base) / base * 100 if base else 0.0
        curve[f"{immuno:.2f}"] = points
    return {"n": n, "seed": seed, "untreated_relapse": base, "curve": curve}


def load(path: Path = CACHE) -> dict:
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Build it with: PYTHONPATH=. python -m backtest.response_curve")
    return json.loads(path.read_text())


def predict(treat: float, immunogenic: float = 0.0, table: dict | None = None) -> float:
    """Simulated relapse change (%) for an arm at this strength, by interpolation."""
    table = table or load()
    key = f"{immunogenic:.2f}"
    if key not in table["curve"]:
        raise KeyError(f"no grid at immunogenic={key}; have {sorted(table['curve'])}")
    points = table["curve"][key]
    xs = np.array(sorted(float(k) for k in points))
    ys = np.array([points[f"{x:.2f}"] for x in xs])
    return float(np.interp(treat, xs, ys))


def main() -> int:
    table = build()
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    CACHE.write_text(json.dumps(table, indent=2) + "\n")
    print(f"wrote {CACHE}")
    for immuno, points in table["curve"].items():
        shown = {k: round(v, 1) for k, v in points.items() if k in ("0.00", "0.25", "0.50", "0.75", "1.00")}
        print(f"  immunogenic={immuno}: {shown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
