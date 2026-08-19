"""B9 — Virtual population generator (THE WEDGE).

Per ms-twin/docs/RESEARCH_FINDINGS.md: no open Python implementation of the
rigorous plausible-patient method exists on GitHub, and nobody has run it on a
neuroinflammation model. This is a first, honest v1 of exactly that.

Method (Allen-Rieger flavour of plausible-patient generation):
  1. Latin-hypercube sample physiological parameter sets over plausibility bounds.
  2. PLAUSIBILITY FILTER: simulate each candidate (untreated) through the toy QSP
     and ACCEPT only patients whose outcome falls in a physiologically plausible
     window -- here, they must actually develop disease (some demyelination) but
     not collapse instantly. Implausible parameter combinations are rejected.
  3. Return the accepted virtual cohort as initial states for the spine.

This is the real method; the model it filters against is a TOY, so the resulting
cohort is method-real, biology-illustrative (validated=False). MAPEL prevalence-
weighting (reweight the accepted cohort to a target biomarker prevalence) is the
documented next step, not yet implemented -- see TODO below.

What makes this NOT the make_cohort() stand-in it replaces: that only jittered a
seed. This samples a parameter space and rejects the implausible region, which
is the entire point of plausible-patient generation.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import qmc

from bricks.qsp import simulate as qsp_simulate

# Physiological parameters sampled per patient, with plausibility bounds.
# (name, low, high). r_CA/k_dmg feed the QSP; bbb_disruption feeds the barrier.
PARAM_BOUNDS: list[tuple[str, float, float]] = [
    ("r_CA", 0.10, 1.80),          # immune-cytokine feedback strength
    ("k_dmg", 0.05, 1.00),         # myelin damage rate
    ("bbb_disruption", 0.00, 0.80),  # blood-brain-barrier leakiness
]

# Plausibility window on untreated final demyelination (1 - myelin).
# A plausible MS-relevant virtual patient shows disease but is not instantly
# fully demyelinated. Values are illustrative, not fitted. The bounds above are
# deliberately WIDE enough to sample implausible regions (low r_CA/k_dmg -> no
# disease develops; extreme combos -> instant collapse), so this filter actually
# REJECTS -- which is the entire point of plausible-patient generation. A filter
# that accepts everything is decoration.
PLAUSIBLE_DAMAGE = (0.20, 0.95)


def _is_plausible(r_CA: float, k_dmg: float) -> bool:
    traj = qsp_simulate(params={"r_CA": r_CA, "k_dmg": k_dmg}, treat=0.0)
    final_damage = 1.0 - float(traj["M"][-1])
    lo, hi = PLAUSIBLE_DAMAGE
    return bool(np.isfinite(final_damage) and lo <= final_damage <= hi)


def sample_vpop(n: int = 20, seed: int = 0, arm: str | None = None,
                oversample: int = 8, max_rounds: int = 20) -> list[dict]:
    """Generate n plausible virtual patients. Returns spine-ready initial states."""
    names = [b[0] for b in PARAM_BOUNDS]
    l_bounds = [b[1] for b in PARAM_BOUNDS]
    u_bounds = [b[2] for b in PARAM_BOUNDS]
    sampler = qmc.LatinHypercube(d=len(PARAM_BOUNDS), seed=seed)

    accepted: list[dict] = []
    n_tried = 0
    for _ in range(max_rounds):
        raw = sampler.random(max(n * oversample, n))
        cand = qmc.scale(raw, l_bounds, u_bounds)
        for row in cand:
            n_tried += 1
            vals = dict(zip(names, (float(x) for x in row)))
            if not _is_plausible(vals["r_CA"], vals["k_dmg"]):
                continue
            i = len(accepted)
            accepted.append({
                "patient_id": f"vp{i:03d}",
                "seed": i,
                "bbb_disruption": round(vals["bbb_disruption"], 3),
                "qsp_params": {"r_CA": round(vals["r_CA"], 3),
                               "k_dmg": round(vals["k_dmg"], 3)},
                **({"intervention_name": arm} if arm else {}),
                "vpop_meta": {
                    "validated": False,
                    "method": "LHS + plausibility filter (Allen-Rieger flavour)",
                    "note": "method real; plausibility model is a TOY QSP. "
                            "MAPEL prevalence-weighting = TODO.",
                },
            })
            if len(accepted) >= n:
                accepted[0]["vpop_meta"]["acceptance_rate"] = round(len(accepted) / n_tried, 3)
                return accepted
    # ran out of rounds; return what we have, labelled
    if accepted:
        accepted[0]["vpop_meta"]["acceptance_rate"] = round(len(accepted) / max(n_tried, 1), 3)
    return accepted


class VPopSampler:
    """Convenience wrapper so a caller can hold config and sample repeatedly."""

    name = "B9 virtual-population (LHS + plausibility filter)"

    def __init__(self, seed: int = 0) -> None:
        self.seed = seed

    def sample(self, n: int, arm: str | None = None) -> list[dict]:
        return sample_vpop(n=n, seed=self.seed, arm=arm)


if __name__ == "__main__":
    print("B9 VPop — generating a plausible virtual cohort (the wedge)...")
    cohort = sample_vpop(n=8, seed=0, arm="IFN-beta")
    ar = cohort[0]["vpop_meta"].get("acceptance_rate")
    print(f"  accepted {len(cohort)} plausible patients  (acceptance rate ~{ar})")
    print("  sample of sampled parameter sets:")
    for m in cohort[:5]:
        p = m["qsp_params"]
        print(f"    {m['patient_id']}: r_CA={p['r_CA']}, k_dmg={p['k_dmg']}, "
              f"bbb={m['bbb_disruption']}")
    print("  method real (LHS + plausibility rejection); model is a toy; validated=False")
    print("  TODO: MAPEL prevalence-weighting to a target biomarker distribution")
