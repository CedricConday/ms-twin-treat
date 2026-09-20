"""B9 — Virtual population generator (THE WEDGE).

Per ms-twin/docs/RESEARCH_FINDINGS.md: no open Python implementation of the
rigorous plausible-patient method exists on GitHub, and nobody has run it on a
neuroinflammation model. This is a first, honest v1 of exactly that.

Method (Allen-Rieger-Musante plausible-patient generation,
**Allen RJ, Rieger TR, Musante CJ, CPT Pharmacometrics Syst Pharmacol 2016,
doi:10.1002/psp4.12063** — the citation belongs here because this brick's whole
claim is "first open implementation of that method", and a claim to implement a
published method has to name it):
  1. Latin-hypercube sample physiological parameter sets over plausibility bounds.
  2. PLAUSIBILITY FILTER: simulate each candidate (untreated) through the toy QSP
     and ACCEPT only patients whose outcome falls in a physiologically plausible
     window -- here, they must actually develop disease (some demyelination) but
     not collapse instantly. Implausible parameter combinations are rejected.
  3. Return the accepted virtual cohort as initial states for the spine.

This is the real method; the model it filters against is a TOY, so the resulting
cohort is method-real, biology-illustrative (validated=False). MAPEL prevalence-
weighting (reweight the accepted cohort to a target biomarker prevalence) is
implemented in weight_to_prevalence() below (simplified to a single severity axis;
full MAPEL optimizes multiple axes).

What makes this NOT the make_cohort() stand-in it replaces: that only jittered a
seed. This samples a parameter space and rejects the implausible region, which
is the entire point of plausible-patient generation.

STILL MISSING FROM THE STATED WEDGE (recovered from the research repo 2026-09-20)
---------------------------------------------------------------------------------
`ms-twin/docs/RESEARCH_FINDINGS.md:39` defines the differentiator as a Python
port of the prevalence/plausible-patient method **plus modern sampling —
DREAM(ZS) or simulation-based inference**. Only the first half was ever built.
Sampling here is plain Latin hypercube with a rejection filter; there is no
DREAM(ZS), no SBI, and no posterior over the parameter space at all. The
reference implementation named there, `BMSQSP/QSPToolbox`, is MATLAB +
SimBiology and has not been consulted.

That gap is not cosmetic: rejection sampling gives an accepted SET, while the
method's value is a prevalence-WEIGHTED population. `weight_to_prevalence()`
below is a simplified single-axis stand-in for MAPEL, not MAPEL.
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
            vals = dict(zip(names, (float(x) for x in row), strict=True))
            if not _is_plausible(vals["r_CA"], vals["k_dmg"]):
                continue
            i = len(accepted)
            accepted.append({
                "patient_id": f"vp{i:03d}",
                # Per-patient SIMULATION seed, and it must depend on the cohort
                # seed. It used to be plain `i`, so every cohort handed the ABM
                # the same 12 seeds no matter what `seed=` was passed. Combined
                # with two other facts -- the sampled qsp_params never reach the
                # scored path (the readout scores abm_damage), and
                # bbb_disruption is inert while no arm sets cns_required -- that
                # made the entire virtual population decorative: changing
                # `seed=` moved the clinical gate by exactly zero. Measured
                # 2026-09-20: untreated proxy 0.147137 for seed=1 and seed=2
                # alike, versus 0.152334 once the patient seeds actually differ.
                "seed": seed * 100_000 + i,
                "bbb_disruption": round(vals["bbb_disruption"], 3),
                "qsp_params": {"r_CA": round(vals["r_CA"], 3),
                               "k_dmg": round(vals["k_dmg"], 3)},
                **({"intervention_name": arm} if arm else {}),
                "vpop_meta": {
                    "validated": False,
                    "method": "LHS + plausibility filter (Allen-Rieger flavour)",
                    "note": "method real; plausibility model is a TOY QSP. "
                            "MAPEL-style prevalence weighting available via weight_to_prevalence().",
                },
            })
            if len(accepted) >= n:
                accepted[0]["vpop_meta"]["acceptance_rate"] = round(len(accepted) / n_tried, 3)
                return accepted
    # ran out of rounds; return what we have, labelled
    if accepted:
        accepted[0]["vpop_meta"]["acceptance_rate"] = round(len(accepted) / max(n_tried, 1), 3)
    return accepted


# --------------------------------------------------------------------------- #
# MAPEL-flavoured prevalence weighting — the second half of the wedge.
#
# Plausible-patient sampling (above) says WHICH virtual patients are physiologically
# admissible. It does NOT say how COMMON each kind should be. Real cohorts have a
# prevalence structure (more mild than severe disease); an unweighted plausible set
# does not. MAPEL (Schmidt et al.) fixes this by assigning each plausible patient a
# weight so the weighted ensemble matches an observed distribution.
#
# This is a simplified, honest v1: bin patients by disease severity (untreated final
# demyelination), then weight each bin by target/observed density ratio so the
# weighted cohort matches a target prevalence. Full MAPEL optimises the weight
# distribution against multiple output axes; this matches one axis. The TARGET here
# is illustrative, not epidemiological — like everything else, validated=False.
# --------------------------------------------------------------------------- #

SEVERITY_BINS: list[tuple[float, float, str]] = [
    (0.20, 0.50, "mild"), (0.50, 0.75, "moderate"), (0.75, 1.01, "severe"),
]
# Illustrative target prevalence (more mild than severe). NOT fitted to any registry.
DEFAULT_PREVALENCE = {"mild": 0.50, "moderate": 0.35, "severe": 0.15}


# --------------------------------------------------------------------------- #
# The same method, against the GROUNDED model (2026-09-20)
# --------------------------------------------------------------------------- #
# `sample_vpop` above filters candidates through bricks/qsp.py, the TOY. Once
# spine/run_demo.py switched to the Velez port, that left the plausibility
# filter validating patients against a model the pipeline no longer runs -- and
# VelezQSPBrick drops `r_CA`/`k_dmg` because they do not exist in it. The method
# was still real and the coupling was gone.
#
# This is the same Allen-Rieger procedure against the grounded model. What makes
# it better grounded than a re-tuning: the sampling bounds are the PAPER'S OWN
# ranges. Velez de Mendizabal 2011 Table 1 gives alpha_E and alpha_R as
# intervals rather than points -- `[1:2]` and `[0.25:2]` day^-1 -- because those
# are the ranges the authors swept. Sampling inside them is using the paper's
# stated plausible region, not inventing one.

VELEZ_PARAM_BOUNDS: list[tuple[str, float, float]] = [
    ("alpha_E", 1.0, 2.0),    # Table 1: "Maximum Teff proliferation rate [1:2] day-1"
    ("alpha_R", 0.25, 2.0),   # Table 1: "Maximum Treg proliferation and activation
                              #           rate [0.25:2] day-1"
]

# Plausibility window on untreated total damage at the horizon. A patient must
# develop disease and must not leave the model's regime. The bounds are on the
# MEDIAN-scale damage this model produces (see backtest/lomo.py on how wide that
# distribution is) and are illustrative, not fitted -- but the REJECTION is real:
# alpha_R near the top of its range gives a healthy configuration that produces
# almost no damage, which is exactly the implausible region a filter should cut.
VELEZ_PLAUSIBLE_DAMAGE = (0.05, 1000.0)
VELEZ_HORIZON_DAYS = 730.0


def _is_plausible_velez(alpha_E: float, alpha_R: float, seed: int) -> bool:
    """Accept a candidate only if it develops disease and stays in regime."""
    from bricks.qsp_velez import UNTREATED_PROFILE, simulate

    traj = simulate(UNTREATED_PROFILE,
                    params={"alpha_E": alpha_E, "alpha_R": alpha_R},
                    t_end=VELEZ_HORIZON_DAYS, seed=seed)
    if not traj["in_regime"]:
        return False
    damage = float(traj["total_damage"][-1])
    lo, hi = VELEZ_PLAUSIBLE_DAMAGE
    return bool(np.isfinite(damage) and lo <= damage <= hi)


def sample_vpop_velez(n: int = 20, seed: int = 0, arm: str | None = None,
                      oversample: int = 6, max_rounds: int = 20) -> list[dict]:
    """Plausible virtual patients for the GROUNDED QSP.

    Same method as `sample_vpop`, different model and different parameters. The
    returned `qsp_params` use Velez names, so `VelezQSPBrick` honours them
    instead of dropping them -- which is the whole point.

    Note what this cohort does NOT yet feed: the clinical gate scores
    `abm_damage`, so these patients vary the QSP and not the number the gate
    reads. Wiring the readout onto `qsp_damage` is the remaining step
    (BUILD_PLAN §8.4, "wiring").
    """
    names = [b[0] for b in VELEZ_PARAM_BOUNDS]
    l_bounds = [b[1] for b in VELEZ_PARAM_BOUNDS]
    u_bounds = [b[2] for b in VELEZ_PARAM_BOUNDS]
    sampler = qmc.LatinHypercube(d=len(VELEZ_PARAM_BOUNDS), seed=seed)

    accepted: list[dict] = []
    n_tried = 0
    for _ in range(max_rounds):
        raw = sampler.random(max(n * oversample, n))
        cand = qmc.scale(raw, l_bounds, u_bounds)
        for row in cand:
            n_tried += 1
            vals = dict(zip(names, (float(x) for x in row), strict=True))
            patient_seed = seed * 100_000 + len(accepted)
            if not _is_plausible_velez(vals["alpha_E"], vals["alpha_R"], patient_seed):
                continue
            i = len(accepted)
            accepted.append({
                "patient_id": f"vv{i:03d}",
                "seed": patient_seed,
                "qsp_params": {"alpha_E": round(vals["alpha_E"], 3),
                               "alpha_R": round(vals["alpha_R"], 3)},
                **({"intervention_name": arm} if arm else {}),
                "vpop_meta": {
                    "validated": False,
                    "method": "LHS + plausibility filter (Allen-Rieger flavour)",
                    "model": "velez2011",
                    "bounds": "Velez de Mendizabal 2011 Table 1 sweep ranges",
                    "n_tried": n_tried,
                    "note": "method real; bounds are the paper's own ranges; the "
                            "model is a transcription, not a reproduction.",
                },
            })
            if len(accepted) == n:
                return accepted
    raise RuntimeError(
        f"only {len(accepted)}/{n} plausible patients after {n_tried} candidates; "
        "widen VELEZ_PLAUSIBLE_DAMAGE or the bounds")


def _severity(patient: dict) -> float:
    """Untreated final demyelination for this patient's sampled parameters."""
    traj = qsp_simulate(params=patient.get("qsp_params", {}), treat=0.0)
    return 1.0 - float(traj["M"][-1])


def _severity_bin(dmg: float) -> str:
    for lo, hi, name in SEVERITY_BINS:
        if lo <= dmg < hi:
            return name
    return SEVERITY_BINS[-1][2] if dmg >= SEVERITY_BINS[-1][0] else SEVERITY_BINS[0][2]


def weight_to_prevalence(cohort: list[dict],
                         target: dict[str, float] | None = None) -> list[dict]:
    """Assign MAPEL-style prevalence weights so the weighted cohort matches `target`.

    Adds `weight` and `severity_bin` to each patient (weights average to 1.0). A bin
    the target does not want gets weight 0; a rare-but-wanted bin gets up-weighted.
    Mutates and returns the cohort.
    """
    target = target or DEFAULT_PREVALENCE
    bins = [_severity_bin(_severity(p)) for p in cohort]
    n = len(cohort)
    counts = {b: bins.count(b) for b in set(bins)}
    raw = []
    for b in bins:
        observed_frac = counts[b] / n if n else 0.0
        raw.append((target.get(b, 0.0) / observed_frac) if observed_frac > 0 else 0.0)
    scale = (n / sum(raw)) if sum(raw) > 0 else 1.0  # normalise mean weight to 1.0
    for p, w, b in zip(cohort, raw, bins, strict=True):
        p["weight"] = round(w * scale, 4)
        p["severity_bin"] = b
    return cohort


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

    print("\n  MAPEL-style prevalence weighting (the second half of the wedge):")
    big = sample_vpop(n=40, seed=2)
    weight_to_prevalence(big)  # match DEFAULT_PREVALENCE (mild 0.50 / mod 0.35 / severe 0.15)
    bins = [p["severity_bin"] for p in big]
    for name in ("mild", "moderate", "severe"):
        obs = bins.count(name) / len(big)
        wt = sum(p["weight"] for p in big if p["severity_bin"] == name) / len(big)
        print(f"    {name:<9} observed {obs:>5.0%}  ->  weighted {wt:>5.0%}   (target {DEFAULT_PREVALENCE[name]:.0%})")
    print("  weighted cohort now matches the target prevalence. Target illustrative; validated=False.")
