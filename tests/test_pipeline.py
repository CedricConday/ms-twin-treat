"""Invariant tests for ms-twin-treat.

These lock in the behaviors the project's credibility rests on: the ruler
discriminates, the nulls behave, the toy bricks are directionally coherent, the
barrier is consequential, the plausibility filter actually rejects, and the
whole pipeline composes. They do NOT test that any model is biologically correct
(it isn't — see RESULTS.md §6). They test that the honest machinery works.

Run:  python -m pytest tests/ -q      (repo root, venv active)
"""

from __future__ import annotations

import numpy as np
import pytest

from backtest.harness import PerturbationBenchmark, score_delta
from bricks.baselines import GlobalMeanShiftNull, IdentityNull


# --------------------------------------------------------------------------- #
# harness + nulls
# --------------------------------------------------------------------------- #

def _toy_bench(seed=0, n_genes=200, n_types=5):
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n_genes)]
    shared = np.zeros(n_genes); shared[rng.choice(n_genes, 30, replace=False)] = rng.normal(2, .3, 30)
    ctrl, pert = {}, {}
    for c in range(n_types):
        base = rng.normal(5, 1, n_genes)
        spec = np.zeros(n_genes); spec[rng.choice(n_genes, 15, replace=False)] = rng.normal(1.5, .3, 15)
        ctrl[f"ct{c}"] = base
        pert[f"ct{c}"] = base + shared + spec
    return PerturbationBenchmark(genes, ctrl, pert)


def test_identity_null_scores_zero():
    """Predicting no change captures none of a real perturbation."""
    bench = _toy_bench()
    r = bench.evaluate(IdentityNull())["delta_pearson"].mean()
    assert abs(r) < 1e-6


def test_meanshift_beats_identity():
    bench = _toy_bench()
    r_id = bench.evaluate(IdentityNull())["delta_pearson"].mean()
    r_ms = bench.evaluate(GlobalMeanShiftNull(bench.global_mean_delta()))["delta_pearson"].mean()
    assert r_ms > r_id + 0.05


def test_score_delta_perfect_and_null():
    true = np.array([1.0, -2.0, 3.0, 0.0])
    assert score_delta(true, true)["delta_pearson"] == pytest.approx(1.0, abs=1e-9)
    assert score_delta(true, np.zeros_like(true))["delta_pearson"] == 0.0  # constant -> 0, not NaN


def test_score_delta_shape_mismatch_raises():
    with pytest.raises(ValueError):
        score_delta(np.zeros(3), np.zeros(4))


# --------------------------------------------------------------------------- #
# bricks are directionally coherent (toy, but not backwards)
# --------------------------------------------------------------------------- #

def test_qsp_treatment_preserves_more_myelin():
    from bricks.qsp import simulate
    assert simulate(treat=0.6)["M"][-1] > simulate(treat=0.0)["M"][-1]


def test_qsp_does_not_diverge():
    """The logistic cap must keep the untreated trajectory finite (regression)."""
    from bricks.qsp import simulate
    assert np.isfinite(simulate(treat=0.0)["M"][-1])


def test_abm_treatment_reduces_damage():
    from bricks.abm import simulate
    assert simulate(treat=0.8, seed=0)[-1] < simulate(treat=0.0, seed=0)[-1]


def test_immunogenic_intervention_causes_harm():
    """The harm mechanism: an immunogenic therapy makes things WORSE than untreated
    (more demyelination) — the opposite of a suppressive one. Without this the stack
    cannot reproduce a therapy that harmed patients (e.g. APL CGP77116)."""
    from bricks.qsp import simulate as qsp
    from bricks.abm import simulate as abm
    assert qsp(treat=0.0, immuno=0.4)["M"][-1] < qsp(treat=0.0, immuno=0.0)["M"][-1]
    assert abm(treat=0.0, immuno=0.4, seed=0)[-1] > abm(treat=0.0, immuno=0.0, seed=0)[-1]


def test_barrier_peripheral_full_exposure():
    """A drug that need not cross the barrier is not penalized."""
    from bricks.barrier import BarrierStage
    s = BarrierStage().run({"intervention": {"dose": 1.0, "cns_required": False}, "bbb_disruption": 0.3})
    assert s["cns_exposure"]["effective"] == 1.0


def test_barrier_cns_required_is_gated():
    from bricks.barrier import BarrierStage
    s = BarrierStage().run({"intervention": {"dose": 1.0, "cns_required": True}, "bbb_disruption": 0.1})
    assert 0.0 <= s["cns_exposure"]["effective"] < 1.0


def test_readout_barrier_makes_undelivered_drug_worse():
    """Same sim + same drug, worse delivery -> worse clinical proxy (B6 is not decoration)."""
    from bricks.readout import ReadoutStage
    dmg = np.linspace(0, 0.3, 20)
    peripheral = ReadoutStage().run({"abm_damage": dmg, "cns_exposure": {"effective": 1.0},
                                     "intervention": {"treat": 0.5}})["readout"]
    locked = ReadoutStage().run({"abm_damage": dmg, "cns_exposure": {"effective": 0.08},
                                 "intervention": {"treat": 0.5}})["readout"]
    assert locked["lesion_proxy"] > peripheral["lesion_proxy"]


# --------------------------------------------------------------------------- #
# the wedge: plausibility filter must actually reject
# --------------------------------------------------------------------------- #

def test_vpop_filter_discriminates():
    """The plausibility filter must REJECT implausible params and ACCEPT diseased
    ones -- deterministic, not a stochastic acceptance rate. A filter that accepts
    everything is decoration; this proves it does not."""
    from bricks.vpop import _is_plausible
    assert not _is_plausible(0.10, 0.05)   # no feedback/damage -> no disease -> reject
    assert _is_plausible(0.90, 0.70)        # disease develops -> accept


def test_vpop_prevalence_weighting_matches_target():
    """MAPEL step: after weighting, the weighted mass per severity bin matches the
    target prevalence -- the correction the raw plausible set does not have."""
    from bricks.vpop import sample_vpop, weight_to_prevalence, DEFAULT_PREVALENCE
    cohort = sample_vpop(n=40, seed=2)
    weight_to_prevalence(cohort)
    n = len(cohort)
    present = {p["severity_bin"] for p in cohort}
    for name, target in DEFAULT_PREVALENCE.items():
        if name not in present:
            continue  # a bin absent from the draw can't be matched into existence
        weighted_mass = sum(p["weight"] for p in cohort if p["severity_bin"] == name) / n
        assert abs(weighted_mass - target) < 0.02


def test_vpop_patients_have_pipeline_keys():
    from bricks.vpop import sample_vpop
    for m in sample_vpop(n=5, seed=1):
        assert {"patient_id", "seed", "bbb_disruption", "qsp_params"} <= set(m)
        assert m["vpop_meta"]["validated"] is False


# --------------------------------------------------------------------------- #
# the whole spine composes
# --------------------------------------------------------------------------- #

def test_pipeline_end_to_end_produces_readout():
    from spine.pipeline import Pipeline
    from spine.run_demo import build_stages
    from bricks.vpop import sample_vpop
    stages = build_stages(with_data=False, arm="IFN-beta")
    result = Pipeline(stages).run(sample_vpop(n=1, seed=0)[0], verbose=False)
    assert "readout" in result and "lesion_proxy" in result["readout"]
    assert result["readout_meta"]["validated"] is False  # never claims validation
