"""Tests for plausible-patient generation against the GROUNDED QSP.

`sample_vpop` filters candidates through `bricks/qsp.py`, the toy. Once the
pipeline switched to the Vélez port that left the filter validating patients
against a model nothing runs, and `VelezQSPBrick` silently dropped the toy
parameters it produced. `sample_vpop_velez` is the same Allen-Rieger method
against the model actually in use.

Two properties are worth guarding: the bounds are the paper's own, and the
filter genuinely rejects. A plausibility filter that accepts everything is
decoration, which the module docstring says and these tests enforce.
"""

from __future__ import annotations

import pytest

from bricks.qsp_velez import VELEZ_PARAMS, VelezQSPBrick
from bricks.vpop import (
    VELEZ_PARAM_BOUNDS,
    _is_plausible_velez,
    sample_vpop,
    sample_vpop_velez,
)


def test_bounds_are_the_papers_own_sweep_ranges():
    """Vélez de Mendizábal 2011 Table 1 gives these as intervals, not points.

    "Maximum Teff proliferation rate [1:2] day-1"
    "Maximum Treg proliferation and activation rate [0.25:2] day-1"

    Sampling inside them uses the authors' stated plausible region rather than
    a range invented here.
    """
    bounds = dict((name, (lo, hi)) for name, lo, hi in VELEZ_PARAM_BOUNDS)
    assert bounds["alpha_E"] == (1.0, 2.0)
    assert bounds["alpha_R"] == (0.25, 2.0)
    for name in bounds:
        assert name in VELEZ_PARAMS, f"{name} is not a parameter of the ported model"


def test_the_papers_operating_point_sits_inside_the_bounds():
    for name, lo, hi in VELEZ_PARAM_BOUNDS:
        assert lo <= VELEZ_PARAMS[name] <= hi


def test_sampled_parameters_reach_the_model_instead_of_being_dropped():
    """The defect this function exists to fix.

    A toy-sampled cohort hands VelezQSPBrick names it does not have, so the
    per-patient variation vanishes. A Vélez-sampled cohort must be honoured.
    """
    brick = VelezQSPBrick(t_end=365.0)
    velez_patient = sample_vpop_velez(n=1, seed=1)[0]
    assert brick.run(dict(velez_patient))["qsp_traj"]["ignored_params"] == []

    toy_patient = sample_vpop(n=1, seed=1)[0]
    dropped = brick.run(dict(toy_patient))["qsp_traj"]["ignored_params"]
    assert dropped == ["k_dmg", "r_CA"], (
        "the toy cohort's parameters should still be visibly dropped, not silently used")


def test_the_filter_actually_rejects():
    """A filter that accepts everything is decoration.

    Recorded in `vpop_meta["n_tried"]`, so the rejection rate is inspectable
    rather than assumed.
    """
    cohort = sample_vpop_velez(n=8, seed=1)
    assert len(cohort) == 8
    assert cohort[-1]["vpop_meta"]["n_tried"] > len(cohort)


def test_a_healthy_configuration_is_rejected():
    """High alpha_R is the paper's HEALTHY regime and produces almost no damage.

    That is precisely the implausible region for an MS virtual patient, so the
    filter must cut it. If this ever passes, the plausibility window has been
    widened until it stopped filtering.
    """
    assert _is_plausible_velez(alpha_E=1.5, alpha_R=2.0, seed=1) is False


def test_an_autoimmune_configuration_is_accepted():
    assert _is_plausible_velez(alpha_E=1.5, alpha_R=0.3, seed=1) is True


def test_the_accepted_cohort_lands_in_the_autoimmune_regime():
    """Not asserted as a design choice — it is what the filter selects.

    The paper identifies low alpha_R as the autoimmune configuration. The filter
    is told only "must develop disease and stay in regime", and it recovers that
    region on its own, which is a small independent check that the port behaves
    the way the paper describes.
    """
    cohort = sample_vpop_velez(n=8, seed=1)
    alpha_Rs = [p["qsp_params"]["alpha_R"] for p in cohort]
    assert max(alpha_Rs) < 1.0, f"expected the autoimmune regime, got {alpha_Rs}"


def test_patients_are_distinct_and_seeded_from_the_cohort_seed():
    """The bug fixed in sample_vpop must not be reintroduced here.

    A cohort seed that does not reach the per-patient simulation seed produces
    identical cohorts, which is how the clinical gate ended up unable to see its
    own virtual population at all.
    """
    a = sample_vpop_velez(n=6, seed=1)
    b = sample_vpop_velez(n=6, seed=2)
    seeds_a = [p["seed"] for p in a]
    seeds_b = [p["seed"] for p in b]
    assert len(set(seeds_a)) == len(seeds_a)
    assert not set(seeds_a) & set(seeds_b)
    assert [p["qsp_params"] for p in a] != [p["qsp_params"] for p in b]


def test_patients_actually_differ_in_simulated_damage():
    brick = VelezQSPBrick(t_end=730.0)
    damages = [float(brick.run(dict(p))["qsp_traj"]["total_damage"][-1])
               for p in sample_vpop_velez(n=6, seed=1)]
    assert len(set(round(d, 6) for d in damages)) > 1


def test_metadata_names_the_model_and_stays_unvalidated():
    meta = sample_vpop_velez(n=2, seed=1)[0]["vpop_meta"]
    assert meta["validated"] is False
    assert meta["model"] == "velez2011"
    assert "Allen-Rieger" in meta["method"]


def test_an_impossible_request_raises_rather_than_returning_a_short_cohort():
    """Silently returning fewer patients than asked for would quietly shrink a
    cohort, which is exactly the class of failure this session kept finding."""
    with pytest.raises(RuntimeError, match="plausible patients"):
        sample_vpop_velez(n=4, seed=1, oversample=1, max_rounds=1)
