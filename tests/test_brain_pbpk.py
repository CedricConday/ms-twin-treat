"""Tests for the Verscheijden 2019 brain PBPK port.

Two of these guard mistakes that were actually made while writing the module,
which is the reason they are worth keeping: the constant-input boundary
condition that cannot express penetration, and the efflux-vs-perfusion mix-up in
the antibody conversion. Both produced plausible-looking numbers.
"""

from __future__ import annotations

import numpy as np
import pytest

from bricks.brain_pbpk import (
    BRAIN_BLOOD_FRACTION,
    DEFAULT_ADULT,
    PSB_PER_BRAIN_L,
    PSC_PER_BRAIN_L,
    PSE_FLAT,
    VCCSF_L,
    BrainPBPK,
    antibody_psb,
    simulate,
)

# --------------------------------------------------------------------------- #
# transcribed constants — each is a citation to pcbi.1007117.s001.R
# --------------------------------------------------------------------------- #

def test_ps_products_match_the_published_r_source():
    # L565/L569: PSb = Vbrain/(1.36/1.04)*1.875 ; PSc = Vbrain/(1.36/1.04)*0.9375
    assert PSB_PER_BRAIN_L == pytest.approx(1.875 / (1.36 / 1.04))
    assert PSC_PER_BRAIN_L == pytest.approx(0.9375 / (1.36 / 1.04))
    assert PSE_FLAT == 300.0          # L567
    assert PSB_PER_BRAIN_L == pytest.approx(2 * PSC_PER_BRAIN_L)


def test_volumes_match_the_published_r_source():
    assert VCCSF_L == 0.143           # L187
    assert BRAIN_BLOOD_FRACTION == 0.05   # L186


def test_the_four_compartments_partition_the_brain():
    m = DEFAULT_ADULT
    total = m.v_bb + m.v_bm + m.v_ccsf + m.v_scsf
    assert total == pytest.approx(m.v_brain)
    assert m.v_scsf / (m.v_ccsf + m.v_scsf) == pytest.approx(0.2)  # L189 cap


def test_a_brain_too_small_for_its_csf_is_refused():
    tiny = BrainPBPK(v_brain=0.1, q_brain=42.0, q_bulk=0.01, q_sin=0.02,
                     q_sout=0.01, q_csink=0.01, q_ssink=0.01)
    with pytest.raises(ValueError, match="too small"):
        _ = tiny.v_bm


def test_adult_substitution_is_flagged_not_hidden():
    """The paper is pediatric; the adult volumes and flows are this repo's."""
    assert DEFAULT_ADULT.meta["pediatric_source"] is True
    assert DEFAULT_ADULT.meta["adult_values_substituted"] is True
    assert "pediatric" in DEFAULT_ADULT.source


# --------------------------------------------------------------------------- #
# the boundary-condition trap
# --------------------------------------------------------------------------- #

def test_constant_input_equilibrates_regardless_of_permeability():
    """The mistake this module was first written with.

    Brain mass has no elimination but CSF turnover, so under a CONSTANT arterial
    concentration every compartment approaches that concentration whatever PSb
    is. Permeability sets the rate, not the level. Reporting that as "CNS
    penetration" gave an antibody an 80% brain:plasma ratio.
    """
    permeable = simulate(DEFAULT_ADULT, hours=240.0)
    ab = BrainPBPK(**{**DEFAULT_ADULT.__dict__,
                      "ps_b": antibody_psb(0.0015, DEFAULT_ADULT.q_bulk)})
    barely = simulate(ab, hours=240.0)

    assert permeable["equilibrium_fraction"] > 0.9
    # Both climb; the near-impermeable one is slower, not lower-plateauing.
    assert barely["equilibrium_fraction"] > simulate(ab, hours=24.0)["equilibrium_fraction"]


def test_auc_ratio_is_the_penetration_quantity_and_it_is_ps_dependent():
    ab_psb = antibody_psb(0.0015, DEFAULT_ADULT.q_bulk)
    ab = BrainPBPK(**{**DEFAULT_ADULT.__dict__, "ps_b": ab_psb, "ps_c": ab_psb / 2})
    small = simulate(DEFAULT_ADULT, hours=10.0, plasma_half_life_h=2.5)
    large = simulate(ab, hours=1920.0, plasma_half_life_h=480.0)
    assert small["auc_ratio"] > 100 * large["auc_ratio"]


def test_constant_input_reports_no_meaningful_auc_ratio_label():
    """`equilibrium_fraction` and `auc_ratio` must stay distinct keys.

    Collapsing them into one "effective" number is how the wrong one gets read.
    """
    r = simulate(DEFAULT_ADULT, hours=48.0)
    assert "equilibrium_fraction" in r
    assert "auc_ratio" in r
    assert r["plasma_half_life_h"] is None


# --------------------------------------------------------------------------- #
# the efflux-vs-perfusion trap, and the check it enables
# --------------------------------------------------------------------------- #

def test_antibody_conversion_competes_against_efflux_not_perfusion():
    """The second mistake made here, wrong by a factor of ~4,000.

    Entry is PSb*Cbb, the only exit is q_efflux*Cbm, so
    Cbm/Cbb = PSb/(PSb + q_efflux). Using perfusion (~42 L/h) in place of CSF
    turnover (~0.0105 L/h) returns a PSb that puts an antibody at 80% of plasma.
    """
    r = 0.0015
    psb = antibody_psb(r, DEFAULT_ADULT.q_bulk)
    assert psb / (psb + DEFAULT_ADULT.q_bulk) == pytest.approx(r, rel=1e-9)
    # and the wrong version would be far larger
    assert psb < 0.001 * antibody_psb(r, DEFAULT_ADULT.q_brain)


def test_pardridge_ratio_is_reproduced_rather_than_assumed():
    """An independent check: the ported structure lands in the published band.

    Pardridge 2019 puts therapeutic antibody CNS:serum at ~0.1-0.2%. That number
    is an INPUT to the PSb conversion, but the AUC ratio that comes back out is
    produced by the four-compartment model with its own CSF turnover — so
    landing in the band is a consistency check on the port, not a tautology.
    """
    psb = antibody_psb(0.0015, DEFAULT_ADULT.q_bulk)
    ab = BrainPBPK(**{**DEFAULT_ADULT.__dict__, "ps_b": psb, "ps_c": psb / 2})
    auc = simulate(ab, hours=1920.0, plasma_half_life_h=480.0)["auc_ratio"]
    assert 0.0005 < auc < 0.003, f"expected the 0.1-0.2% band, got {auc:.5f}"


def test_a_ratio_outside_zero_one_is_refused():
    for bad in (0.0, 1.0, -0.1, 2.0):
        with pytest.raises(ValueError):
            antibody_psb(bad, DEFAULT_ADULT.q_bulk)


def test_zero_efflux_is_refused_with_the_reason():
    with pytest.raises(ValueError, match="no route out of the brain"):
        antibody_psb(0.0015, 0.0)


# --------------------------------------------------------------------------- #
# integration sanity
# --------------------------------------------------------------------------- #

def test_trajectories_are_finite_and_non_negative():
    r = simulate(DEFAULT_ADULT, hours=72.0)
    for key in ("brain_blood", "brain_mass", "cranial_csf", "spinal_csf"):
        arr = np.asarray(r[key])
        assert np.all(np.isfinite(arr))
        assert np.all(arr >= -1e-9)


def test_nothing_is_marked_validated():
    assert simulate(DEFAULT_ADULT)["validated"] is False
