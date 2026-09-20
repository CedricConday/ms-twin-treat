"""Tests for MRI-channel potency fitting.

The single property that makes this module worth having is that the fit never
touches a drug's relapse number. Most of what follows guards that, and the rest
guards the two things a reader would otherwise have to take on trust: that every
magnitude carries a source, and that the one available cross-check disagrees.
"""

from __future__ import annotations

import numpy as np
import pytest

from backtest import potency
from backtest.lomo import POTENCY_GRID
from bricks.profiles import PROFILES, touched_points
from bricks.sormani import invert, predict_relapse_ratio


def _table(arm: str, ratios: dict[float, float]) -> dict:
    """A response table where `arm` produces the given lesion ratios.

    Stored the way the real table stores them — as post-Sormani relapse changes
    — so the inversion under test is actually exercised.
    """
    key = "|".join(touched_points(PROFILES[arm]))
    return {
        "potency_grid": POTENCY_GRID,
        "table": {key: {
            f"{s:g}": {"mean": predict_relapse_ratio(ratios[s]).percent_change,
                       "sd": 0.0, "n_ok": 1}
            for s in POTENCY_GRID
        }},
    }


# --------------------------------------------------------------------------- #
# the non-circularity property
# --------------------------------------------------------------------------- #

def test_fitting_never_reads_a_relapse_outcome(monkeypatch):
    """The whole point. If KNOWN_OUTCOMES is touched during a fit, the gate is
    no longer testing anything — it would be scoring a number it helped choose.
    """
    import backtest.clinical as clinical

    class Tripwire(list):
        def __iter__(self):
            raise AssertionError(
                "potency fitting read KNOWN_OUTCOMES — that is the circularity "
                "this module exists to avoid")

    monkeypatch.setattr(clinical, "KNOWN_OUTCOMES", Tripwire())
    ratios = {s: max(1.0 - s, 0.02) for s in POTENCY_GRID}
    fit = potency.fit_potency("ocrelizumab", 0.5, "test", _table("ocrelizumab", ratios))
    assert fit.potency >= 0.0


def test_the_only_trial_input_is_a_lesion_ratio():
    """`fit_potency`'s signature is the guarantee: a ratio and a source, no ARR."""
    import inspect
    params = set(inspect.signature(potency.fit_potency).parameters)
    assert params == {"arm", "observed_lesion_ratio", "source", "table"}


# --------------------------------------------------------------------------- #
# the fit itself
# --------------------------------------------------------------------------- #

def test_fit_recovers_the_potency_that_reproduces_an_observed_ratio():
    ratios = {s: float(np.exp(-3.0 * s)) for s in POTENCY_GRID}
    target = float(np.exp(-3.0 * 0.4))
    fit = potency.fit_potency("ocrelizumab", target, "test", _table("ocrelizumab", ratios))
    assert fit.potency == pytest.approx(0.40, abs=0.05)
    assert fit.achieved_lesion_ratio == pytest.approx(target, rel=0.15)


def test_multipliers_keep_the_sourced_direction():
    """Ocrelizumab's ke goes DOWN. A fitted magnitude may not flip that."""
    ratios = {s: max(1.0 - s, 0.02) for s in POTENCY_GRID}
    fit = potency.fit_potency("ocrelizumab", 0.3, "test", _table("ocrelizumab", ratios))
    assert fit.points == ("ke",)
    assert fit.multipliers["ke"] < 1.0


def test_a_ratio_the_model_cannot_reach_is_flagged_not_clamped():
    """Reporting the nearest grid point silently would be a fabricated fit."""
    ratios = {s: max(1.0 - 0.5 * s, 0.6) for s in POTENCY_GRID}
    fit = potency.fit_potency("ocrelizumab", 0.01, "test", _table("ocrelizumab", ratios))
    assert fit.in_range is False


def test_a_small_lesion_effect_inherits_sormanis_blind_spot():
    """The lenercept case: no MRI signal is indistinguishable from no drug."""
    ratios = {s: max(1.0 - s, 0.02) for s in POTENCY_GRID}
    fit = potency.fit_potency("ocrelizumab", 0.98, "test", _table("ocrelizumab", ratios))
    assert fit.blind_spot is True


def test_a_large_lesion_effect_is_not_flagged():
    ratios = {s: max(1.0 - s, 0.02) for s in POTENCY_GRID}
    fit = potency.fit_potency("ocrelizumab", 0.2, "test", _table("ocrelizumab", ratios))
    assert fit.blind_spot is False


# --------------------------------------------------------------------------- #
# sourcing
# --------------------------------------------------------------------------- #

def test_a_magnitude_without_a_source_is_refused():
    ratios = {s: max(1.0 - s, 0.02) for s in POTENCY_GRID}
    with pytest.raises(ValueError, match="needs a source"):
        potency.fit_potency("ocrelizumab", 0.3, "   ", _table("ocrelizumab", ratios))


def test_every_observed_ratio_names_its_trial_and_an_identifier():
    """A PMID or an NCT number — the registry results are a primary source too.

    Most of these now come from ClinicalTrials.gov posted results rather than
    journal text, which is why NCT counts. What must never be acceptable is a
    ratio with no identifier at all.
    """
    for arm, (ratio, source) in potency.OBSERVED_LESION_RATIOS.items():
        assert 0.0 < ratio < 1.5, f"{arm}: implausible lesion ratio {ratio}"
        assert ("PMID" in source) or ("NCT" in source), (
            f"{arm}: source must carry a PMID or an NCT identifier")


def test_every_observed_ratio_names_its_mri_metric():
    """new-T2, Gd-T1/scan, CUAL/year and active-T2 are different measurements.

    Two arms measured on different metrics must not be compared as if they were
    the same quantity, so each ratio records which one it is.
    """
    known = ("[new-T2]", "[Gd-T1/scan]", "[CUAL/year]", "[active-T2]")
    for arm, (_, source) in potency.OBSERVED_LESION_RATIOS.items():
        assert any(k in source for k in known), (
            f"{arm}: source must name its MRI metric, one of {known}")


def test_out_of_range_arms_are_exactly_the_ones_on_damage_raising_dials():
    """Not a fitting failure — a restatement of the model's structural limit.

    gamma_E and a lowered alpha_R both RAISE damage, so no potency reproduces a
    lesion reduction on them. If an arm on those dials ever fits in range, the
    model changed and bricks/qsp_velez.py's docstring is stale.
    """
    fits = potency.fit_all()
    out = {a for a, f in fits.items() if not f.in_range}
    for arm in out:
        points = touched_points(PROFILES[arm])
        raises_damage = ("gamma_E" in points
                         or (PROFILES[arm].alpha_R < 1.0 and points == ("alpha_R",)))
        assert raises_damage, f"{arm} is out of range for an unexplained reason"


def test_the_extraction_gap_is_countable_and_complete():
    """Every quantified arm is either fitted or explicitly listed as pending.

    An arm that is silently in neither list is the failure mode this guards:
    a gap nobody can see is a gap nobody closes.
    """
    from backtest.lomo import mechanism_groups
    quantified = {a for arms in mechanism_groups().values() for a in arms}
    covered = set(potency.OBSERVED_LESION_RATIOS) | set(potency.PENDING_EXTRACTION)
    assert quantified <= covered
    assert not (set(potency.OBSERVED_LESION_RATIOS) & set(potency.PENDING_EXTRACTION))


def test_non_positive_ratio_is_rejected():
    with pytest.raises(ValueError):
        potency.fit_potency("ocrelizumab", 0.0, "test")


def test_unknown_arm_is_rejected():
    with pytest.raises(KeyError):
        potency.fit_potency("aspirin", 0.5, "test")


# --------------------------------------------------------------------------- #
# the cross-check, asserted because it disagrees
# --------------------------------------------------------------------------- #

def test_the_two_independent_sources_for_ocrelizumab_disagree():
    """Pinned so nobody quietly reconciles them.

    ke = 0.85 from Martinez-Pasamar's EAE fit, ke ~ 0.40 from OPERA's MRI
    ratio. The gap is the best estimate this repo has of how far the potency
    layer can be trusted, and it is more informative than either number.
    """
    fits = potency.fit_all()
    mri_ke = fits["ocrelizumab"].multipliers["ke"]
    eae_ke = PROFILES["ocrelizumab"].ke
    assert eae_ke == 0.85
    assert mri_ke != eae_ke
    assert max(eae_ke, mri_ke) / min(eae_ke, mri_ke) > 1.5


def test_profiles_keeps_the_direct_measurement_as_the_default():
    """The EAE value measures the parameter; the MRI value is inverted through
    two models. Until that is argued out, the direct measurement is the default.
    """
    assert PROFILES["ocrelizumab"].ke == 0.85


# --------------------------------------------------------------------------- #
# the Sormani inverse this module depends on
# --------------------------------------------------------------------------- #

def test_sormani_round_trips():
    for rr in (0.05, 0.3, 0.9, 1.0, 1.4):
        assert invert(predict_relapse_ratio(rr).rr_relapse) == pytest.approx(rr, rel=1e-9)


def test_invert_rejects_a_non_positive_ratio():
    with pytest.raises(ValueError):
        invert(0.0)
