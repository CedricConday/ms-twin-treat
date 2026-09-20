"""Tests for the Sormani lesions→relapses map.

Protects three things:

1. The published constants, so a later edit cannot quietly retune them.
2. The boundary conditions the relation must obey to be a relation at all.
3. The blind spot, which is the one case where this map is confidently wrong
   and looks right.
"""

from __future__ import annotations

import math

import pytest

from bricks.sormani import (
    INTERCEPT,
    INTERCEPT_IS_ASSUMED,
    N_PATIENTS,
    N_TRIALS,
    R_SQUARED,
    SLOPE,
    lesion_ratio,
    predict_relapse_ratio,
)

# --------------------------------------------------------------------------- #
# published constants
# --------------------------------------------------------------------------- #

def test_published_constants_match_the_paper():
    """Sormani & Bruzzi 2013, Lancet Neurol 12(7):669-76, PMID 23743084."""
    assert SLOPE == 0.52
    assert R_SQUARED == 0.71
    assert N_TRIALS == 31
    assert N_PATIENTS == 18901


def test_intercept_is_flagged_as_assumed_not_cited():
    """The one number here that is not from the paper must say so.

    If someone reads the PDF and fills in the real intercept, this test should
    be updated deliberately — not pass silently because the flag was dropped.
    """
    assert INTERCEPT == 0.0
    assert INTERCEPT_IS_ASSUMED is True
    assert predict_relapse_ratio(0.5).intercept_assumed is True


# --------------------------------------------------------------------------- #
# the relation's own boundary conditions
# --------------------------------------------------------------------------- #

def test_no_lesion_effect_means_no_relapse_effect():
    """The condition that forces intercept = 0. If this fails the map is broken."""
    pred = predict_relapse_ratio(1.0)
    assert pred.rr_relapse == pytest.approx(1.0)
    assert pred.percent_change == pytest.approx(0.0)


def test_the_map_is_monotone_and_damped():
    """Fewer lesions must mean fewer relapses, and by LESS than proportionally.

    Slope 0.52 < 1 is the substantive content of the meta-analysis: a treatment
    buys roughly half as much relapse reduction, on the log scale, as it buys
    lesion reduction. A slope >= 1 would mean MRI understates clinical benefit,
    which is the opposite of what 31 trials found.
    """
    assert 0.0 < SLOPE < 1.0
    ratios = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.2]
    preds = [predict_relapse_ratio(r).rr_relapse for r in ratios]
    assert preds == sorted(preds)
    for r, p in zip(ratios, preds, strict=True):
        if r < 1.0:
            assert p > r, "relapse reduction must be smaller than lesion reduction"


def test_log_linearity_is_what_is_implemented():
    """The paper regresses log on log; check the implementation actually does."""
    a, b = 0.4, 0.8
    pa = predict_relapse_ratio(a).rr_relapse
    pb = predict_relapse_ratio(b).rr_relapse
    assert math.log(pa) - math.log(pb) == pytest.approx(
        SLOPE * (math.log(a) - math.log(b)))


def test_worsening_lesions_predicts_worsening_relapses():
    pred = predict_relapse_ratio(1.5)
    assert pred.rr_relapse > 1.0
    assert pred.percent_change > 0.0


# --------------------------------------------------------------------------- #
# an independent sanity check against a real trial
# --------------------------------------------------------------------------- #

def test_natalizumab_is_reproduced_within_a_few_points():
    """AFFIRM is not in this repo's fit and is a genuine out-of-sample check.

    AFFIRM (NEJM 2006, PMID 16510744) reported roughly a 90% reduction in
    gadolinium-enhancing lesions and a 68% reduction in relapse rate at one
    year. Putting the lesion effect through the map should land near the
    observed relapse effect.

    This is ONE trial and a loose check — it is not validation, and the map
    stays validated=False. It is here because a map that could not reproduce
    the best-characterised arm in the set would not be worth wiring in.
    """
    pred = predict_relapse_ratio(0.10)
    assert pred.percent_change == pytest.approx(-68.0, abs=5.0)


# --------------------------------------------------------------------------- #
# the blind spot
# --------------------------------------------------------------------------- #

def test_small_lesion_effects_are_flagged_as_the_blind_spot():
    """Where the map cannot separate 'does nothing' from 'harms'."""
    assert predict_relapse_ratio(1.0).blind_spot is True
    assert predict_relapse_ratio(0.95).blind_spot is True
    assert predict_relapse_ratio(1.05).blind_spot is True


def test_large_lesion_effects_are_not_flagged():
    assert predict_relapse_ratio(0.3).blind_spot is False
    assert predict_relapse_ratio(1.5).blind_spot is False


def test_lenercept_is_the_documented_failure_of_this_map():
    """The reason `blind_spot` exists at all.

    The lenercept trial (Neurology 1999;53:457, PMID 10449104) reported NO
    significant MRI difference while the relapse rate rose significantly
    (p=0.006) with more severe and longer relapses. Fed a null lesion effect,
    the map says 'no change' — confidently, and wrongly. The flag is what stops
    a consumer reading that as a clean negative.
    """
    pred = predict_relapse_ratio(1.0)
    assert pred.percent_change == pytest.approx(0.0)
    assert pred.blind_spot is True


def test_predictions_are_never_marked_validated():
    """The regression is published; feeding it simulated lesions is not."""
    assert predict_relapse_ratio(0.5).validated is False


# --------------------------------------------------------------------------- #
# guards
# --------------------------------------------------------------------------- #

def test_non_positive_ratio_is_rejected():
    for bad in (0.0, -1.0):
        with pytest.raises(ValueError):
            predict_relapse_ratio(bad)


def test_lesion_ratio_refuses_a_zero_comparator():
    with pytest.raises(ValueError, match="comparator arm has no lesions"):
        lesion_ratio(1.0, 0.0)


def test_lesion_ratio_refuses_a_zero_treated_arm():
    """'Perfect' is not expressible on a log scale and must not be faked."""
    with pytest.raises(ValueError, match="no lesions at all"):
        lesion_ratio(0.0, 10.0)


def test_lesion_ratio_is_a_plain_ratio():
    assert lesion_ratio(3.0, 12.0) == pytest.approx(0.25)
