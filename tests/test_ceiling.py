"""The oracle ceiling's contract — mostly that it cannot flatter itself.

The failure mode this file guards is a ceiling that reads as a calibration
result while being a sign error, so the tests are about the SPLIT: an arm the
model cannot reach must never be counted as a fitted one, and the headline must
be the reachable subset.
"""

from __future__ import annotations

import dataclasses

import pytest

from gate.ceiling import GRID_HI, GRID_LO, OracleFit, recoverability, run_ceiling


@pytest.fixture(scope="module")
def ceiling():
    return run_ceiling()


def test_every_quantified_arm_is_classified(ceiling):
    """12 quantified arms, each either reachable or not. No arm silently dropped."""
    assert len(ceiling["fits"]) == 12
    assert ceiling["n_reachable"] + ceiling["n_unreachable"] == 12


def test_grid_edge_is_never_counted_as_a_fit(ceiling):
    for f in ceiling["fits"]:
        at_edge = f.potency <= GRID_LO + 1e-6 or f.potency >= GRID_HI - 1e-6
        assert f.reachable is not at_edge, f"{f.arm} classified against its own potency"


def test_reachable_arms_are_fitted_essentially_exactly(ceiling):
    """One free parameter, one target. If this drifts, the curve stopped passing
    through the trial numbers and the ceiling means something different."""
    for f in ceiling["fits"]:
        if f.reachable:
            assert f.error < 3.0, f"{f.arm} oracle error {f.error:.1f}pp"


def test_unreachable_arms_predict_no_effect(ceiling):
    """Their best fit is the boundary, which on these dials predicts ~0%."""
    for f in ceiling["fits"]:
        if not f.reachable:
            assert abs(f.predicted) < 5.0


def test_the_combined_figure_is_not_the_headline(ceiling):
    """Mixing the two subsets must move the number a lot -- if it did not, the
    split would be cosmetic and this module would not need to make it."""
    assert ceiling["mae_all"] > ceiling["mae_reachable"] * 5


def test_recoverability_scores_only_arms_with_both_estimates():
    rec = recoverability()
    ceil = run_ceiling()
    reachable = {f.arm for f in ceil["fits"] if f.reachable}
    for row in rec["rows"]:
        assert row["arm"] in reachable, "scored an arm with no oracle potency"


def test_recoverability_is_an_out_of_fit_prediction():
    """The MRI potency never saw the relapse number, so its error must be a real
    error -- not the ~0 the oracle gets by construction."""
    rec = recoverability()
    assert rec["rows"], "no arm scorable; the comparison has gone empty"
    assert rec["mae"] > 1.0


def test_oracle_fit_is_frozen():
    f = OracleFit(arm="x", known=-30.0, potency=0.5, predicted=-30.0,
                  error=0.0, null_error=10.0, reachable=True)
    with pytest.raises(dataclasses.FrozenInstanceError):
        f.potency = 0.6  # type: ignore[misc]
