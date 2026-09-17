"""Guards on the guards (backtest/validation.py).

These are software-correctness tests. They assert that malformed scoring inputs
raise instead of producing a number, which is the failure mode that would put a
meaningless score into an aggregate. They assert nothing about biology.
"""

import numpy as np
import pytest

from backtest.harness import score_delta
from backtest.validation import validate_profile_pair, validate_top_k


def test_profile_guard_rejects_shape_mismatch():
    with pytest.raises(ValueError, match="shape mismatch"):
        validate_profile_pair(np.zeros(2), np.zeros(3))


def test_profile_guard_rejects_non_finite():
    with pytest.raises(ValueError, match="finite"):
        validate_profile_pair(np.array([0.0, np.nan]), np.zeros(2))
    with pytest.raises(ValueError, match="finite"):
        validate_profile_pair(np.zeros(2), np.array([0.0, np.inf]))


def test_profile_guard_rejects_non_numeric():
    with pytest.raises(TypeError, match="numeric"):
        validate_profile_pair(np.array(["a", "b"]), np.zeros(2))


def test_profile_guard_rejects_non_1d():
    with pytest.raises(ValueError, match="one-dimensional"):
        validate_profile_pair(np.zeros((2, 2)), np.zeros((2, 2)))


def test_profile_guard_accepts_numeric_vectors():
    validate_profile_pair(np.zeros(2), np.ones(2))


@pytest.mark.parametrize("value", [0, -1, True, 1.5, "50", None])
def test_top_k_must_be_positive_integer(value):
    with pytest.raises(ValueError):
        validate_top_k(value)


def test_top_k_accepts_positive_integer():
    validate_top_k(1)
    validate_top_k(np.int64(50))


# --- the guards are wired into the scorer, not decoration --------------------

def test_score_delta_rejects_nan_instead_of_scoring_it():
    """A NaN delta used to travel through and land in the aggregate."""
    with pytest.raises(ValueError, match="finite"):
        score_delta(np.array([1.0, np.nan, 3.0]), np.array([1.0, 2.0, 3.0]))


def test_score_delta_rejects_meaningless_top_k():
    with pytest.raises(ValueError, match="top_k"):
        score_delta(np.zeros(5), np.zeros(5), top_k=0)


def test_score_delta_still_scores_valid_input():
    out = score_delta(np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0]), top_k=2)
    assert out["delta_pearson"] == pytest.approx(1.0)
    assert out["delta_mse"] == pytest.approx(0.0)
