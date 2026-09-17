"""Invariants on the response curve and the leave-one-arm-out test.

These run off the cached curve (`results/response_curve.json`), so they cost
milliseconds rather than the twelve minutes the grid itself takes. What they
protect is that the curve is the model tabulated -- not a fitted stand-in for it
-- and that the LOO keeps reporting a null beside its own score.
"""

from __future__ import annotations

import math

import pytest

from backtest.loo import _predicted_change, _quantified_arms, run_loo
from backtest.response_curve import CACHE, load, predict

pytestmark = pytest.mark.skipif(
    not CACHE.exists(),
    reason="response curve cache absent; build with `PYTHONPATH=. python -m backtest.response_curve`")


def test_curve_is_exact_at_its_grid_points():
    """Interpolation must not disturb a point the simulation actually produced."""
    table = load()
    points = table["curve"]["0.00"]
    for key in ("0.00", "0.50", "1.00"):
        assert predict(float(key), 0.0, table) == pytest.approx(points[key])


def test_untreated_is_the_zero_of_the_curve():
    table = load()
    assert predict(0.0, 0.0, table) == pytest.approx(0.0, abs=1e-9)


def test_placebo_controlled_arm_reads_straight_off_the_curve():
    table = load()
    ifnb = next(o for o in _quantified_arms() if o.arm == "IFN-beta")
    assert _predicted_change(ifnb, 0.5, table) == pytest.approx(predict(0.5, 0.0, table))


def test_active_comparator_arm_is_exactly_zero_under_a_class_strength():
    """The structural failure, pinned as a test rather than left as a remark.

    Both arms of OPERA get the same `treat` while one strength serves the whole
    class, so the predicted difference is 0 at every strength. When per-drug
    potency lands (blocker 4), this test should start failing -- that is the
    signal that the rule has gained drug-specific information, and the assertion
    is where to record it.
    """
    table = load()
    ocre = next(o for o in _quantified_arms() if o.arm == "ocrelizumab")
    for k in (0.2, 0.5, 0.78):
        assert _predicted_change(ocre, k, table) == pytest.approx(0.0, abs=1e-9)


def test_full_suppression_leaves_an_active_comparison_undefined():
    """At treat=1.0 the comparator arm has no relapses left, so the ratio is undefined.

    Reported as NaN and excluded from fits rather than clamped to a number that
    would look like a prediction.
    """
    table = load()
    ocre = next(o for o in _quantified_arms() if o.arm == "ocrelizumab")
    assert math.isnan(_predicted_change(ocre, 1.0, table))


def test_loo_always_reports_a_null_beside_its_score():
    """docs/QUALITY.md: a benchmark without a null is not a benchmark."""
    r = run_loo()
    assert r["null_mae"] > 0
    assert len(r["rows"]) == len(_quantified_arms())
    for row in r["rows"]:
        assert "null_error" in row
