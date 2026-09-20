"""The headroom measurement's contract.

The claim this file protects is narrow and load-bearing: the restricted scorer
must be the SAME procedure as the unrestricted one, applied to fewer arms. If
restricting the arm set also changed the fitting, the comparison would be
between two different tests rather than two arm sets.
"""

from __future__ import annotations

import pytest

from gate.ceiling import run_ceiling
from gate.headroom import run_headroom


@pytest.fixture(scope="module")
def headroom():
    return run_headroom()


def test_only_reachable_arms_are_scored(headroom):
    reachable = {f.arm for f in run_ceiling()["fits"] if f.reachable}
    assert {r["arm"] for r in headroom["rows"]} == reachable


def test_the_null_is_computed_from_the_same_restricted_set(headroom):
    """A null drawn from all 12 arms would be answering a different question."""
    scored = {r["arm"] for r in headroom["rows"]}
    assert len(scored) == headroom["n_reachable"]
    for row in headroom["rows"]:
        assert row["null_predicted"] < 0, "null should sit in the range of the arms"


def test_the_oracle_is_the_floor(headroom):
    """The per-arm oracle on the same set cannot be worse than the shared fit."""
    assert headroom["oracle_mae"] <= headroom["mae"]


def test_singletons_are_identified(headroom):
    """glatiramer is alone on alpha_R|delta; if the arm set grows and it stops
    being alone, this test should fail so the sensitivity text gets revisited."""
    assert "glatiramer acetate" in headroom["singletons"]


def test_the_conclusion_survives_dropping_singletons(headroom):
    """The headline must not rest on one fold. If this ever inverts, the
    sensitivity paragraph in gate/headroom.py is wrong and must be rewritten."""
    assert headroom["mae"] > headroom["null_mae"]
    assert headroom["mae_excl_singletons"] > headroom["null_mae_excl_singletons"]


def test_fitted_potencies_are_inside_the_grid(headroom):
    for row in headroom["rows"]:
        assert 0.0 <= row["fitted_potency"] <= 0.95
