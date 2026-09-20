"""The headroom measurement's contract.

The claim this file protects is narrow and load-bearing: the restricted scorer
must be the SAME procedure as the unrestricted one, applied to fewer arms. If
restricting the arm set also changed the fitting, the comparison would be
between two different tests rather than two arm sets.
"""

from __future__ import annotations

import pytest

from gate.ceiling import run_ceiling
from gate.headroom import dial_ceiling, run_headroom


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


# --- the dial ceiling: bounding the exam rather than the model ------------


def test_dial_ceiling_reproduces_the_in_sample_figure():
    """6.6pp vs a 10.6pp null. Computed from trial numbers alone, so this is a
    fixed property of the arm set and should only move when an arm is added."""
    d = dial_ceiling()
    assert round(d["in_sample"]["mae"], 1) == 6.6
    assert round(d["in_sample"]["null_mae"], 1) == 10.6


def test_singleton_groups_are_fitted_exactly_in_sample():
    """Which is why the in-sample row cannot be the honest one."""
    d = dial_ceiling()
    assert set(d["singletons"]) == {"glatiramer acetate", "daclizumab"}
    assert d["in_sample"]["n"] - d["multi_only"]["n"] == len(d["singletons"])


def test_the_headroom_shrinks_when_scored_the_way_everything_else_is():
    """In sample 4.0pp, out of sample 0.8pp. If this inverts, the exam got wider
    and the 'growing the arm set is the cheapest move' conclusion needs redoing."""
    d = dial_ceiling()
    assert d["out_of_sample"]["headroom"] < d["in_sample"]["headroom"]
    assert d["out_of_sample"]["headroom"] < 1.0


def test_a_perfect_dial_model_still_beats_the_null():
    """The representation is not incapable — it is nearly worthless here, which
    is a different claim and the one the plan should carry."""
    d = dial_ceiling()
    assert d["out_of_sample"]["mae"] < d["out_of_sample"]["null_mae"]


def test_the_null_pool_matches_the_scored_set():
    """The restricted rows must not be scored against a null drawn from arms the
    model was never examined on — that hands the null free information. Scoring
    the 10 multi-member arms against a 12-arm null reads 11.9pp; the matched
    convention reads 11.7pp, and the repo carries only the second."""
    d = dial_ceiling()
    assert round(d["out_of_sample"]["null_mae"], 1) == 11.7
    assert round(d["multi_only"]["null_mae"], 1) == 10.5
