"""Tests for the candidate kill filter.

The load-bearing test is `test_ranking_is_refused`. Everything else in this repo
is a measurement; that one is a guardrail, and it exists because a ranked list of
candidates is the artifact somebody would screenshot, and the model has not
earned one.
"""

from __future__ import annotations

import pytest

from bricks.profiles import PROFILES, touched_points
from bricks.qsp_velez import INTERVENTION_POINTS, MechanismProfile
from screen import kill_filter
from screen.kill_filter import (
    NOISE_FLOOR,
    KillReason,
    enumerate_candidates,
    rank_candidates,
    screen,
)

# --------------------------------------------------------------------------- #
# the guardrail
# --------------------------------------------------------------------------- #

def test_ranking_is_refused():
    """Gated on backtest/lomo.py beating its predict-the-mean null. It does not.

    The refusal lives in code rather than only in a docstring, because a
    docstring does not stop anyone.
    """
    with pytest.raises(NotImplementedError, match="not supported"):
        rank_candidates([])


def test_the_refusal_names_the_number_that_would_lift_it():
    try:
        rank_candidates()
    except NotImplementedError as exc:
        msg = str(exc)
    assert "lomo" in msg.lower()
    assert "null" in msg.lower()
    assert "45.9" in msg


# --------------------------------------------------------------------------- #
# enumeration
# --------------------------------------------------------------------------- #

def test_enumeration_is_exhaustive_not_heuristic():
    """A search heuristic over a space you can enumerate hides which corners
    were never visited. Seven dials is enumerable."""
    singles = enumerate_candidates(max_points=1)
    assert len(singles) == 2 * len(INTERVENTION_POINTS)
    pairs = enumerate_candidates(max_points=2)
    n = len(INTERVENTION_POINTS)
    assert len(pairs) == 2 * n + 4 * (n * (n - 1) // 2)


def test_every_candidate_moves_only_its_own_dials():
    for cand in enumerate_candidates(max_points=2, potency=0.5):
        moved = touched_points(cand)
        assert 1 <= len(moved) <= 2
        for p in INTERVENTION_POINTS:
            if p not in moved:
                assert getattr(cand, p) == 1.0


def test_both_directions_are_enumerated():
    labels = {c.label for c in enumerate_candidates(max_points=1)}
    assert "alpha_R+" in labels and "alpha_R-" in labels


# --------------------------------------------------------------------------- #
# the filters
# --------------------------------------------------------------------------- #

def _fake_damage(monkeypatch, mapping, default=1.0):
    def fake(profile, carrying_capacity=None):
        return mapping.get(profile.label, default)
    monkeypatch.setattr(kill_filter, "_median_damage", fake)


def test_a_diverging_candidate_is_killed_as_out_of_regime(monkeypatch):
    def fake(profile, carrying_capacity=None):
        return None if profile.label == "boom" else 1.0
    monkeypatch.setattr(kill_filter, "_median_damage", fake)
    res = screen([MechanismProfile(label="boom", gamma_R=4.0, source="t")])
    assert res[0].killed_by is KillReason.OUT_OF_REGIME
    assert not res[0].survived


def test_a_candidate_that_never_helps_is_killed_as_unreachable(monkeypatch):
    # untreated 1.0, candidate always worse
    _fake_damage(monkeypatch, {"untreated": 1.0}, default=2.0)
    res = screen([MechanismProfile(label="useless", alpha_E=1.5, source="t")])
    assert res[0].killed_by is KillReason.UNREACHABLE


def test_an_improvement_smaller_than_the_noise_floor_does_not_count(monkeypatch):
    """A difference under the measured noise floor is not a result."""
    _fake_damage(monkeypatch, {"untreated": 1.0}, default=1.0 - NOISE_FLOOR / 2)
    res = screen([MechanismProfile(label="marginal", alpha_E=0.5, source="t")])
    assert res[0].killed_by is KillReason.UNREACHABLE


def test_a_candidate_matching_an_existing_drug_is_killed_as_degenerate(monkeypatch):
    _fake_damage(monkeypatch, {"untreated": 1.0}, default=0.1)
    # alpha_E alone is IFN-beta / teriflunomide / dimethyl fumarate
    res = screen([MechanismProfile(label="me-too", alpha_E=0.5, source="t")])
    assert res[0].killed_by is KillReason.DEGENERATE
    assert res[0].like_existing


def test_a_genuinely_new_and_reachable_candidate_survives(monkeypatch):
    _fake_damage(monkeypatch, {"untreated": 1.0}, default=0.1)
    novel = MechanismProfile(label="novel", gamma_R=0.5, ke=0.5, source="t")
    assert touched_points(novel) not in {touched_points(p) for p in PROFILES.values()}
    res = screen([novel])
    assert res[0].survived
    assert res[0].killed_by is None


# --------------------------------------------------------------------------- #
# reporting honesty
# --------------------------------------------------------------------------- #

def test_the_reported_damage_is_the_one_that_passed_the_filter(monkeypatch):
    """A presentation bug caught in review: the screen probes several potencies
    and passes on the BEST, but first reported the candidate's default-potency
    damage — so some survivors printed worse-than-untreated numbers."""
    def fake(profile, carrying_capacity=None):
        if profile.label == "untreated":
            return 1.0
        return 0.05 if profile.source == "probe" else 3.0
    monkeypatch.setattr(kill_filter, "_median_damage", fake)

    res = screen([MechanismProfile(label="late-bloomer", gamma_R=0.5, ke=0.5,
                                   source="t")])[0]
    assert res.survived
    assert res.best_damage == 0.05
    assert res.median_damage == 3.0
    assert res.best_damage < res.untreated_damage


def test_nothing_claims_to_be_validated(monkeypatch):
    _fake_damage(monkeypatch, {"untreated": 1.0}, default=0.1)
    res = screen([MechanismProfile(label="novel", gamma_R=0.5, ke=0.5, source="t")])
    assert res[0].validated is False


def test_a_kill_always_carries_its_reason(monkeypatch):
    _fake_damage(monkeypatch, {"untreated": 1.0}, default=2.0)
    res = screen([MechanismProfile(label="useless", alpha_E=1.5, source="t")])
    assert res[0].detail.strip()


def test_screen_refuses_when_the_untreated_arm_diverges(monkeypatch):
    monkeypatch.setattr(kill_filter, "_median_damage",
                        lambda p, carrying_capacity=None: None)
    with pytest.raises(RuntimeError, match="nothing to screen against"):
        screen([MechanismProfile(label="x", alpha_E=0.5, source="t")])


def test_results_are_returned_in_input_order(monkeypatch):
    _fake_damage(monkeypatch, {"untreated": 1.0}, default=0.1)
    cands = [MechanismProfile(label=f"c{i}", gamma_R=0.5, ke=0.5, source="t")
             for i in range(3)]
    assert [r.profile.label for r in screen(cands)] == ["c0", "c1", "c2"]
