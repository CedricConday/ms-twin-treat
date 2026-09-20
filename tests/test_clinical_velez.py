"""Tests for the direction gate scored on the grounded stack.

The expensive part is the simulation, so these tests exercise the scoring logic
against injected damages rather than running 128 histories per arm. What they
protect is the part that would quietly flatter the result: that both gates use
the SAME rules and the SAME cited outcomes, that an unscoreable arm counts as a
miss rather than being dropped, and that a row passing for the wrong reason
stays labelled.
"""

from __future__ import annotations

import numpy as np
import pytest

from backtest import clinical, clinical_velez


def _inject(monkeypatch, damages_by_arm: dict[str, float]):
    """Give every arm a flat per-seed damage, so ratios are exact."""
    seeds = list(clinical_velez.SEEDS)

    def fake(arm: str) -> dict[int, float]:
        if arm not in damages_by_arm:
            return {}
        return {s: damages_by_arm[arm] for s in seeds}

    monkeypatch.setattr(clinical_velez, "_arm_damages", fake)


# --------------------------------------------------------------------------- #
# same rules, same outcomes — or it is not a comparison
# --------------------------------------------------------------------------- #

def test_both_gates_score_the_same_arms_and_outcomes():
    assert clinical_velez.KNOWN_OUTCOMES is clinical.KNOWN_OUTCOMES


def test_both_gates_use_the_same_thresholds():
    assert clinical_velez.NEUTRAL_BAND == clinical.NEUTRAL_BAND
    assert clinical_velez.MAG_TOLERANCE == clinical.MAG_TOLERANCE
    assert clinical_velez._direction is clinical._direction


def test_active_comparator_arms_are_scored_against_their_own_comparator():
    """OPERA ran against IFN beta-1a, not placebo. Scoring it against untreated
    would invent a comparison the trial never made."""
    by_arm = {o.arm: o for o in clinical.KNOWN_OUTCOMES}
    assert by_arm["ocrelizumab"].comparator == "IFN-beta"
    assert by_arm["alemtuzumab"].comparator == "IFN-beta"
    assert by_arm["ponesimod"].comparator == "teriflunomide"


# --------------------------------------------------------------------------- #
# scoring behaviour
# --------------------------------------------------------------------------- #

def test_a_drug_that_lowers_damage_is_scored_as_an_improvement(monkeypatch):
    _inject(monkeypatch, {"untreated": 10.0, "IFN-beta": 2.0})
    row = next(r for r in clinical_velez.run_gate()["rows"] if r["arm"] == "IFN-beta")
    assert row["rr_lesion"] == pytest.approx(0.2)
    assert row["predicted"] < 0
    assert row["predicted_direction"] == "improves"


def test_a_drug_that_raises_damage_is_scored_as_harm(monkeypatch):
    _inject(monkeypatch, {"untreated": 2.0, "IFN-beta": 10.0})
    row = next(r for r in clinical_velez.run_gate()["rows"] if r["arm"] == "IFN-beta")
    assert row["predicted"] > 0
    assert row["predicted_direction"] == "harms"


def test_the_ratio_is_paired_per_seed_not_a_ratio_of_means(monkeypatch):
    """Paired first, then median — the reverse would be a different statistic."""
    seeds = list(clinical_velez.SEEDS)

    def fake(arm):
        if arm == "untreated":
            return {s: 10.0 + s for s in seeds}
        if arm == "IFN-beta":
            return {s: (10.0 + s) * 0.25 for s in seeds}
        return {}

    monkeypatch.setattr(clinical_velez, "_arm_damages", fake)
    row = next(r for r in clinical_velez.run_gate()["rows"] if r["arm"] == "IFN-beta")
    assert row["rr_lesion"] == pytest.approx(0.25)


def test_an_arm_with_no_usable_history_is_undefined_and_counted_as_a_miss(monkeypatch):
    """Lenercept's real behaviour: every history leaves the model's regime.

    An arm nothing can be said about must not quietly leave the denominator —
    that would turn a failure to simulate into a better-looking score.
    """
    _inject(monkeypatch, {"untreated": 10.0})   # every drug arm returns {}
    result = clinical_velez.run_gate()
    assert all(r["predicted_direction"] == "undefined" for r in result["rows"])
    assert result["direction_hits"] == 0
    assert result["n_arms"] == len(clinical.KNOWN_OUTCOMES) - 1


def test_undefined_arms_stay_in_the_denominator(monkeypatch):
    """The headline is X/n_arms, not X/n_scored."""
    _inject(monkeypatch, {"untreated": 10.0, "IFN-beta": 2.0})
    result = clinical_velez.run_gate()
    assert result["n_arms"] > result["n_scored"]
    assert result["n_arms"] == len(clinical.KNOWN_OUTCOMES) - 1


def test_the_blind_spot_flag_rides_along(monkeypatch):
    """A near-null lesion effect is where the map cannot tell 'nothing' from harm."""
    _inject(monkeypatch, {"untreated": 10.0, "IFN-beta": 10.0})
    row = next(r for r in clinical_velez.run_gate()["rows"] if r["arm"] == "IFN-beta")
    assert row["blind_spot"] is True


def test_magnitude_is_only_scored_where_a_trial_number_exists(monkeypatch):
    _inject(monkeypatch, {"untreated": 10.0, "IFN-beta": 2.0})
    result = clinical_velez.run_gate()
    quantified = [o for o in clinical.KNOWN_OUTCOMES
                  if o.relapse_change_pct is not None and o.arm != "untreated"]
    assert result["n_quantified"] <= len(quantified)


def test_dropped_histories_are_reported_not_hidden(monkeypatch):
    seeds = list(clinical_velez.SEEDS)

    def fake(arm):
        if arm == "untreated":
            return {s: 10.0 for s in seeds}
        if arm == "IFN-beta":
            return {s: 2.0 for s in seeds[:10]}
        return {}

    monkeypatch.setattr(clinical_velez, "_arm_damages", fake)
    row = next(r for r in clinical_velez.run_gate()["rows"] if r["arm"] == "IFN-beta")
    assert row["n_seeds"] == 10
    assert row["n_dropped"] == len(seeds) - 10


def test_predicted_change_is_finite_for_a_scored_arm(monkeypatch):
    _inject(monkeypatch, {"untreated": 10.0, "IFN-beta": 2.0})
    row = next(r for r in clinical_velez.run_gate()["rows"] if r["arm"] == "IFN-beta")
    assert np.isfinite(row["predicted"])
