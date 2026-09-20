"""Tests for the regulatory-cell harm channel.

The point of this channel is that it reads NO trial outcome, so the first test
is a tripwire on exactly that. The rest guard the expression floor — which
throws out the most impressive-looking number in the analysis — and pin the
statistical claim at its real, unimpressive size so it cannot drift upward in
the retelling.
"""

from __future__ import annotations

import pytest

from bricks import harm_channel
from bricks.harm_channel import (
    ARM_TARGET,
    EFFECTOR_SUBSETS,
    HPA_NTPM,
    MIN_EFFECTOR_NTPM,
    ranking,
    score,
    top_k_probability,
)

# --------------------------------------------------------------------------- #
# the property that makes the channel worth having
# --------------------------------------------------------------------------- #

def test_scoring_reads_no_trial_outcome(monkeypatch):
    """If KNOWN_OUTCOMES is touched while scoring, this is not an independent
    channel — it is the gate's own answer wearing a different hat."""
    import backtest.clinical as clinical

    class Tripwire(list):
        def __iter__(self):
            raise AssertionError("the harm channel read KNOWN_OUTCOMES")

    monkeypatch.setattr(clinical, "KNOWN_OUTCOMES", Tripwire())
    assert score("daclizumab").treg_ratio > 1.0
    assert len(ranking()) >= 5


def test_the_channel_sees_the_two_cases_nothing_else_here_can():
    """Lenercept is inert through MRI; daclizumab is a success through ARR.

    Both are ranked as liabilities here, which is the entire reason this module
    exists.
    """
    top_two = [s.arm for s in ranking()[:2]]
    assert set(top_two) == {"daclizumab", "lenercept"}


def test_daclizumab_scores_highest():
    """Its target IS the canonical Treg marker, and it was withdrawn for fatal
    encephalitis while REDUCING relapses 45%."""
    assert ranking()[0].arm == "daclizumab"
    assert ranking()[0].target == "IL2RA"


# --------------------------------------------------------------------------- #
# the expression floor
# --------------------------------------------------------------------------- #

def test_atacicept_is_excluded_rather_than_ranked_first():
    """The best-looking result in the analysis, thrown out on purpose.

    TACI's raw ratio is 88x only because its effector-T denominator is ~0.025
    nTPM. It is a B-cell gene (52-131 nTPM in B cells). A ratio of two
    near-zero numbers is division by noise, not a signal about regulation.
    """
    s = score("atacicept")
    assert s.classifiable is False
    assert s.treg_ratio is None
    assert "below the" in s.reason
    assert "atacicept" not in [r.arm for r in ranking()]


def test_the_floor_is_the_conventional_not_expressed_threshold():
    assert MIN_EFFECTOR_NTPM == 1.0


def test_an_unclassifiable_arm_never_leaks_into_the_ranking():
    for s in ranking():
        assert s.classifiable is True
        assert s.treg_ratio is not None


# --------------------------------------------------------------------------- #
# the statistical claim, pinned at its real size
# --------------------------------------------------------------------------- #

def test_the_probability_is_the_unimpressive_one():
    """1/C(6,2) = 0.0667. It does NOT reach 0.05 and must not be restated as if
    it did — an earlier draft of the docstring claimed 1/21 = 0.048 by counting
    seven targets instead of six."""
    rows = ranking()
    assert len(rows) == 6
    p = top_k_probability(len(rows), 2)
    assert p == pytest.approx(1 / 15)
    assert p > 0.05


def test_top_k_probability_is_exact_combinatorics():
    assert top_k_probability(6, 2) == pytest.approx(1 / 15)
    assert top_k_probability(8, 3) == pytest.approx(1 / 56)
    assert top_k_probability(5, 5) == pytest.approx(1.0)


def test_top_k_probability_rejects_nonsense():
    for n, k in ((5, 0), (5, 6), (0, 1)):
        with pytest.raises(ValueError):
            top_k_probability(n, k)


# --------------------------------------------------------------------------- #
# data integrity
# --------------------------------------------------------------------------- #

def test_every_mapped_arm_has_expression_data():
    for arm, target in ARM_TARGET.items():
        assert target in HPA_NTPM, f"{arm} maps to {target} with no data"


def test_every_target_carries_the_effector_subsets():
    for target, v in HPA_NTPM.items():
        missing = [c for c in EFFECTOR_SUBSETS if c not in v]
        assert not missing, f"{target} missing {missing}"
        assert "T-reg" in v


def test_separation_is_not_just_drug_efficacy():
    """CD52 sits at 1.81 and alemtuzumab is among the most effective DMTs there
    is. The claim is a regulatory liability, not a potency ranking."""
    by_arm = {s.arm: s for s in ranking()}
    assert by_arm["alemtuzumab"].treg_ratio > 1.0
    assert by_arm["alemtuzumab"].treg_ratio < by_arm["lenercept"].treg_ratio


def test_nothing_is_marked_validated():
    assert score("daclizumab").validated is False


def test_unmapped_arm_is_refused():
    with pytest.raises(KeyError):
        score("aspirin")


def test_one_row_per_target_not_per_arm():
    """Fingolimod and ponesimod share S1PR1; ocrelizumab and ofatumumab share
    MS4A1. Counting them twice would inflate the denominator of the test above."""
    targets = [s.target for s in ranking()]
    assert len(targets) == len(set(targets))


def test_the_transcribed_values_are_reproducible():
    """scripts/derive_harm_channel.py --check re-downloads and compares.

    Not run here (it is a network call); this asserts the script exists, because
    a constant whose provenance script has been deleted is an orphan number.
    """
    from pathlib import Path
    script = Path(harm_channel.__file__).resolve().parent.parent / "scripts" / \
        "derive_harm_channel.py"
    assert script.is_file()
