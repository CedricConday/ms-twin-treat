"""Invariants on the clinical gate's arm set (BUILD_PLAN §8 blocker 3).

These are cheap structural checks: they never run the pipeline. What they protect
is the separation the gate rests on — every arm exists as a simulable
intervention, every outcome carries a citation, and the two counterexamples that
make the gate falsifiable are still in it.
"""

from __future__ import annotations

from backtest.clinical import KNOWN_OUTCOMES
from bricks.grounding import SUPPRESSIVE
from bricks.intervention import LIBRARY


def test_every_anchor_arm_and_comparator_is_simulable():
    """An outcome naming an arm the library cannot build is scored against nothing."""
    for o in KNOWN_OUTCOMES:
        assert o.arm in LIBRARY, f"{o.arm} has a trial anchor but no intervention"
        assert o.comparator in LIBRARY, f"{o.arm} is scored against unknown arm {o.comparator}"


def test_every_outcome_cites_its_trial():
    """The numbers are the exam; an uncited one cannot be checked by a reader."""
    for o in KNOWN_OUTCOMES:
        if o.arm == "untreated":
            continue
        assert "PMID" in o.source, f"{o.arm} outcome has no PMID"


def test_active_comparator_arms_are_not_scored_against_untreated():
    """OPERA, CARE-MS and OPTIMUM ran against an active drug. Scoring them against
    the untreated arm would silently compare a number to a different experiment."""
    expected = {"ocrelizumab": "IFN-beta", "alemtuzumab": "IFN-beta",
                "ponesimod": "teriflunomide"}
    for o in KNOWN_OUTCOMES:
        if o.arm in expected:
            assert o.comparator == expected[o.arm]


def test_suppressive_but_harmful_counterexamples_are_present():
    """lenercept and atacicept are immunosuppressive by mechanism and harmed patients.

    They are the arms the current mechanism-class rule cannot get right, and the
    reason the gate can fail at all. Removing either one to make the gate greener
    would remove the test, not fix the model.
    """
    by_arm = {o.arm: o for o in KNOWN_OUTCOMES}
    for arm in ("lenercept", "atacicept"):
        assert arm in by_arm, f"{arm} counterexample removed from the gate"
        assert by_arm[arm].direction == "harms"
        assert LIBRARY[arm].mechanism == SUPPRESSIVE
