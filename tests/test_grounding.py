"""Citation guards for the intervention class strengths (bricks/grounding.py).

These do not test that the strengths are biologically right. They test the
claims the grounding rests on — the ones that would rot silently:

  1. SUPPRESSIVE_STRENGTH is the number its own cited arithmetic produces. A
     hand-nudged constant with the Kang citation still sitting above it is the
     easiest way for "grounded" to quietly become "fitted";
  2. the recorded derivation inputs still match the data, when the Kang cache is
     present (the slow, real check);
  3. IMMUNOGENIC_STRENGTH is still the ungrounded, reasoned value — it has no
     independent source, and the moment it acquires one this test should be the
     thing that fails;
  4. the strengths stay CLASS-level: two suppressive drugs cannot drift apart.

Run:  python -m pytest tests/test_grounding.py -q
"""

from __future__ import annotations

import os

import pytest

from bricks.grounding import (
    IMMUNOGENIC,
    IMMUNOGENIC_STRENGTH,
    KANG_IDENTITY_DISTANCE,
    KANG_IFNB_DELTA_NORM,
    NEUTRAL,
    SUPPRESSIVE,
    SUPPRESSIVE_STRENGTH,
    mechanism_to_params,
)

_KANG_CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           "data", "cache", "kang_2018.h5ad")


def test_suppressive_strength_is_its_own_arithmetic():
    """0.78 must BE 15.902 / 20.325, not a number sitting next to that citation."""
    assert KANG_IFNB_DELTA_NORM == pytest.approx(15.902)
    assert KANG_IDENTITY_DISTANCE == pytest.approx(20.325)
    assert SUPPRESSIVE_STRENGTH == 0.78
    assert SUPPRESSIVE_STRENGTH == round(KANG_IFNB_DELTA_NORM / KANG_IDENTITY_DISTANCE, 2)


@pytest.mark.skipif(not os.path.exists(_KANG_CACHE),
                    reason="Kang cache absent; run backtest.run_kang once to fetch it")
def test_derivation_still_reproduces_from_the_kang_data():
    """The slow guard: recompute from Kang and demand the module still matches.

    If the loader, the reliability bar or the cached dataset changes, the
    constant above stops being derived and starts being folklore. This is the
    test that notices.
    """
    from scripts.derive_suppressive_strength import derive

    d = derive()
    assert d["dropped_unreliable"] == ["Megakaryocytes"]
    assert len(d["cell_types"]) == 7
    assert d["mean_delta_norm"] == pytest.approx(KANG_IFNB_DELTA_NORM, abs=0.01)
    assert d["mean_identity_distance"] == pytest.approx(KANG_IDENTITY_DISTANCE, abs=0.01)
    assert d["strength"] == SUPPRESSIVE_STRENGTH


def test_immunogenic_strength_is_still_the_ungrounded_one():
    """No independent source exists for it (Kang has no encephalitogenic arm).

    It is deliberately unchanged and deliberately labelled. If someone grounds
    it, this assertion is where they will land — update it WITH the derivation,
    never without one.
    """
    assert IMMUNOGENIC_STRENGTH == 0.4
    assert IMMUNOGENIC_STRENGTH != round(KANG_IFNB_DELTA_NORM / KANG_IDENTITY_DISTANCE, 2)


def test_strengths_stay_class_level_not_per_drug():
    """The honesty property of the rule: a class has ONE number, whatever the arm."""
    assert mechanism_to_params(SUPPRESSIVE) == (SUPPRESSIVE_STRENGTH, 0.0)
    assert mechanism_to_params(IMMUNOGENIC) == (0.0, IMMUNOGENIC_STRENGTH)
    assert mechanism_to_params(NEUTRAL) == (0.0, 0.0)

    from bricks.intervention import APL_CGP77116, GLATIRAMER, IFN_BETA
    assert IFN_BETA.treat == GLATIRAMER.treat == SUPPRESSIVE_STRENGTH
    assert APL_CGP77116.immunogenic == IMMUNOGENIC_STRENGTH


def test_regulation_disrupting_drug_gets_both_channels_and_no_new_constant():
    """A drug whose target also carries a regulatory function suppresses AND provokes.

    The harm channel must be exactly IMMUNOGENIC_STRENGTH. A third constant, tuned
    until the known-harmful arms come out harmful, would be fitted to the outcome
    the clinical gate exists to test — this assertion is the tripwire for that.
    """
    assert mechanism_to_params(SUPPRESSIVE, disrupts_regulation=True) == (
        SUPPRESSIVE_STRENGTH, IMMUNOGENIC_STRENGTH)

    from bricks.intervention import LENERCEPT
    assert LENERCEPT.treat == SUPPRESSIVE_STRENGTH
    assert LENERCEPT.immunogenic == IMMUNOGENIC_STRENGTH


def test_atacicept_is_left_unflagged_as_the_honest_gap():
    """No independent source pins BAFF/APRIL blockade to a lost regulatory function.

    The anti-CD20 Breg evidence (PMID 18802481) is a different target with a timing
    dependence this model cannot express. Flagging atacicept anyway would be
    reasoning backwards from ATAMS, so it stays unflagged and the gate keeps failing
    it. Ground it properly and this test is where to update — with the citation.
    """
    from bricks.intervention import ATACICEPT
    assert ATACICEPT.immunogenic == 0.0
