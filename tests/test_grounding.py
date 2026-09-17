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
