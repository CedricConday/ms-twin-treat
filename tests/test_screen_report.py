"""Tests for the screen's report half.

The load-bearing one is `test_survivors_are_not_ordered_by_benefit`. Everything
else here is shape; that one is the guardrail, because an ordered list of
candidates with a damage column beside it IS a ranking whatever the header says,
and this model loses to a predict-the-mean null on unseen mechanisms.

None of these run the real screen — that is ~3,000 simulations. `screen` and
`enumerate_candidates` are monkeypatched so the ordering logic is tested
directly.
"""

from __future__ import annotations

import pytest

from bricks.qsp_velez import MechanismProfile
from screen import report
from screen.kill_filter import KillReason, ScreenResult


def _profile(label: str) -> MechanismProfile:
    return MechanismProfile(label=label, source="test", alpha_E=0.5)


@pytest.fixture
def fake(monkeypatch):
    """Three survivors whose alphabetical order DISAGREES with their damage."""
    rows = [
        ScreenResult(_profile("zeta+"), None, "survived all four filters",
                     median_damage=0.10, best_damage=0.10, best_potency=0.2,
                     untreated_damage=1.0),
        ScreenResult(_profile("alpha+"), None, "survived all four filters",
                     median_damage=0.90, best_damage=0.90, best_potency=0.5,
                     untreated_damage=1.0),
        ScreenResult(_profile("mid+"), None, "survived all four filters",
                     median_damage=0.50, best_damage=0.50, best_potency=0.8,
                     untreated_damage=1.0),
        ScreenResult(_profile("dead+"), KillReason.OUT_OF_REGIME, "diverged",
                     untreated_damage=1.0),
    ]
    monkeypatch.setattr(report, "enumerate_candidates",
                        lambda **kw: [r.profile for r in rows])
    monkeypatch.setattr(report, "screen", lambda c, **kw: rows)
    return rows


def test_survivors_are_not_ordered_by_benefit(fake):
    """Alphabetical, even though zeta+ is by far the 'best' simulated candidate."""
    rep = report.run()
    labels = [s["label"] for s in rep["survivors"]]
    assert labels == ["alpha+", "mid+", "zeta+"]
    best_first = [s["label"] for s in
                  sorted(rep["survivors"], key=lambda s: s["best_damage"])]
    assert labels != best_first, "ordering coincides with a ranking; pick a harder case"
    assert "NOT by predicted benefit" in rep["ordering"]


def test_the_caveats_name_the_number_that_gates_ranking(fake):
    joined = " ".join(rep_c for rep_c in report.run()["caveats"])
    assert "lomo" in joined.lower()
    assert "45.9" in joined and "12.3" in joined


def test_the_blind_spot_is_stated_not_implied(fake):
    joined = " ".join(report.run()["caveats"]).lower()
    for drug in ("natalizumab", "fingolimod", "alemtuzumab"):
        assert drug in joined


def test_killed_candidates_keep_their_reason(fake):
    rep = report.run()
    assert rep["counts"]["OUT_OF_REGIME"] == 1
    assert rep["killed"][0]["detail"] == "diverged"


def test_the_report_never_claims_validation(fake):
    rep = report.run()
    assert rep["validated"] is False
    assert "validated=False" in report.to_markdown(rep)


def test_the_report_states_how_many_survivors_beat_the_best_real_drug(fake):
    """The survivor table's most misleading number, computed rather than implied.

    32 of 38 real survivors claim a bigger effect than natalizumab's -68% in
    AFFIRM, and 22 claim better than -90%, from a model that cannot reproduce
    interferon beta's -30%. If that comparison is not in the artifact, the
    artifact invites the misreading.
    """
    rep = report.run()
    imp = rep["implausibility"]
    assert imp["n"] == 3
    assert imp["beating_best_real_arm"] == 1   # only zeta+ at -90% clears -68%
    md = report.to_markdown(rep)
    assert "not credible" in md.lower()
    assert "natalizumab" in md.lower()
