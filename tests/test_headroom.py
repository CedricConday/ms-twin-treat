"""The headroom measurement's contract.

The claim this file protects is narrow and load-bearing: the restricted scorer
must be the SAME procedure as the unrestricted one, applied to fewer arms. If
restricting the arm set also changed the fitting, the comparison would be
between two different tests rather than two arm sets.
"""

from __future__ import annotations

import statistics

import pytest

from backtest.clinical import KNOWN_OUTCOMES
from backtest.lomo import mechanism_groups
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


def _expected_in_sample():
    """Recompute the in-sample ceiling here, independently of gate.headroom.

    Deliberately a second implementation rather than a pinned literal. Pinning
    6.6pp caught real drift once, but it churns every time an arm is wired --
    and a test that cries wolf gets deleted, taking the invariant with it. A
    parallel implementation keeps the check and survives the arm set growing.
    """
    known = {o.arm: o.relapse_change_pct for o in KNOWN_OUTCOMES
             if o.relapse_change_pct is not None and o.arm != "untreated"}
    grand = statistics.fmean(known.values())
    errs, nulls = [], []
    for arms in mechanism_groups().values():
        vals = [known[a] for a in arms]
        mu = statistics.fmean(vals)
        errs += [abs(v - mu) for v in vals]
        nulls += [abs(v - grand) for v in vals]
    return statistics.fmean(errs), statistics.fmean(nulls), len(errs)


def test_dial_ceiling_matches_an_independent_computation():
    """The ceiling is arithmetic on the trial numbers; two implementations of it
    must agree to the decimal at whatever arm count the table currently holds."""
    mae, null_mae, n = _expected_in_sample()
    d = dial_ceiling()
    assert round(d["in_sample"]["mae"], 6) == round(mae, 6)
    assert round(d["in_sample"]["null_mae"], 6) == round(null_mae, 6)
    assert d["in_sample"]["n"] == n


def test_singleton_groups_are_fitted_exactly_in_sample():
    """Which is why the in-sample row cannot be the honest one. Membership is
    derived rather than listed, so wiring an arm onto a singleton dial -- the
    cheapest way to widen the exam -- does not fail this test spuriously."""
    d = dial_ceiling()
    expected = {arms[0] for arms in mechanism_groups().values() if len(arms) == 1}
    assert set(d["singletons"]) == expected
    assert d["in_sample"]["n"] - d["multi_only"]["n"] == len(d["singletons"])


def test_the_headroom_shrinks_when_scored_the_way_everything_else_is():
    """The invariant, not the number: scoring a perfect dial model the way every
    other scorer here is scored must leave less headroom than the in-sample
    figure, and must leave very little. If the out-of-sample headroom ever
    exceeds 2pp the 'the exam is the binding limit' conclusion needs redoing."""
    d = dial_ceiling()
    assert d["out_of_sample"]["headroom"] < d["in_sample"]["headroom"]
    assert d["out_of_sample"]["headroom"] < 2.0


def test_a_perfect_dial_model_still_beats_the_null():
    """The representation is not incapable — it is nearly worthless here, which
    is a different claim and the one the plan should carry."""
    d = dial_ceiling()
    assert d["out_of_sample"]["mae"] < d["out_of_sample"]["null_mae"]


def test_the_null_pool_matches_the_scored_set():
    """The restricted rows must not be scored against a null drawn from arms the
    model was never examined on -- that hands the null information the model was
    not given. Both conventions are computed here and the code must carry the
    matched one, at any arm count."""
    known = {o.arm: o.relapse_change_pct for o in KNOWN_OUTCOMES
             if o.relapse_change_pct is not None and o.arm != "untreated"}
    multi = {k: v for k, v in mechanism_groups().items() if len(v) > 1}
    scored = [a for arms in multi.values() for a in arms]

    matched_oos, wide_oos, matched_ins = [], [], []
    scored_mean = statistics.fmean(known[a] for a in scored)
    for a in scored:
        matched_oos.append(abs(known[a] - statistics.fmean(known[x] for x in scored if x != a)))
        wide_oos.append(abs(known[a] - statistics.fmean(v for k, v in known.items() if k != a)))
        matched_ins.append(abs(known[a] - scored_mean))

    d = dial_ceiling()
    assert round(d["out_of_sample"]["null_mae"], 6) == round(statistics.fmean(matched_oos), 6)
    assert round(d["multi_only"]["null_mae"], 6) == round(statistics.fmean(matched_ins), 6)
    # And the two conventions must actually differ, or this test proves nothing.
    assert round(statistics.fmean(matched_oos), 3) != round(statistics.fmean(wide_oos), 3)
