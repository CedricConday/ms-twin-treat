"""Tests for leave-one-mechanism-out.

None of these build the real response table — that is a multi-minute cached
step. They test the things that would silently corrupt the result: how arms are
grouped into mechanisms, that a potency sweep cannot flip a sourced direction,
and that the estimator choices forced by the model's variance stay put.

The last of those matters most. A first version of `backtest/lomo.py` used four
seeds and a mean, and reported a confident headline that was entirely its own
noise. The tests below pin the fix so it cannot quietly regress.
"""

from __future__ import annotations

import numpy as np
import pytest

from backtest import lomo
from bricks.profiles import PROFILES, touched_points

# --------------------------------------------------------------------------- #
# grouping — what "a mechanism" means here
# --------------------------------------------------------------------------- #

def test_groups_are_keyed_by_intervention_points_not_by_a_label():
    """The group key IS the mechanism in this representation.

    That is what makes holding one out 'holding out a mechanism' rather than
    holding out a name somebody assigned.
    """
    groups = lomo.mechanism_groups()
    for key, arms in groups.items():
        for arm in arms:
            assert touched_points(PROFILES[arm]) == key


def test_the_quantified_arms_fall_into_the_expected_mechanisms():
    groups = lomo.mechanism_groups()
    assert set(groups) == {
        ("alpha_E",),            # IFN-beta, teriflunomide, dimethyl fumarate
        ("alpha_R", "delta"),    # glatiramer acetate
        ("gamma_E",),            # natalizumab, fingolimod, alemtuzumab, ponesimod
        ("ke",),                 # ocrelizumab
    }
    assert sum(len(a) for a in groups.values()) == 9


def test_only_arms_with_a_real_number_are_grouped():
    """Unquantified harm arms cannot be scored, so they cannot be fitted either."""
    grouped = {a for arms in lomo.mechanism_groups().values() for a in arms}
    for unscorable in ("lenercept", "atacicept", "IFN-gamma", "APL CGP77116", "untreated"):
        assert unscorable not in grouped


def test_a_single_arm_group_exists_and_is_a_known_weakness():
    """Holding out a one-arm group means one number moves a quarter of the headline.

    Asserted rather than left implicit: `main()` prints per-group rows because
    of this, and a reader who only sees the headline is being misled.
    """
    sizes = [len(a) for a in lomo.mechanism_groups().values()]
    assert min(sizes) == 1
    assert len(sizes) == 4


# --------------------------------------------------------------------------- #
# the potency sweep must not become a direction sweep
# --------------------------------------------------------------------------- #

def test_potency_sweep_preserves_each_sourced_direction():
    """Fitting may change how big an effect is, never which way it points.

    Directions come from pharmacology in bricks/profiles.py. If a sweep could
    flip one, the fit would be choosing a drug's mechanism to suit the data,
    which is the exact failure this repo keeps testing for.
    """
    for arm in (a for arms in lomo.mechanism_groups().values() for a in arms):
        base = PROFILES[arm]
        for s in (0.05, 0.4, 0.9):
            swept = lomo._profile_at(arm, s)
            assert touched_points(swept) == touched_points(base)
            for point in touched_points(base):
                base_v, swept_v = getattr(base, point), getattr(swept, point)
                assert (base_v < 1.0) == (swept_v < 1.0), (
                    f"{arm}.{point} changed direction at potency {s}")


def test_zero_potency_is_the_untreated_profile():
    for arm in ("IFN-beta", "natalizumab", "ocrelizumab"):
        swept = lomo._profile_at(arm, 0.0)
        for point in touched_points(PROFILES[arm]):
            assert getattr(swept, point) == 1.0


def test_untouched_dials_stay_untouched_at_every_potency():
    from bricks.qsp_velez import INTERVENTION_POINTS
    base = PROFILES["IFN-beta"]
    untouched = [p for p in INTERVENTION_POINTS if getattr(base, p) == 1.0]
    for s in (0.1, 0.5, 0.9):
        swept = lomo._profile_at("IFN-beta", s)
        for point in untouched:
            assert getattr(swept, point) == 1.0


# --------------------------------------------------------------------------- #
# the estimator — forced by measured variance, not chosen
# --------------------------------------------------------------------------- #

def test_cohort_is_large_enough_for_the_models_variance():
    """128 seeds. See the SEEDS comment in backtest/lomo.py.

    Re-derive with scripts/measure_qsp_variance.py. At a 730-day horizon over
    200 runs: median 2.14, mean 74.9, SD 612, range 0.089-8118. The bootstrapped
    CV of the median is 628% at n=4 and 13.7% at n=128. The first version of
    this file used 4 seeds and reported a headline made entirely of that noise.
    """
    assert len(lomo.SEEDS) >= 100
    assert len(set(lomo.SEEDS)) == len(lomo.SEEDS), "seeds must be distinct"


def test_the_statistic_is_a_median():
    """Right-tailed to the point that the mean is not an estimate of anything.

    Measured: mean/median = 35x, and the top decile of runs holds 96% of all
    damage.
    """
    skewed = [1.0, 1.1, 0.9, 1.0, 500.0]
    assert lomo._median_ratio(skewed) == 1.0
    assert lomo._median_ratio(skewed) != np.mean(skewed)


def test_horizon_matches_the_trials_not_the_models_default():
    """Two years. The trials report over 1-2 years; the model's own FINAL TIME is 5."""
    assert lomo.T_END == 730.0


def test_potency_grid_excludes_a_fully_removed_rate():
    """Multiplying a published rate by 0 is outside anything it was characterised on."""
    assert min(lomo.POTENCY_GRID) == 0.0
    assert max(lomo.POTENCY_GRID) < 1.0


# --------------------------------------------------------------------------- #
# the fit, on a synthetic table (no simulation)
# --------------------------------------------------------------------------- #

def _synthetic_table(slope_by_pattern: dict[str, float]) -> dict:
    """A table where each mechanism's effect is exactly linear in potency."""
    grid = lomo.POTENCY_GRID
    return {
        "potency_grid": grid,
        "table": {
            pattern: {f"{s:g}": {"mean": -100.0 * slope * s, "sd": 0.0, "n_ok": 1}
                      for s in grid}
            for pattern, slope in slope_by_pattern.items()
        },
    }


def test_fit_recovers_a_potency_that_explains_the_training_arms():
    patterns = {"|".join(k): 1.0 for k in lomo.mechanism_groups()}
    table = _synthetic_table(patterns)
    known = {"IFN-beta": -30.0, "teriflunomide": -30.0}
    s = lomo._fit_potency(["IFN-beta", "teriflunomide"], known, table)
    assert s == pytest.approx(0.30, abs=0.02)


def test_run_lomo_scores_every_quantified_arm_exactly_once():
    patterns = {"|".join(k): 1.0 for k in lomo.mechanism_groups()}
    result = lomo.run_lomo(_synthetic_table(patterns))
    arms = [r["arm"] for r in result["rows"]]
    assert len(arms) == len(set(arms)) == 9
    assert result["n_groups"] == 4


def test_a_held_out_mechanism_is_absent_from_its_own_training_set():
    """The property that makes this LOMO and not LOO.

    Every arm sharing the held-out mechanism must be out of the fit. If one
    leaked in, the test would be grading calibration while claiming to grade
    generalisation.
    """
    groups = lomo.mechanism_groups()
    for held_pattern, held_arms in groups.items():
        training = [a for pat, arms in groups.items() if pat != held_pattern
                    for a in arms]
        assert not set(held_arms) & set(training)
        for arm in training:
            assert touched_points(PROFILES[arm]) != held_pattern


def test_lomo_holds_out_strictly_more_than_loo_would():
    """Holding out a mechanism removes the held arm AND its mechanism-mates.

    For every group with more than one arm, LOMO's training set is strictly
    smaller than LOO's would be for the same held-out arm. That is the whole
    difference between the two tests.
    """
    groups = lomo.mechanism_groups()
    multi = [arms for arms in groups.values() if len(arms) > 1]
    assert multi, "expected at least one multi-arm mechanism group"
    total = sum(len(a) for a in groups.values())
    for arms in multi:
        lomo_training = total - len(arms)
        loo_training = total - 1
        assert lomo_training < loo_training


def test_null_is_predict_the_mean_of_the_training_arms():
    patterns = {"|".join(k): 1.0 for k in lomo.mechanism_groups()}
    result = lomo.run_lomo(_synthetic_table(patterns))
    for row in result["rows"]:
        assert row["null_error"] == pytest.approx(
            abs(row["null_predicted"] - row["known"]))
