"""Tests for LOMO under the capacity extension.

None of these build a response table — that is the multi-minute cached step. They
pin the two things that would make the extended headline a lie: that the capacity
is fitted inside the fold rather than chosen by looking at the answer, and that
the published (uncapped) model stays in the search so the extension has to earn
its place.
"""

from __future__ import annotations

import numpy as np

from backtest import lomo, lomo_capacity
from bricks.profiles import PROFILES, touched_points
from bricks.qsp_velez import MechanismProfile

# --------------------------------------------------------------------------- #
# the capacity is a different model, and is cached as one
# --------------------------------------------------------------------------- #

def test_the_capacity_table_never_overwrites_the_transcription_curve():
    assert lomo.cache_path(None) == lomo.CACHE
    capped = lomo.cache_path(2000.0)
    assert capped != lomo.CACHE
    assert "K2000" in capped.name


def test_the_published_model_stays_in_the_search():
    """None must be a candidate capacity, or the extension cannot lose."""
    assert None in lomo_capacity.K_GRID


def test_the_capacity_reaches_the_solver():
    """A plumbing test: the kwarg must change the simulation, not be swallowed."""
    depleting = MechanismProfile(label="gamma_E x3", source="test", gamma_E=3.0)
    uncapped = lomo._damage(depleting, seed=0)
    capped = lomo._damage(depleting, seed=0, carrying_capacity=2000.0)
    assert uncapped is None or capped is None or uncapped != capped


# --------------------------------------------------------------------------- #
# the fit must be blind to the held-out mechanism
# --------------------------------------------------------------------------- #

def _flat_table(value_by_pattern: dict[str, float]) -> dict:
    """A response table whose every potency predicts the same per-pattern value."""
    grid = [0.0, 0.5, 0.9]
    return {
        "potency_grid": grid,
        "table": {pattern: {f"{s:g}": {"mean": v} for s in grid}
                  for pattern, v in value_by_pattern.items()},
    }


def test_the_fit_ignores_the_held_out_arms_outcome():
    """Fitting on N-1 arms must give the same (K, s) whatever the Nth arm did."""
    groups = lomo.mechanism_groups()
    patterns = ["|".join(p) for p in groups]
    tables = {
        None: _flat_table({p: -10.0 for p in patterns}),
        2000.0: _flat_table({p: -30.0 for p in patterns}),
    }
    held = next(iter(groups))
    training = [a for pat, arms in groups.items() if pat != held for a in arms]

    known = {a: -30.0 for a in training}
    for outcome in (0.0, -90.0, +50.0):
        with_held = dict(known)
        for arm in groups[held]:
            with_held[arm] = outcome
        assert lomo_capacity._fit(training, with_held, tables)[:2] == (2000.0, 0.0)


def test_predictions_come_from_the_fitted_capacitys_own_table():
    """Mixing a fold's potency with another capacity's curve would be a silent bug."""
    groups = lomo.mechanism_groups()
    arm = next(a for arms in groups.values() for a in arms)
    key = "|".join(touched_points(PROFILES[arm]))
    tables = {None: _flat_table({key: -10.0}), 2000.0: _flat_table({key: -30.0})}
    assert np.isclose(lomo._predicted(arm, 0.4, tables[2000.0]), -30.0)
    assert np.isclose(lomo._predicted(arm, 0.4, tables[None]), -10.0)
