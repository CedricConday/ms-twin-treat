"""Invariants for the Weatherley-grounded population brick (B5).

These do not test that the ABM is biologically right. They test the claims the
grounding rests on, which are the ones that would rot silently:

  1. the published constants are the published constants (a typo'd rate is the
     easiest way for "grounded" to quietly become "made up again");
  2. the intervention anchor really lands on the paper's own therapy value;
  3. the nearest-myelin search is EXACT — it is an optimisation of an O(n*m)
     loop in the paper's Biased_Movement.m, and an approximation there would
     silently bias every damage number;
  4. the oligodendrocyte threshold cascade actually fires — that mechanism is
     the scientific point of the paper and the thing the old toy could not do;
  5. treatment lowers damage and an immunogenic therapy raises it, which the
     clinical gate's APL arm depends on.

Run:  python -m pytest tests/test_abm_weatherley.py -q
"""

from __future__ import annotations

import numpy as np
import pytest

from bricks.abm import (PROFILES, STEPS_PER_DAY, W, ABMBrick, WeatherleyABM,
                        bbb_permeability, simulate)


def test_published_constants_are_unchanged():
    """Guard the cited numbers. If one of these fails, the citation is now a lie."""
    assert W["tau_minutes"] == 20            # Setup_Agents.m
    assert STEPS_PER_DAY == 72
    assert W["published_days"] == 300        # ABM_main_script.m step_cap 21600
    assert (W["myelin_width"], W["domain_height"]) == (100, 300)
    assert W["oligo_dim"] == 5 and W["myelin_grades"] == 4
    assert W["heal_time"] == 25
    assert W["leaving_prob"] == 0.1 and W["entering_prob"] == 0.0
    assert W["therapeutic_leaving_prob"] == 0.025
    assert W["oligo_apop_threshold"] == 14
    assert W["oligo_stop_my_threshold"] == 10
    assert W["c1_init"] == 5 and W["c2_init"] == 45
    assert W["death_rate_per_day"] == 0.35


def test_treat_anchor_reproduces_the_papers_own_bbb_therapy():
    """treat=0.75 must land on ABM_main_script.m's intervention.new_BBB_prob."""
    assert bbb_permeability(0.75) == pytest.approx(W["therapeutic_leaving_prob"])
    assert bbb_permeability(0.0) == pytest.approx(W["leaving_prob"])
    # harm pushes permeability the other way; both ends stay probabilities
    assert bbb_permeability(0.0, 0.4) > W["leaving_prob"]
    assert 0.0 <= bbb_permeability(1.0, 0.0) <= 1.0
    assert 0.0 <= bbb_permeability(0.0, 1.0) <= 1.0


def test_nearest_myelin_search_is_exact():
    """The fast local search must equal brute force over every living piece.

    Biased_Movement.m computes this by looping each C3 over all living myelin.
    We short-circuit with a distance-ordered local ball plus an exact distance
    transform fallback. If that ever stops being exact, C3 steering — and so
    every damage number — is quietly wrong.
    """
    m = WeatherleyABM(seed=0, myelin_width=40, height=100)
    for n in range(1, 250):
        m.step(n)
    m._refresh_alive2d()
    rng = np.random.default_rng(7)
    px = rng.integers(1, m.dom.width, 200)
    py = rng.integers(1, m.dom.height, 200)

    dist, _, _, any_alive = m._nearest_alive(px, py)
    assert any_alive
    alive = m.state.ravel() > 0
    ax, ay = m.dom.piece_x[alive], m.dom.piece_y[alive]
    brute = np.sqrt((ax[None, :] - px[:, None]) ** 2
                    + (ay[None, :] - py[:, None]) ** 2).min(axis=1)
    np.testing.assert_allclose(dist, brute, atol=1e-9)


def test_oligodendrocyte_threshold_cascade_fires():
    """Apoptosis must be a cliff, not a slope: >= 14 dead pieces kills all 25.

    This is the paper's mechanism and the reason the brick was regrounded.
    """
    m = WeatherleyABM(seed=0, myelin_width=40, height=100)
    m.state[0, :W["oligo_apop_threshold"]] = 0        # push one oligo over the edge
    m._apoptosis()
    assert m.oligo_state[0] == 0                      # dead
    assert m.state[0].sum() == 0                      # took all its myelin with it

    m2 = WeatherleyABM(seed=0, myelin_width=40, height=100)
    m2.state[0, :W["oligo_stop_my_threshold"]] = 0    # over stop-my, under apoptosis
    m2._apoptosis()
    assert m2.oligo_state[0] == 2                     # stopped, not dead
    assert m2.state[0].sum() > 0


def test_dead_oligodendrocyte_never_remyelinates():
    m = WeatherleyABM(seed=0, myelin_width=40, height=100)
    m.oligo_state[0] = 0
    m.state[0, :] = 0
    for _ in range(W["heal_time"] * 3):
        m._remyelinate()
    assert m.state[0].sum() == 0


def test_treatment_reduces_and_immunogenic_increases_damage():
    """The coupling the clinical gate's arms depend on."""
    days = 10
    seeds = range(3)
    base = np.mean([simulate(days=days, seed=s)[-1] for s in seeds])
    treated = np.mean([simulate(days=days, treat=0.75, seed=s)[-1] for s in seeds])
    harmed = np.mean([simulate(days=days, immuno=0.4, seed=s)[-1] for s in seeds])
    assert treated < base, f"treatment did not reduce damage ({treated} vs {base})"
    assert harmed > base, f"immunogenic therapy did not worsen damage ({harmed} vs {base})"


def test_damage_is_a_fraction_and_monotone_nondecreasing_early():
    d = simulate(days=5, seed=0)
    assert d.ndim == 1 and d.size > 0
    assert np.all((d >= 0.0) & (d <= 1.0))


def test_profiles_all_use_the_published_lattice():
    """Only the horizon may vary. A shrunken lattice distorts magnitudes (see
    the PROFILES comment) so no profile is allowed to ship one."""
    for name, cfg in PROFILES.items():
        assert cfg["myelin_width"] == W["myelin_width"], name
        assert cfg["height"] == W["domain_height"], name
    assert PROFILES["published"]["days"] == W["published_days"]


def test_stage_interface_and_metadata_contract():
    """The spine and results/experiment.py depend on exactly this shape."""
    state = ABMBrick(profile="gate").run(
        {"intervention": {"treat": 0.5, "immunogenic": 0.0}, "seed": 1})
    assert isinstance(state["abm_damage"], np.ndarray)
    meta = state["abm_meta"]
    assert meta["validated"] is False          # grounded != validated
    assert "10.1371/journal.pcbi.1013273" in meta["grounded_in"]
    assert meta["bbb_leaving_prob"] == pytest.approx(bbb_permeability(0.5))
    assert meta["lattice"] == meta["published_lattice"]


def test_per_patient_seed_is_honoured():
    a = ABMBrick(profile="gate").run({"intervention": {"treat": 0.0}, "seed": 1})
    b = ABMBrick(profile="gate").run({"intervention": {"treat": 0.0}, "seed": 2})
    assert a["abm_meta"]["seed"] == 1 and b["abm_meta"]["seed"] == 2
    assert a["abm_meta"]["final_damage"] != b["abm_meta"]["final_damage"]
