"""Tests for the Jenner 2026 two-equation port and its arm map.

1. Pin the transcribed parameter sets and the paper's analytic landmarks
   (branch point, Hopf locus, disease equilibrium) against the numbers the
   paper prints in its captions.
2. Reproduce the paper's result: raising phi carries the system from the
   healthy state, through a stable disease state, into limit cycles (section
   4, figures 2 and 3), with the period the captions state.
3. Guard the arm map: every non-unit multiplier is the shared stub, the
   pattern key is signed, and the arm set is the exam's.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from backtest.clinical import KNOWN_OUTCOMES
from bricks.profiles import ENHANCE, SUPPRESS
from bricks.profiles_minimal import LUMPED, PROFILES, UNASSIGNABLE, profile_at, touched_points
from bricks.qsp_minimal import (
    DT,
    FIG2D_PARAMS,
    FIG2D_PERIOD_MONTHS,
    FIG3_PARAMS,
    FIG3_PHI_HOPF,
    FIG3_PHI_STAR,
    RANGES,
    RATES,
    UNTREATED_PARAMS,
    UNTREATED_PROFILE,
    Y0,
    MinimalProfile,
    compare,
    disease_equilibrium,
    integrate,
    oscillation_amplitude,
    oscillation_period,
    phi_hopf,
    phi_star,
    simulate,
)


# --------------------------------------------------------------------------- #
# 1. Transcription
# --------------------------------------------------------------------------- #
def test_untreated_arm_is_the_papers_figure_2d_patient():
    assert UNTREATED_PARAMS == {"r": 0.5, "phi": 0.7, "eta": 0.5, "delta": 0.2}
    for k, v in UNTREATED_PARAMS.items():
        lo, hi = RANGES[k]
        assert lo < v <= hi, k


def test_landmarks_match_the_figure_3_caption():
    p = dict(FIG3_PARAMS, phi=1.0)
    assert phi_star(p) == pytest.approx(FIG3_PHI_STAR, abs=1e-9)
    assert phi_hopf(p) == pytest.approx(FIG3_PHI_HOPF, abs=0.005)


def test_figure_2_landmarks_bracket_the_three_panels():
    p = dict(FIG2D_PARAMS)
    bp, hb = phi_star(p), phi_hopf(p)
    assert bp == pytest.approx(0.3) and hb == pytest.approx(0.6)
    assert 0.2 < bp < 0.45 < hb < 0.7   # panels B, C, D of figure 2


def test_disease_equilibrium_is_where_the_flow_settles_below_the_hopf():
    p = dict(FIG2D_PARAMS, phi=0.45)     # figure 2C
    tr = integrate(p, Y0, 600.0)
    m, i = disease_equilibrium(p)
    assert tr["in_regime"]
    assert tr["M"][-1] == pytest.approx(m, abs=2e-3)
    assert tr["I"][-1] == pytest.approx(i, abs=2e-3)
    assert oscillation_amplitude(tr) < 1e-3


# --------------------------------------------------------------------------- #
# 2. Reproduction: phi carries the system through the Hopf
# --------------------------------------------------------------------------- #
def test_figure_2b_healthy_state_is_recovered():
    tr = integrate(dict(FIG2D_PARAMS, phi=0.2), Y0, 300.0)
    assert tr["M"][-1] == pytest.approx(1.0, abs=1e-3)
    assert tr["I"][-1] == pytest.approx(0.0, abs=1e-6)


def test_figure_2d_limit_cycle_has_the_captions_period():
    tr = integrate(FIG2D_PARAMS, Y0, 600.0)
    assert oscillation_amplitude(tr) > 0.3
    # The caption says "a period of t = 30"; read off a plot. The integrated
    # value is 28.3 months; a 10% band holds the reproduction without pretending
    # the caption's figure carries more digits than it does.
    assert oscillation_period(tr) == pytest.approx(FIG2D_PERIOD_MONTHS, rel=0.10)


def test_figure_3_hopf_separates_a_fixed_point_from_a_cycle():
    below = integrate(dict(FIG3_PARAMS, phi=3.0), Y0, 600.0)
    above = integrate(dict(FIG3_PARAMS, phi=7.0), Y0, 600.0)
    assert oscillation_amplitude(below) < 1e-3
    assert oscillation_amplitude(above) > 0.1
    # caption: periods "from a minimum of about T = 10 months to T = 20 months"
    assert 10.0 <= oscillation_period(above) <= 20.0


def test_amplitude_grows_with_phi_past_the_hopf():
    amps = [oscillation_amplitude(integrate(dict(FIG2D_PARAMS, phi=phi), Y0, 600.0))
            for phi in (0.65, 0.7, 0.8, 1.0)]
    assert all(a < b for a, b in zip(amps[:-1], amps[1:], strict=True))


def test_step_size_is_converged():
    coarse = simulate(UNTREATED_PROFILE)
    fine = simulate(UNTREATED_PROFILE, dt=DT / 2)
    assert oscillation_period(fine) == pytest.approx(oscillation_period(coarse), rel=1e-3)
    assert fine["integrated_I"] == pytest.approx(coarse["integrated_I"], rel=1e-4)


def test_untreated_patient_relapses_and_stays_in_regime():
    c = compare(UNTREATED_PROFILE)
    assert c["in_regime"]
    assert c["lesion_ratio"] == pytest.approx(1.0)
    assert c["untreated_relapses_per_year"] > 0.0


def test_removing_the_disease_strength_entirely_heals():
    c = compare(MinimalProfile(label="phi off", phi=0.05))
    assert c["in_regime"]
    assert c["lesion_ratio"] < 0.05
    assert c["relapses_per_year"] == 0.0
    assert c["mean_M"] > 0.95


# --------------------------------------------------------------------------- #
# 3. The arm map
# --------------------------------------------------------------------------- #
def test_arm_set_is_the_exams():
    assert set(PROFILES) == {o.arm for o in KNOWN_OUTCOMES}


def test_every_magnitude_is_the_shared_stub():
    for name, prof in PROFILES.items():
        for rate in RATES:
            v = getattr(prof, rate)
            assert v in (1.0, SUPPRESS, ENHANCE), f"{name}.{rate} = {v} is not the stub"


def test_lumped_and_unassignable_are_disjoint_and_complete():
    assert not set(LUMPED) & set(UNASSIGNABLE)
    for name in LUMPED:
        assert touched_points(PROFILES[name]) == ("phi-",), name
    for name in UNASSIGNABLE:
        assert touched_points(PROFILES[name]) == (), name
        assert PROFILES[name].is_untreated
    rest = set(PROFILES) - set(LUMPED) - set(UNASSIGNABLE) - {"untreated"}
    assert rest == {"IFN-gamma", "APL CGP77116", "lenercept"}


def test_pattern_key_carries_the_sign():
    assert touched_points(PROFILES["IFN-gamma"]) == ("phi+",)
    assert touched_points(PROFILES["natalizumab"]) == ("phi-",)
    assert touched_points(PROFILES["lenercept"]) == ("r-", "phi-")


def test_profile_at_keeps_direction_and_moves_only_size():
    for arm in ("natalizumab", "IFN-gamma", "lenercept", "abatacept"):
        base = PROFILES[arm]
        for s in (0.0, 0.3, 0.9):
            swept = profile_at(arm, s)
            assert touched_points(swept) == (touched_points(base) if s > 0 else ())
            for rate in RATES:
                v = getattr(base, rate)
                expect = 1.0 if v == 1.0 else ((1.0 - s) if v < 1.0 else (1.0 + s))
                assert getattr(swept, rate) == pytest.approx(expect)


def test_profiles_reject_unknown_rates_and_negative_values():
    with pytest.raises(TypeError):
        MinimalProfile(label="x", gamma_E=0.5)  # a Velez dial, not one of these
    with pytest.raises(ValueError):
        MinimalProfile(label="x", phi=-0.1)
    with pytest.raises(ValueError):
        MinimalProfile(label="x", r=math.nan)


def test_integrated_readouts_are_finite_for_every_arm_on_the_grid():
    u = simulate(UNTREATED_PROFILE)
    for arm in PROFILES:
        for s in (0.0, 0.5, 0.95):
            c = compare(profile_at(arm, s), u)
            assert c["in_regime"], (arm, s)
            assert np.isfinite(c["lesion_ratio"]), (arm, s)
