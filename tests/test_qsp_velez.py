"""Tests for the Vélez de Mendizábal 2011 QSP port.

Three jobs, in order of what they protect:

1. **Pin every transcribed constant** against the authors' Additional file 2.
   These tests exist so that a later "tuning" commit cannot quietly move a rate
   that came from a paper. If one fails, either the port drifted or someone
   decided to fit the model — and the second one needs to be an explicit,
   argued change to the docstring, not a number nudge.
2. **Reproduce the paper's Figure 3 claim** — alpha_R is the health/autoimmunity
   axis. This is the only test here that checks the port against a published
   RESULT rather than a published number.
3. **Guard the two traps found while porting**: the relapse detector's baseline,
   and the model's lack of a carrying capacity on E.
"""

from __future__ import annotations

import numpy as np
import pytest

from bricks.qsp_velez import (
    INTERVENTION_POINTS,
    UNTREATED_PROFILE,
    VELEZ_DT,
    VELEZ_PARAMS,
    VELEZ_T_END,
    VELEZ_Y0,
    MechanismProfile,
    VelezQSPBrick,
    baseline_from,
    relapse_events,
    simulate,
)

# A short horizon keeps the suite fast. Long enough for the cross-regulation
# loop to settle and for several stochastic immune events to land.
SHORT = dict(t_end=365.0, seed=7)


# --------------------------------------------------------------------------- #
# 1. transcribed constants — every one of these is a citation, not a choice
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("key,value", [
    ("delta", 1.0),      # MDL AntigenPresentation = 1
    ("beta", 0.01),      # MDL TeAnergy = Resting Te*0.01
    ("eta", 0.01),       # MDL TeMemory = Activated Te*0.01
    ("alpha_E", 2.0),    # MDL Max Te Proliferation Rate = 2
    ("alpha_R", 0.25),   # MDL Max Tr activation and proliferation Rate = 0.25
    ("gamma_E", 0.2),    # MDL Max Te Death Anergy Migration Rate = 0.2
    ("gamma_R", 0.2),    # MDL Tr Death Anergy Migration Rate = 0.2
    ("ke", 1000.0),      # MDL ke = 1000
    ("kr", 200.0),       # MDL kr = 200
    ("h", 5.0),          # MDL h = 5
    ("d1", 1.0),         # MDL d1 = 1
    ("d2", 0.002),       # MDL d2 = 0.002  (Table 1 says 0.02 — documented discrepancy)
    ("r", 0.1),          # MDL r = 0.1
    ("A", 22800.0),      # MDL A = 22800
    ("n", 2.0),          # MDL n = 2
])
def test_parameter_matches_published_model(key, value):
    assert VELEZ_PARAMS[key] == value, (
        f"{key} no longer matches Vélez de Mendizábal 2011 Additional file 2. "
        "These are transcribed, not tuned — see the module docstring."
    )


@pytest.mark.parametrize("key,value", [
    ("Er", 7.5), ("Rr", 2.4), ("E", 1000.0), ("R", 200.0), ("l", 0.0), ("L", 0.0),
])
def test_initial_conditions_match_published_model(key, value):
    assert VELEZ_Y0[key] == value


def test_simulation_control_matches_published_model():
    # MDL .Control: TIME STEP = 0.1 Day, FINAL TIME = 1825 Day
    assert VELEZ_DT == 0.1
    assert VELEZ_T_END == 1825.0


def test_d2_discrepancy_is_the_mdl_value_not_the_table_value():
    """Documented DISCREPANCY 1. Table 1 prints 0.02; the runnable model says 0.002.

    Pinned explicitly so that "fixing" it to the table value has to be a
    deliberate, argued change rather than a plausible-looking typo repair.
    """
    assert VELEZ_PARAMS["d2"] == 0.002
    assert VELEZ_PARAMS["d2"] != 0.02


# --------------------------------------------------------------------------- #
# 2. the paper's Figure 3 result
# --------------------------------------------------------------------------- #

def test_alpha_R_is_the_health_autoimmunity_axis():
    """Vélez de Mendizábal 2011, Figure 3 / Results.

    "by decreasing the maximum activation and proliferation rate of Treg
     (alpha_R) ... immune homeostasis was lost and spontaneous immune responses
     emerged in the absence of infectious agents"

    So raising alpha_R from the MDL's operating point must move the system
    toward homeostasis: lower effector load and less tissue damage. This is a
    check against a published RESULT, and it is the test that would catch a
    port whose numbers are all individually right and whose dynamics are wrong.
    """
    damages, mean_E = [], []
    for mult in (1.0, 2.0, 4.0, 8.0):
        traj = simulate(MechanismProfile(label=f"aR x{mult}", alpha_R=mult), **SHORT)
        assert traj["in_regime"]
        damages.append(float(traj["total_damage"][-1]))
        mean_E.append(float(np.mean(traj["E"])))

    assert damages == sorted(damages, reverse=True), (
        f"damage should fall monotonically as alpha_R rises, got {damages}")
    assert mean_E == sorted(mean_E, reverse=True), (
        f"effector load should fall monotonically as alpha_R rises, got {mean_E}")
    # Not a marginal effect: the paper describes two qualitatively different regimes.
    assert damages[0] > 10 * damages[-1]


def test_relapses_emerge_rather_than_being_scheduled():
    """The reason this model replaces the toy: relapses are not put in by hand.

    Nothing in the equations schedules an event. Peaks come from stochastic
    immune input hitting a regulated loop, so different infection histories
    must give different relapse timings.
    """
    a = simulate(UNTREATED_PROFILE, t_end=730.0, seed=1)
    b = simulate(UNTREATED_PROFILE, t_end=730.0, seed=2)

    ta, tb = relapse_events(a), relapse_events(b)
    assert ta, "untreated autoimmune configuration should produce effector peaks"
    assert ta != tb, "relapse timings must depend on the infection history, not a schedule"


def test_same_seed_is_the_same_infection_history():
    """Arm-vs-arm comparison is only meaningful on a shared stochastic history."""
    a = simulate(UNTREATED_PROFILE, **SHORT)
    b = simulate(UNTREATED_PROFILE, **SHORT)
    assert np.array_equal(a["E"], b["E"])


# --------------------------------------------------------------------------- #
# 3. the two traps
# --------------------------------------------------------------------------- #

def test_relapse_baseline_must_come_from_the_untreated_arm():
    """The trap that made the detector blind across arms.

    A treated arm's own median E falls with the treatment, so a self-referenced
    threshold slides down with it and reports business as usual for a cohort
    whose effector load collapsed. Scored against the UNTREATED arm's baseline,
    a strongly Treg-supporting arm must show strictly fewer relapses.
    """
    untreated = simulate(UNTREATED_PROFILE, t_end=730.0, seed=3)
    treated = simulate(MechanismProfile(label="Treg support", alpha_R=4.0),
                       t_end=730.0, seed=3)
    ref = baseline_from(untreated)

    shared = relapse_events(treated, baseline=ref)
    self_referenced = relapse_events(treated)

    assert len(shared) < len(relapse_events(untreated, baseline=ref))
    assert len(shared) < len(self_referenced), (
        "self-referenced thresholding hides the treatment effect — that is the "
        "whole reason `baseline` exists"
    )


def test_stripping_regulation_leaves_the_model_regime_and_says_so():
    """The published model has NO carrying capacity on E.

    Remove enough regulation and dE/dt -> (alpha_E - eta)*E, i.e. unbounded
    growth. That is the model's property, not a porting bug, but the resulting
    damage number is meaningless and must never be ranked against an in-regime
    arm. The run has to report that rather than return a very large number.
    """
    traj = simulate(MechanismProfile(label="regulation stripped", gamma_R=4.0),
                    t_end=730.0, seed=5)
    assert traj["in_regime"] is False
    assert traj["left_regime_at"] is not None
    assert traj["left_regime_at"] > 0.0


def test_in_regime_arms_are_reported_as_such():
    traj = simulate(UNTREATED_PROFILE, **SHORT)
    assert traj["in_regime"] is True
    assert traj["left_regime_at"] is None


# --------------------------------------------------------------------------- #
# interface + guards
# --------------------------------------------------------------------------- #

def test_untreated_profile_touches_nothing():
    assert UNTREATED_PROFILE.is_untreated
    for point in INTERVENTION_POINTS:
        assert getattr(UNTREATED_PROFILE, point) == 1.0


def test_negative_multiplier_is_rejected():
    with pytest.raises(ValueError):
        MechanismProfile(label="nonsense", alpha_E=-1.0)


def test_unknown_damage_form_is_rejected():
    with pytest.raises(ValueError):
        simulate(UNTREATED_PROFILE, damage_form="whatever", **SHORT)


def test_damage_forms_differ_only_by_scale():
    """DISCREPANCY 3: the two readings differ by a constant factor, nothing more.

    This is why the resolution could not have changed any arm-vs-arm ratio the
    clinical gate scores — and why it still had to be resolved before any
    absolute damage number is quoted.
    """
    paper = simulate(UNTREATED_PROFILE, damage_form="paper", **SHORT)
    literal = simulate(UNTREATED_PROFILE, damage_form="mdl_literal", **SHORT)
    ratio = literal["total_damage"][-1] / paper["total_damage"][-1]
    assert ratio == pytest.approx(VELEZ_PARAMS["A"], rel=1e-6)


def test_brick_writes_the_keys_the_spine_expects():
    brick = VelezQSPBrick(t_end=365.0)
    state = brick.run({"seed": 11})
    assert set(state["qsp_traj"]) >= {"t", "E", "R", "l", "L", "total_damage"}
    assert state["qsp_traj"]["validated"] is False
    assert len(state["qsp_damage"]) == len(state["qsp_traj"]["t"])
    assert isinstance(state["qsp_relapses"], list)


def test_brick_reads_a_mechanism_profile_from_state():
    brick = VelezQSPBrick(t_end=365.0)
    state = brick.run({
        "seed": 11,
        "mechanism_profile": {"label": "Treg support", "alpha_R": 4.0},
    })
    assert state["qsp_traj"]["profile"]["alpha_R"] == 4.0


def test_scalar_shim_is_labelled_as_a_shim():
    """The legacy fallback must never be mistaken for a grounded mechanism."""
    brick = VelezQSPBrick(t_end=365.0)
    state = brick.run({
        "seed": 11,
        "intervention": {"name": "legacy", "treat": 0.5, "immunogenic": 0.0},
    })
    profile = state["qsp_traj"]["profile"]
    assert profile["alpha_E"] == 0.5
    assert "SHIM" in profile["source"]


def test_damage_is_driven_by_peaks_not_by_median_effector_load():
    """The structural finding that explains the whole gamma_E failure.

    Damage is driven by (E/a)^2, so it is set by effector peak excursions. A
    treatment can barely move the median effector load and still change damage
    by an order of magnitude — which is why "it lowers effector numbers" is not
    a safe proxy for "it helps" in this model.

    Measured over 48 histories at 730 days: median E lands at 1020-1126 for
    untreated, damped and killing arms alike, while peak E spans 30k-192k.
    """
    seeds = range(24)
    out = {}
    for label, profile in (
        ("untreated", UNTREATED_PROFILE),
        ("damped", MechanismProfile(label="damped", alpha_E=0.5)),
        ("killing", MechanismProfile(label="killing", gamma_E=1.5)),
    ):
        med, peak, dmg = [], [], []
        for seed in seeds:
            traj = simulate(profile, t_end=730.0, seed=seed)
            if not traj["in_regime"]:
                continue
            med.append(float(np.median(traj["E"])))
            peak.append(float(np.max(traj["E"])))
            dmg.append(float(traj["total_damage"][-1]))
        out[label] = (float(np.median(med)), float(np.median(peak)), float(np.median(dmg)))

    medians = [v[0] for v in out.values()]
    peaks = [v[1] for v in out.values()]
    assert max(medians) / min(medians) < 1.5, (
        f"median effector load should barely move, got {medians}")
    assert max(peaks) / min(peaks) > 3.0, (
        f"peak effector load should be what differs, got {peaks}")


def test_killing_effectors_is_self_defeating_in_this_model():
    """Raising gamma_E raises damage. Structural, not a calibration gap.

    Effectors recruit their own regulators, so removing them releases the
    proliferation brake and the population rebounds into a larger excursion.
    An additive regulation-independent depletion term was tested as an
    extension and did not fix this either; it was removed rather than tuned.

    This is why the depleting/sequestering class cannot come out beneficial
    here, and why backtest/clinical_velez.py scores 5/13.
    """
    seeds = range(24)

    def median_damage(profile):
        vals = [float(simulate(profile, t_end=730.0, seed=s)["total_damage"][-1])
                for s in seeds
                if simulate(profile, t_end=730.0, seed=s)["in_regime"]]
        return float(np.median(vals))

    untreated = median_damage(UNTREATED_PROFILE)
    killing = median_damage(MechanismProfile(label="killing", gamma_E=1.5))
    damping = median_damage(MechanismProfile(label="damping", alpha_E=0.5))

    assert killing > untreated, "raising gamma_E should RAISE damage in this model"
    assert damping < untreated, "damping proliferation should lower it"
