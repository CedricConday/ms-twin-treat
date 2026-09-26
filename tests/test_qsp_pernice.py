"""Tests for the Pernice 2020 port.

1. Pin every transcribed constant against Table S1, Table S2 and Table 2, and
   the paper's statement that HD and MS differ in exactly two parameters.
2. Reproduce Figure S2: the deterministic 30-day HD-vs-MS solution, with the
   irreversibly damaged ODC count rising in MS and staying near zero in HD.
3. Guard the exam-facing readouts and the arm map.
"""

from __future__ import annotations

import numpy as np
import pytest

from backtest.clinical import KNOWN_OUTCOMES
from bricks.profiles import ENHANCE, SUPPRESS
from bricks.profiles_pernice import (
    DAC_MAX_DOSE,
    LUMPED,
    PROFILES,
    UNASSIGNABLE,
    profile_at,
    touched_points,
)
from bricks.qsp_pernice import (
    DEFAULT_READINGS,
    DIALS,
    FIG_S2_LANDMARKS,
    HEALTHY_PARAMS,
    HOURLY_TURNOVER,
    INITIAL_MARKING,
    MS_PARAMS,
    PLACES,
    READINGS_DOC,
    UNTREATED_PROFILE,
    PerniceProfile,
    calibration_run,
    landmark_report,
    simulate,
)


# --------------------------------------------------------------------------- #
# 1. Transcription
# --------------------------------------------------------------------------- #
def test_the_net_has_26_places_with_odc_unfolded():
    assert len(PLACES) == 26
    assert sum(1 for p in PLACES if p.startswith("ODC_L")) == 5


def test_table_s1_calibrated_values():
    s1 = {"pTeff_Activation": 0.015, "pTreg_Activation": 4e-04, "pTreg_Dup": 0.006,
          "pTeff_Dup": 0.04, "pTeff_KillsODC": 6e-04, "pTrkTe": 0.02, "pTekA": 6e-04,
          "pPass_BBB_treg": 0.45, "pPass_BBB_teff": 0.005, "pNKkillsTeff": 0.01,
          "pNK_prod_IFNg": 0.03, "pNK_prod_IL10": 0.045, "pIL17_BBB": 0.0115,
          "pIL10_BBB": 0.0765, "pRemyelinization": 0.01, "pIL10Consuption": 0.09,
          "pIL17Consuption": 0.03, "pIFNgConsuption": 0.05, "Cifn": 20.0, "CIL10": 10.0}
    for k, v in s1.items():
        assert HEALTHY_PARAMS[k] == v, k


def test_table_s2_fixed_values():
    assert HEALTHY_PARAMS["FromTimoREG"] == 0.317
    assert HEALTHY_PARAMS["FromTimoEFF"] == 0.296
    for k in ("NKdup", "NKDegradation", "Teff_death", "Teff_to_NLT", "Treg_death", "Treg_to_NLT"):
        assert HEALTHY_PARAMS[k] == HOURLY_TURNOVER == pytest.approx(1 / 24)
    assert HEALTHY_PARAMS["Treg_prod_IL10"] == 0.05556
    assert HEALTHY_PARAMS["Teff_prod_IL17"] == 0.00895
    assert HEALTHY_PARAMS["Teff_prod_IFNg"] == 0.0466
    assert HEALTHY_PARAMS["DACDegradation"] == 0.001444057


def test_ms_differs_from_healthy_in_exactly_two_parameters():
    diff = {k for k in HEALTHY_PARAMS if HEALTHY_PARAMS[k] != MS_PARAMS[k]}
    assert diff == {"pTeff_Activation", "pTreg_Activation"}
    assert MS_PARAMS["pTeff_Activation"] == 0.018
    assert MS_PARAMS["pTreg_Activation"] == 7e-05


def test_table_2_initial_marking():
    assert INITIAL_MARKING == {"Resting_Teff_out": 1689.0, "Resting_Treg_out": 63.0,
                               "NK_out": 30.0, "IL17_out": 8.0, "IL10_out": 13.0,
                               "IFNg_out": 42.0, "IL17_in": 1.0, "IL10_in": 1.0,
                               "IFNg_in": 1.0, "ODC_L5": 500.0}


def test_every_reading_is_documented_with_both_alternatives():
    assert set(DEFAULT_READINGS) == set(READINGS_DOC)
    for k, (a, b) in READINGS_DOC.items():
        assert a and b and a != b, k


# --------------------------------------------------------------------------- #
# 2. Reproduction: Figure S2
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def report():
    return landmark_report()


def test_figure_s2_ms_damages_odc_and_hd_does_not(report):
    ms_irr, (lo, hi), ok = report["MS"]["ODC_L1"]
    assert ok, (ms_irr, lo, hi)
    hd_irr = report["HD"]["ODC_L1"][0]
    assert hd_irr < 1.0
    assert ms_irr > 10 * max(hd_irr, 0.1)


def test_figure_s2_barrier_opens_in_ms_only(report):
    assert report["MS"]["BBB"][2] and report["HD"]["BBB"][2]
    assert report["MS"]["BBB"][0] > 3 * report["HD"]["BBB"][0]


def test_figure_s2_effectors_reach_the_cns_in_ms_only(report):
    assert report["MS"]["Teff_in"][2] and report["HD"]["Teff_in"][2]


def test_figure_s2_landmark_count_is_the_recorded_one(report):
    hits = sum(1 for lab in ("MS", "HD") for k in FIG_S2_LANDMARKS if report[lab][k][2])
    # 10 of 12 on 2026-09-26; the two misses are the Teff_out peak in both
    # configurations, recorded in the module. A drop below 10 means the port moved.
    assert hits >= 10
    assert not report["MS"]["Teff_out"][2] and not report["HD"]["Teff_out"][2]


def test_calibration_runs_stay_in_regime():
    for ms in (True, False):
        tr = calibration_run(ms)
        assert tr["in_regime"]
        assert np.isclose(tr[[p for p in PLACES if p.startswith("ODC")][0]].size, tr["t_hours"].size)
        odc_total = sum(tr[f"ODC_L{k}"] for k in range(1, 6))
        assert np.allclose(odc_total, 500.0, atol=1e-3)    # ODC tokens are conserved


# --------------------------------------------------------------------------- #
# 3. Readouts and the arm map
# --------------------------------------------------------------------------- #
def test_two_year_untreated_run_is_finite_and_damages():
    tr = simulate(UNTREATED_PROFILE, sample_hours=24.0)
    assert tr["in_regime"]
    assert tr["lesion_load"] > 0.0
    assert 0.0 < tr["irreversibly_damaged"] <= 500.0


def test_arm_set_is_the_exams():
    assert set(PROFILES) == {o.arm for o in KNOWN_OUTCOMES}


def test_every_multiplier_is_the_shared_stub():
    for name, prof in PROFILES.items():
        for k, v in prof.multipliers.items():
            assert k in DIALS
            assert v in (SUPPRESS, ENHANCE), f"{name}.{k} = {v} is not the stub"


def test_lumped_and_unassignable_are_disjoint():
    assert not set(LUMPED) & set(UNASSIGNABLE)
    for name in UNASSIGNABLE:
        assert PROFILES[name].is_untreated, name


def test_the_port_separates_what_velez_lumped():
    pats = {a: touched_points(PROFILES[a]) for a in ("natalizumab", "alemtuzumab", "Tovaxin",
                                                      "ocrelizumab", "IFN-beta")}
    assert len(set(pats.values())) == 5, pats
    assert touched_points(PROFILES["fingolimod"]) == touched_points(PROFILES["natalizumab"])


def test_profile_at_moves_size_only_and_dac_dose_scales():
    for arm in ("natalizumab", "IFN-gamma", "lenercept", "alemtuzumab"):
        base = PROFILES[arm]
        for s in (0.0, 0.3, 0.9):
            swept = profile_at(arm, s)
            assert touched_points(swept) == (touched_points(base) if s > 0 else ())
            for k, v in base.multipliers.items():
                expect = (1.0 - s) if v < 1.0 else (1.0 + s)
                assert swept.multipliers[k] == pytest.approx(expect)
    d = profile_at("daclizumab", 0.5)
    assert d.dac_dose == pytest.approx(0.5 * DAC_MAX_DOSE)
    assert profile_at("daclizumab", 0.0).is_untreated


def test_profiles_reject_unknown_dials():
    with pytest.raises(ValueError):
        PerniceProfile(label="x", multipliers={"gamma_E": 0.5})
