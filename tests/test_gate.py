"""The device's contract, including the one property that matters most:

    no input can produce PASS while the certificate does not unlock it.

Tests that only checked "decide returns something sensible" would pass equally
well against a device that had hard-coded ABSTAIN, so the PASS path is exercised
here with a synthetic certificate. If it were dead code, this file would fail --
which is the point: the day the real measurement inverts, the device must change
its answer with nobody editing it.
"""

from __future__ import annotations

import dataclasses
import json

import pytest

from bricks.qsp_velez import MechanismProfile
from gate.criterion import CRITERION
from gate.device import Verdict, VerdictKind, decide, decide_all
from gate.evidence import (
    MEASURED,
    UNMEASURED,
    EvidenceCertificate,
    ScorerEvidence,
    _paired_bootstrap,
    certify,
)
from gate.provenance import (
    PUBLISHED,
    ModelMismatch,
    ModelProvenance,
    combined,
    from_simulator,
    from_table,
)
from screen.kill_filter import screen


@pytest.fixture(scope="module")
def survivor() -> MechanismProfile:
    """A candidate pattern nothing currently occupies, found at test time.

    Hardcoding one breaks whenever an arm is wired onto that dial -- `delta-`
    was free until abatacept landed on `delta`, at which point the device
    correctly began returning DEGENERATE and three tests failed for the right
    reason. Asking the screen which patterns survive today keeps the fixture
    valid across arm-set changes instead of deferring the same break.
    """
    from screen.kill_filter import enumerate_candidates

    for result in screen(enumerate_candidates(max_points=2)):
        if result.survived:
            return result.profile
    pytest.skip("no surviving candidate on this model; the device has nothing to abstain on")


@pytest.fixture(scope="module")
def doomed() -> MechanismProfile:
    """A candidate the model kills, found the same way and for the same reason."""
    from screen.kill_filter import enumerate_candidates

    for result in screen(enumerate_candidates(max_points=1)):
        if not result.survived:
            return result.profile
    pytest.skip("nothing is killed on this model; the kill filters have stopped working")


def _passing_scorer(name: str) -> ScorerEvidence:
    """A scorer that clears every clause of the criterion, by construction."""
    return ScorerEvidence(
        name=name, status=MEASURED, n_units=12,
        mae=5.0, null_mae=15.0,            # 67% relative margin, well past 20%
        ci_low=-14.0, ci_high=-6.0,        # CI entirely below zero
        mae_subset=5.0, null_mae_subset=15.0,
        detail="synthetic certificate for the PASS-path test",
    )


def _passing_certificate() -> EvidenceCertificate:
    return EvidenceCertificate(scorers={
        "leave-one-arm-out": _passing_scorer("leave-one-arm-out"),
        "leave-one-mechanism-out": _passing_scorer("leave-one-mechanism-out"),
    })


# --- the criterion itself ------------------------------------------------


def test_margin_requires_the_full_fifth():
    assert CRITERION.margin_met(mae=12.0, null_mae=15.0)        # exactly 20%
    assert not CRITERION.margin_met(mae=12.1, null_mae=15.0)    # 19.3%
    assert not CRITERION.margin_met(mae=15.0, null_mae=15.0)    # ties lose


def test_a_null_with_no_error_cannot_be_beaten():
    """Guards a degenerate arm set dividing by zero into a free pass."""
    assert not CRITERION.margin_met(mae=0.0, null_mae=0.0)


def test_criterion_is_frozen_and_dated():
    assert CRITERION.fixed_on == "2026-09-20"
    with pytest.raises(dataclasses.FrozenInstanceError):
        CRITERION.min_relative_margin = 0.0  # type: ignore[misc]


# --- the evidence certificate --------------------------------------------


def test_unmeasured_scorer_is_not_a_passing_scorer():
    sc = ScorerEvidence(name="leave-one-mechanism-out", status=UNMEASURED,
                        detail="no cached run")
    assert not sc.beats_null()


def test_missing_lomo_cache_reads_as_unmeasured(tmp_path):
    cert = certify(lomo_cache=tmp_path / "absent.json")
    lomo = cert.scorers["leave-one-mechanism-out"]
    assert lomo.status == UNMEASURED
    assert not cert.unlocks_pass


def test_malformed_lomo_cache_is_unmeasured_not_a_crash(tmp_path):
    bad = tmp_path / "lomo_certificate.json"
    bad.write_text("{not json at all")
    cert = certify(lomo_cache=bad)
    assert cert.scorers["leave-one-mechanism-out"].status == UNMEASURED

    incomplete = tmp_path / "incomplete.json"
    incomplete.write_text(json.dumps({"mae": 1.0}))     # no "folds"
    cert = certify(lomo_cache=incomplete)
    assert cert.scorers["leave-one-mechanism-out"].status == UNMEASURED


def test_live_certificate_does_not_unlock_pass():
    """The measurement as it stands. If this ever fails, the science moved --
    re-read the number before changing the test."""
    cert = certify()
    loo = cert.scorers["leave-one-arm-out"]
    assert loo.status == MEASURED
    assert loo.mae > loo.null_mae, "LOO now beats its null; re-read gate/__init__.py"
    assert not cert.unlocks_pass
    assert cert.blocking_reasons()


def test_one_unit_cannot_bound_its_own_error():
    lo, hi = _paired_bootstrap([-5.0], resamples=100, confidence=0.95)
    assert lo == float("-inf") and hi == float("inf")


def test_bootstrap_ci_brackets_a_clear_win():
    lo, hi = _paired_bootstrap([-9.0, -10.0, -11.0, -10.5], resamples=2000,
                               confidence=0.95)
    assert hi < 0


def test_a_scorer_that_wins_on_average_but_not_reliably_is_refused():
    """Mean beats the null, resampling says it could be chance. Refused."""
    sc = ScorerEvidence(name="leave-one-arm-out", status=MEASURED, n_units=4,
                        mae=10.0, null_mae=15.0, ci_low=-12.0, ci_high=+2.0,
                        mae_subset=10.0, null_mae_subset=15.0)
    assert not sc.beats_null()


def test_headline_win_with_a_losing_placebo_subset_is_refused():
    """Clause (2): the discriminating arms must improve too."""
    sc = ScorerEvidence(name="leave-one-arm-out", status=MEASURED, n_units=12,
                        mae=5.0, null_mae=15.0, ci_low=-14.0, ci_high=-6.0,
                        mae_subset=14.9, null_mae_subset=15.0)
    assert not sc.beats_null()


def test_arm_holdout_alone_does_not_unlock_pass():
    """Clause (4): a new mechanism cannot inherit a score earned on old ones."""
    cert = EvidenceCertificate(scorers={
        "leave-one-arm-out": _passing_scorer("leave-one-arm-out"),
        "leave-one-mechanism-out": ScorerEvidence(
            name="leave-one-mechanism-out", status=UNMEASURED, detail="not run"),
    })
    assert not cert.unlocks_pass


# --- the device ----------------------------------------------------------


def test_a_doomed_candidate_is_killed_without_consulting_the_certificate(doomed):
    """A KILL must not need a certificate at all."""
    v = decide(doomed, certificate=None)
    assert v.kind is VerdictKind.KILL
    assert v.certificate is None, "a kill computed the certificate it did not need"


def test_survivors_abstain_on_the_live_certificate(survivor):
    v = decide(survivor)
    assert v.kind is VerdictKind.ABSTAIN
    assert any("PASS unavailable" in r for r in v.reasons)


def test_the_pass_path_is_live_code(survivor):
    """With an unlocking certificate the same survivor passes. If this test
    fails, PASS has become unreachable by construction and the device is a
    hard-coded refusal wearing a criterion."""
    v = decide(survivor, certificate=_passing_certificate())
    assert v.kind is VerdictKind.PASS
    assert v.is_pass


def test_an_unlocking_certificate_still_cannot_revive_a_killed_candidate(doomed):
    v = decide(doomed, certificate=_passing_certificate())
    assert v.kind is VerdictKind.KILL


def test_no_verdict_claims_validation(survivor):
    for cert in (None, _passing_certificate()):
        assert decide(survivor, certificate=cert).validated is False


def test_decide_all_preserves_input_order():
    """Sorting verdicts would be ranking by the back door."""
    profiles = [
        MechanismProfile(label="delta-", source="test", delta=0.5),
        MechanismProfile(label="gamma_E+", source="test", gamma_E=1.5),
        MechanismProfile(label="gamma_R-", source="test", gamma_R=0.5),
    ]
    verdicts = decide_all(profiles)
    assert [v.candidate for v in verdicts] == [p.label for p in profiles]


def test_the_device_exposes_no_score_or_rank(survivor):
    """The refusal `screen.rank_candidates` already makes, kept at this layer."""
    v = decide(survivor, certificate=_passing_certificate())
    for banned in ("score", "rank", "effect_size", "predicted_change"):
        assert not hasattr(v, banned), f"Verdict exposes {banned}"
    assert set(Verdict.__dataclass_fields__) == {
        "kind", "candidate", "model", "reasons", "screen_result", "certificate",
        "criterion"}


def test_verdict_explains_itself_without_a_certificate(doomed):
    text = decide(doomed, certificate=None).explain()
    assert "validated=False" in text


# --- model provenance ----------------------------------------------------


def test_every_verdict_names_its_model(survivor):
    """Survival is model-relative: two models the gate cannot tell apart
    disagree on 11 of 40 survivors, so an unlabelled verdict is unreproducible."""
    v = decide(survivor)
    assert v.model.model == PUBLISHED
    assert v.model.label() == "velez2011"
    assert "velez2011" in v.explain()


def test_model_is_required_and_has_no_default():
    """A default would let a verdict be built without naming its model, which is
    the whole failure this field exists to prevent."""
    assert Verdict.__dataclass_fields__["model"].default is dataclasses.MISSING


def test_the_label_carries_the_parameter_not_just_the_family():
    """K=50000 and K=2000 are both 'the extension' and behave oppositely."""
    ext = from_table({"carrying_capacity": 2000.0})
    assert ext.label() == "velez2011+K=2000"
    assert from_table({"carrying_capacity": 50000.0}).label() == "velez2011+K=50000"
    assert ext.is_extension


def test_a_table_without_a_capacity_key_is_the_published_model():
    assert from_table({}).model == PUBLISHED
    assert from_table({}).label() == "velez2011"


def test_an_unnamed_parameter_is_printed_rather_than_dropped():
    p = ModelProvenance(model="x", parameters={"some_new_dial": 3.0})
    assert "some_new_dial=3" in p.label()


def test_mixing_two_models_raises_rather_than_picking_one():
    """A verdict assembled from two models is about neither."""
    with pytest.raises(ModelMismatch):
        combined(from_simulator(), from_table({"carrying_capacity": 2000.0}))


def test_the_live_device_halves_agree_today():
    """If this fails, the cached table was rebuilt under a different model than
    the kill filters run, and every verdict was about to be incoherent."""
    from backtest.lomo import load

    assert combined(from_simulator(), from_table(load())).label() == "velez2011"
