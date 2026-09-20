"""Tests for the per-arm intervention-point map.

The point of this file is to make three things un-quiet:

1. An arm's intervention points are a MECHANISM CLAIM and must stay sourced.
2. Magnitudes are a single shared stub, not potencies, and must not drift into
   per-drug numbers without a source.
3. The model's ceiling — which arms it cannot tell apart — is asserted, not
   noticed later. If a future model change separates them, these tests fail and
   say so, which is the correct way to find out.
"""

from __future__ import annotations

import pytest

from bricks.profiles import (
    FITTED,
    LUMPED,
    PROFILES,
    STUB_MAGNITUDE,
    degeneracy_report,
    touched_points,
    with_magnitudes,
)
from bricks.qsp_velez import INTERVENTION_POINTS


def test_every_arm_carries_a_source():
    for name, prof in PROFILES.items():
        assert prof.source.strip(), f"{name} has no source — mechanism claims must be cited"


def test_every_multiplier_is_the_shared_stub_or_a_cited_fit():
    """No per-drug magnitude may appear without a citation.

    A hand-picked number per drug is indistinguishable from fitting, which is
    the failure the leave-one-arm-out test exists to catch. The only escape is
    an entry in FITTED, which carries the source that supplied the number.
    """
    allowed = {1.0, 1.0 - STUB_MAGNITUDE, 1.0 + STUB_MAGNITUDE}
    for name, prof in PROFILES.items():
        for point in INTERVENTION_POINTS:
            value = getattr(prof, point)
            if value in allowed:
                continue
            assert (name, point) in FITTED, (
                f"{name}.{point}={value} is neither the shared stub nor listed in "
                "FITTED. Per-drug magnitudes need a source."
            )


def test_every_fitted_magnitude_is_cited_on_its_arm():
    """FITTED is a registry, not a bypass: the citation must be on the profile too."""
    for (name, point), note in FITTED.items():
        assert getattr(PROFILES[name], point) != 1.0, (
            f"{name}.{point} is registered as fitted but does not move the dial")
        assert note.strip()
        assert "Martinez-Pasamar" in PROFILES[name].source or "PMC" in PROFILES[name].source


def test_ocrelizumab_escaped_the_gamma_E_lump_via_a_measured_parameter():
    """The one arm pulled out of the ceiling, and the only way that happens.

    Martinez-Pasamar 2013 reproduced post-anti-CD20 T-cell dynamics by moving
    K_eff from 1000 to ~850 cells against EAE data, "independently of the
    alpha_reg parameter". So ocrelizumab acts on ke, not gamma_E, and its
    alpha_R stays at 1.0 because the same paper says B-cell depletion works
    "without strengthening T_reg activation".
    """
    ocr = PROFILES["ocrelizumab"]
    assert touched_points(ocr) == ("ke",)
    assert ocr.ke == 0.85
    assert ocr.gamma_E == 1.0
    assert ocr.alpha_R == 1.0
    assert "ocrelizumab" not in LUMPED


def test_atacicept_was_not_routed_through_ke_by_analogy():
    """Both are B-lineage agents; only one has the measurement.

    Extending the anti-CD20 result to a BAFF/APRIL blocker would be exactly the
    reasoning-by-analogy that turns a sourced map into a fitted one.
    """
    ata = PROFILES["atacicept"]
    assert ata.ke == 1.0
    assert touched_points(ata) == ("gamma_E",)
    assert "atacicept" in LUMPED


def test_untreated_touches_nothing():
    assert touched_points(PROFILES["untreated"]) == ()


def test_the_arm_set_matches_the_intervention_library():
    """profiles.py and intervention.py must describe the same 14 arms."""
    from bricks.intervention import LIBRARY
    assert set(PROFILES) == set(LIBRARY)


# --------------------------------------------------------------------------- #
# mechanism claims that carry the most weight — pinned so they cannot drift
# --------------------------------------------------------------------------- #

def test_glatiramer_acts_at_antigen_presentation():
    """The only arm whose primary mechanism sits at the APC step.

    MHC class II competition (Arnon & Aharoni, PNAS 2004) -> delta down, plus
    restoration of the deficient FoxP3+ regulatory population -> alpha_R up.
    """
    ga = PROFILES["glatiramer acetate"]
    assert ga.delta < 1.0
    assert ga.alpha_R > 1.0
    assert ga.alpha_E == 1.0


def test_lenercept_is_the_one_arm_that_suppresses_and_deregulates():
    """What the two-constant scheme could not express.

    alpha_E down (TNF neutralisation is immunosuppressive) AND alpha_R down
    (TNF/TNFR2 is required for regulation and remyelination — Liu 1998,
    Arnett 2001). Whether that nets to harm is left to the dynamics.
    """
    len_ = PROFILES["lenercept"]
    assert len_.alpha_E < 1.0
    assert len_.alpha_R < 1.0
    assert len_ is not PROFILES["IFN-beta"]
    assert touched_points(len_) != touched_points(PROFILES["IFN-beta"])


def test_alemtuzumab_is_not_given_a_treg_depleting_axis():
    """The sign trap. CD52 is on Tregs, so 'it depletes Tregs too' is tempting.

    It is also wrong on the timescale a trial measures: repopulation after
    alemtuzumab is Treg-BIASED, with regulatory cells returning faster than
    effectors and with increased suppressive capacity (PMC4519957, PMC8581537).
    A naive alpha_R < 1 would have the sign backwards.
    """
    alem = PROFILES["alemtuzumab"]
    assert alem.alpha_R == 1.0
    assert alem.gamma_R == 1.0


def test_teriflunomide_treg_axis_is_left_unset_because_evidence_is_contested():
    """Klotz 2019 reports no Treg change; PMID 40879143 reports impaired function.

    The repo does not get to pick a side. Pinned so that a later 'improvement'
    setting this axis has to argue with the citation in bricks/profiles.py.
    """
    teri = PROFILES["teriflunomide"]
    assert teri.alpha_R == 1.0
    assert teri.alpha_E < 1.0


def test_harm_arms_do_not_share_a_single_harm_constant():
    """IFN-gamma and APL act at antigen presentation and the naive pool;
    lenercept harms by de-regulation. Different mechanisms, different points —
    which is the whole reason for this representation."""
    assert touched_points(PROFILES["IFN-gamma"]) != touched_points(PROFILES["lenercept"])


# --------------------------------------------------------------------------- #
# the ceiling — asserted, not discovered later
# --------------------------------------------------------------------------- #

def test_the_gamma_E_cluster_is_a_model_ceiling_not_a_stub_artifact():
    """Five arms still collapse onto one dial, and potency will not fix it.

    natalizumab, fingolimod, ponesimod, alemtuzumab and atacicept all reduce the
    active effector pool by mechanisms the model has no structure for — no CNS
    compartment, no lymph node, no B cells. Ocrelizumab used to be the sixth;
    it left because someone measured its effect on T-cell dynamics.
    """
    assert set(LUMPED) == {
        "natalizumab", "fingolimod", "ponesimod",
        "alemtuzumab", "atacicept",
    }
    for name in LUMPED:
        assert touched_points(PROFILES[name]) == ("gamma_E",)


def test_a_worked_drug_and_a_harmful_drug_are_the_same_object():
    """The sharpest statement of the ceiling, and the reason it blocks the screen.

    Natalizumab reduced relapses 68% (AFFIRM). Atacicept RAISED them and ATAMS
    was halted. In this model they are identical. No magnitude, no fit and no
    amount of leave-one-out can make one right without making the other wrong —
    only a model with the missing structure can.
    """
    assert PROFILES["natalizumab"].as_dict() == {
        **PROFILES["atacicept"].as_dict(),
        "label": "natalizumab",
        "source": PROFILES["natalizumab"].source,
    }


def test_degeneracy_report_separates_model_ceiling_from_stub_artifact():
    rep = degeneracy_report()
    assert rep["model"], "the gamma_E cluster must be reported as a model ceiling"
    assert rep["stub"], "same-axis arms must be reported as stub-only degeneracy"
    assert ("natalizumab", "atacicept") in rep["model"]
    assert ("teriflunomide", "dimethyl fumarate") in rep["stub"]
    assert not (set(rep["model"]) & set(rep["stub"]))


def test_representation_separates_more_than_the_old_two_classes():
    patterns = {touched_points(p) for p in PROFILES.values()}
    assert len(patterns) > 2


# --------------------------------------------------------------------------- #
# the migration path off the stub
# --------------------------------------------------------------------------- #

def test_with_magnitudes_replaces_a_stub_on_an_axis_the_arm_already_acts_on():
    fitted = with_magnitudes("IFN-beta", alpha_E=0.83)
    assert fitted.alpha_E == 0.83
    assert "fitted" in fitted.source


def test_with_magnitudes_refuses_to_invent_a_new_intervention_point():
    """Turning on a new axis is a mechanism claim and belongs in profiles.py."""
    with pytest.raises(ValueError, match="does not act on"):
        with_magnitudes("IFN-beta", gamma_E=1.4)


def test_with_magnitudes_refuses_to_silently_delete_a_mechanism():
    with pytest.raises(ValueError, match="silently remove"):
        with_magnitudes("lenercept", alpha_R=1.0)


def test_with_magnitudes_rejects_a_non_intervention_point():
    with pytest.raises(ValueError, match="not an intervention point"):
        with_magnitudes("IFN-beta", kr=2.0)
