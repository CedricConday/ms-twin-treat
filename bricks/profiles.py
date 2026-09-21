"""Blocker (4), first half — each arm's INTERVENTION POINTS, from pharmacology.

`bricks/grounding.py` gives a drug one of two class constants. That is what the
leave-one-arm-out test graded and failed: "a single class strength carries no
drug-specific information" (BUILD_PLAN §8.4, 2026-09-17). The fix is not a
better constant, it is a representation in which two drugs in the same class
can differ at all.

`bricks/qsp_velez.py` supplies that: a drug is a `MechanismProfile`, a set of
multipliers on the published model's own named rates. This module assigns those
points, per arm, **from each drug's independent pharmacology and never from its
trial outcome** — the same rule `grounding.py` already follows.

WHAT IS REAL HERE AND WHAT IS STUBBED
--------------------------------------
REAL: **which** parameters each drug touches, and in **which direction**. Every
assignment below is sourced.

STUBBED: **how much** — with one exception. Every non-unit multiplier is built
from one shared constant, `STUB_MAGNITUDE`, because inventing a number per drug
is exactly the failure this repo keeps testing for. Magnitudes arrive from
arm-level MRI lesion outcomes (BUILD_PLAN §8.4, "the MRI channel").

The exception is `ocrelizumab.ke = 0.85`, which is FITTED — Martinez-Pasamar et
al. 2013 reproduced post-anti-CD20 T-cell dynamics by moving K_eff from 1000 to
~850 cells against EAE flow-cytometry data. Fitted magnitudes live in `FITTED`
with their citation, and the test suite allows a non-stub value only for pairs
listed there.

That split is deliberate and it is diagnosable: `degeneracy_report()` separates
arm pairs that are identical because the MODEL cannot tell them apart from ones
identical only because the STUB gives every axis the same size. The first kind
is a ceiling. The second dissolves when potency lands.

WHAT THE MODEL CANNOT EXPRESS — read before trusting any ranking
----------------------------------------------------------------
Vélez de Mendizábal has no CNS compartment and no trafficking (verified against
the full text; see BUILD_PLAN §8.4). Loss from the active effector pool is one
lumped term, `gamma_E`. So **depletion, sequestration and transit blockade all
collapse onto the same dial**, and the model cannot distinguish:

    natalizumab   anti-alpha4-integrin, blocks BBB transit
    fingolimod    S1P modulator, blocks lymph-node egress
    ponesimod     selective S1P1 modulator, same egress block
    alemtuzumab   anti-CD52, depletes T and B lymphocytes
    atacicept     TACI-Ig, blocks BAFF/APRIL and depletes plasma cells

Each of those removes effectors from the pool that does damage, so `gamma_E` is
a defensible home for all five — but it is a LUMP, it is marked as one on every
affected arm, and it caps what a screen over this model can ever say.

**Ocrelizumab escaped that lump**, and how it escaped is the template for the
rest. Not by a better guess about depletion, but because the successor paper
measured what removing B cells does to T-cell dynamics and expressed it as a
parameter this model already has (`ke`). The lump breaks where someone has done
that experiment; it does not break by reasoning harder about the mechanism.
Atacicept is the counter-example held deliberately: also B-lineage, no such
measurement, so it stays lumped.
"""

from __future__ import annotations

from dataclasses import replace

from bricks.qsp_velez import INTERVENTION_POINTS, MechanismProfile

# --------------------------------------------------------------------------- #
# The single stub magnitude. NOT a potency.
# --------------------------------------------------------------------------- #
# One number for every axis of every drug, on purpose. A per-drug number chosen
# by hand would be indistinguishable from fitting, and a per-class number is
# what the LOO already rejected. Replaced axis by axis as MRI-derived potency
# lands; nothing should read this constant as an effect size.
STUB_MAGNITUDE = 0.5

SUPPRESS = 1.0 - STUB_MAGNITUDE   # 0.5 — dial turned down
ENHANCE = 1.0 + STUB_MAGNITUDE    # 1.5 — dial turned up

# Magnitudes that are NOT the stub because a source supplies the number. Keyed by
# (arm, point) so the test that forbids hand-picked magnitudes can tell a cited
# value from an invented one. Adding an entry here requires the citation to be in
# the arm's `source` string.
FITTED: dict[tuple[str, str], str] = {
    ("ofatumumab", "ke"): (
        "Martinez-Pasamar et al. 2013 (PMC3651362), applied at TARGET level: the same "
        "anti-CD20 K_eff shift as ocrelizumab, 850/1000 = 0.85. Not a per-drug fit."
    ),
    ("ublituximab", "ke"): (
        "Martinez-Pasamar et al. 2013 (PMC3651362), applied at TARGET level, exactly "
        "as for ofatumumab: the same anti-CD20 K_eff shift, 850/1000 = 0.85. Not a "
        "per-drug fit, and the point of the arm is that it CANNOT differ from the "
        "other anti-CD20 arms in this model."
    ),
    ("rituximab", "ke"): (
        "Martinez-Pasamar et al. 2013 (PMC3651362), applied at TARGET level as for "
        "the other anti-CD20 arms: 850/1000 = 0.85. Not a per-drug fit."
    ),
    ("ocrelizumab", "ke"): (
        "Martinez-Pasamar et al. 2013, BMC Syst Biol 7:34 (PMC3651362): anti-CD20 "
        "dynamics reproduced by reducing K_eff from 1000 to ~850 cells. 850/1000 = 0.85."
    ),
}

# Arms whose real mechanism the model cannot represent, collapsed onto gamma_E.
# Keyed by arm name so a report can name them without re-deriving the reason.
LUMPED: dict[str, str] = {
    "natalizumab": "BBB transit block — no CNS compartment in the model",
    "fingolimod": "lymph-node egress block — no lymph-node compartment",
    "ponesimod": "lymph-node egress block — no lymph-node compartment",
    "alemtuzumab": "pan-lymphocyte depletion — no cell-type resolution",
    "atacicept": "plasma-cell depletion — the model has no B cells",
    "cladribine": "lymphocyte-depleting purine analogue — no cell-type resolution",
}


def _p(name: str, source: str, **points: float) -> MechanismProfile:
    unknown = set(points) - set(INTERVENTION_POINTS)
    if unknown:
        raise ValueError(f"{name}: not intervention points: {sorted(unknown)}")
    return MechanismProfile(label=name, source=source, **points)


# --------------------------------------------------------------------------- #
# The arms. Sources are mechanism papers and pharmacology, never trial results.
# --------------------------------------------------------------------------- #

UNTREATED = _p("untreated", "control arm")

IFN_BETA = _p(
    "IFN-beta",
    "Type I interferon; reduces T-cell activation and proliferation and shifts "
    "the cytokine balance away from Th1. Mapped to alpha_E only. A Treg axis is "
    "NOT asserted: IFN-beta effects on Treg frequency are reported but this repo "
    "has not pinned an independent source, and guessing the sign of a second axis "
    "is how a mechanism map becomes a fit.",
    alpha_E=SUPPRESS,
)

GLATIRAMER = _p(
    "glatiramer acetate",
    "Binds promiscuously and with high affinity to class II MHC and competes with "
    "myelin antigens for presentation, displacing them from the binding groove "
    "(Arnon & Aharoni, PNAS 2004, doi:10.1073/pnas.0404887101) — that is an effect "
    "at ANTIGEN PRESENTATION, so delta. It also restores the CD4+CD25+FoxP3+ "
    "regulatory population that is deficient in MS, so alpha_R. The only arm in "
    "the set whose primary described mechanism sits at the APC step.",
    delta=SUPPRESS, alpha_R=ENHANCE,
)

NATALIZUMAB = _p(
    "natalizumab",
    "Anti-alpha4-integrin; blocks leukocyte transit across the blood-brain "
    "barrier. LUMPED onto gamma_E — the model has no CNS compartment, so 'cannot "
    "enter the tissue' and 'left the active pool' are the same event here.",
    gamma_E=ENHANCE,
)

FINGOLIMOD = _p(
    "fingolimod",
    "S1P receptor modulator; sequesters lymphocytes in lymph nodes by blocking "
    "egress. LUMPED onto gamma_E — no lymph-node compartment.",
    gamma_E=ENHANCE,
)

PONESIMOD = _p(
    "ponesimod",
    "Selective S1P1 modulator; the same egress-blocking mechanism as fingolimod. "
    "LUMPED onto gamma_E. This pair is genuinely indistinguishable in this model, "
    "which is a real limit and not a stub artifact.",
    gamma_E=ENHANCE,
)

TERIFLUNOMIDE = _p(
    "teriflunomide",
    "DHODH inhibitor; blocks de novo pyrimidine synthesis, which proliferating "
    "lymphocytes require and resting cells do not (they use the salvage pathway). "
    "That is a proliferation block, so alpha_E. The Treg axis is left at 1.0 "
    "because the evidence is CONTESTED and this repo does not get to pick: "
    "Klotz et al. (Sci Transl Med 2019, doi:10.1126/scitranslmed.aao5563) report "
    "no change in CD4 regulatory T-cell frequency or function, while later work "
    "(PMID 40879143) reports impaired FOXP3+ Treg function via mitochondrial "
    "respiration. Asserting either sign here would be a coin flip dressed as "
    "biology. Revisit when one side is settled.",
    alpha_E=SUPPRESS,
)

DIMETHYL_FUMARATE = _p(
    "dimethyl fumarate",
    "Nrf2 pathway activator with immunomodulatory and anti-oxidative effects; "
    "reduces effector lymphocyte expansion. Mapped to alpha_E.",
    alpha_E=SUPPRESS,
)

OCRELIZUMAB = _p(
    "ocrelizumab",
    "Anti-CD20 monoclonal antibody; depletes B cells. NOT lumped onto gamma_E — "
    "this arm has a mechanism and a magnitude from the successor paper. "
    "Martinez-Pasamar et al. 2013 (BMC Syst Biol 7:34, PMC3651362) ran a "
    "sensitivity analysis of these same equations against EAE flow-cytometry data "
    "and found the post-anti-CD20 T-cell dynamics were reproduced by reducing the "
    "K_eff threshold below the healthy standard, 1000 -> ~850 cells, "
    "'independently of the alpha_reg parameter'; their reading is that B-cell "
    "depletion 'prevents uncontrolled activation of T_eff without strengthening "
    "T_reg activation'. So: ke = 850/1000 = 0.85, alpha_R left at 1.0, and the "
    "magnitude is FITTED to mouse T-cell dynamics — never to a human relapse rate. "
    "The model still has no B cells; what it has is a sourced surrogate for what "
    "removing them does to T-cell activation.",
    ke=0.85,
)

ALEMTUZUMAB = _p(
    "alemtuzumab",
    "Anti-CD52; depletes T and B lymphocytes by CDC and ADCC. LUMPED onto "
    "gamma_E. Crucially it is NOT given a Treg-depleting axis, even though CD52 "
    "sits on Tregs too: repopulation after alemtuzumab is Treg-BIASED — "
    "regulatory T cells return faster than effectors, with an increased "
    "proportion and increased suppressive capacity (Cossburn et al., PMC4519957; "
    "PMC8581537). Over a trial horizon the net regulatory effect is upward, not "
    "downward, so a naive 'it depletes Tregs as well' would have the sign wrong. "
    "alpha_R is nonetheless left at 1.0: the direction is documented but the "
    "timescale separation is not something this model represents.",
    gamma_E=ENHANCE,
)

LENERCEPT = _p(
    "lenercept",
    "TNF receptor p55-IgG fusion protein; neutralizes TNF. Suppressive, so "
    "alpha_E. ALSO regulation-stripping, so alpha_R down: TNF-deficient mice "
    "develop severe MOG-induced EAE and TNF treatment reduces severity (Liu, Nat "
    "Med 1998, PMID 9427610), and TNFR2 is required for oligodendrocyte "
    "progenitor proliferation and remyelination (Arnett, Nat Neurosci 2001, PMID "
    "11600888). Both are animal-model mechanism papers, not readings of the "
    "trial. This is the arm the two-constant scheme could not express: whether "
    "suppression or de-regulation wins is now a prediction of the dynamics.",
    alpha_E=SUPPRESS, alpha_R=SUPPRESS,
)

ATACICEPT = _p(
    "atacicept",
    "TACI-Ig; blocks BAFF and APRIL and depletes plasma cells. LUMPED onto "
    "gamma_E — no B cells in the model. **Deliberately NOT routed through `ke` "
    "the way ocrelizumab is**, even though both are B-lineage agents: the "
    "Martinez-Pasamar result is specific to anti-CD20 B-cell depletion, and "
    "extending it to a BAFF/APRIL blocker would be reasoning by analogy, which "
    "is how a sourced map turns into a fitted one. Deliberately NOT given a "
    "regulation-stripping axis: the nearest independent evidence is anti-CD20 "
    "depletion removing IL-10-producing regulatory B cells before disease onset "
    "(Matsushita, J Clin Invest 2008, PMID 18802481), a different target with a "
    "timing dependence this model cannot express. Flagging it would be reasoning "
    "backwards from the ATAMS result. The gate keeps failing this arm; that is "
    "the honest gap, not an oversight.",
    gamma_E=ENHANCE,
)

OFATUMUMAB = _p(
    "ofatumumab",
    "Fully human anti-CD20 monoclonal antibody; depletes B cells. SAME TARGET as "
    "ocrelizumab, so it inherits the same intervention point and the same fitted "
    "magnitude: Martinez-Pasamar et al. 2013 (PMC3651362) reproduced post-anti-CD20 "
    "T-cell dynamics by moving K_eff from 1000 to ~850 cells. The EAE experiment "
    "behind that number used a specific anti-CD20 antibody, so applying it to a "
    "second anti-CD20 agent is a TARGET-LEVEL extension, not a per-drug "
    "measurement -- which is exactly why this arm is valuable: it and ocrelizumab "
    "now share a dial and a magnitude, so any difference the trials show between "
    "them is a difference this model provably cannot produce.",
    ke=0.85,
)

DACLIZUMAB = _p(
    "daclizumab",
    "Anti-CD25 (IL-2 receptor alpha). CD25 blockade cuts regulatory T-cell numbers "
    "by roughly 50% over 52 weeks (Bielekova and colleagues; see also PMID 25416807 "
    "on maintained Treg function under CD25 blockade), so alpha_R DOWN. That is a "
    "pharmacodynamic observation about the drug's target, independent of DECIDE. "
    "NOT REPRESENTABLE, and it matters here: the same blockade raises IL-2 "
    "bioavailability and expands CD56bright regulatory NK cells, which this model "
    "has no cell type for. So the model sees only the Treg loss and none of the "
    "compensation -- it should predict harm, and the trial reported a 45% relapse "
    "reduction. This arm is in the set precisely because it is a case the model is "
    "expected to get wrong for a stateable reason.",
    alpha_R=SUPPRESS,
)

CLADRIBINE = _p(
    "cladribine",
    "Purine nucleoside analogue; selectively cytotoxic to lymphocytes, giving "
    "sustained depletion then reconstitution. LUMPED onto gamma_E with the rest of "
    "the depleting class -- and see the qsp_velez docstring: this model cannot "
    "produce benefit from removing effector cells at all.",
    gamma_E=ENHANCE,
)

IFN_GAMMA = _p(
    "IFN-gamma",
    "Recombinant type II interferon; pro-inflammatory, activates macrophages and "
    "upregulates MHC class II. Upregulated antigen presentation is delta UP; the "
    "expanded pool of activatable autoreactive cells is naive_E UP. Both from its "
    "immunology (Panitch 1987, PMID 2882294), not from the exacerbations.",
    delta=ENHANCE, naive_E=ENHANCE,
)

APL_CGP77116 = _p(
    "APL CGP77116",
    "Altered peptide ligand of MBP 83-99; encephalitogenic in T-cell assays — an "
    "independent in-vitro property (Bielekova et al., Nat Med 2000, "
    "doi:10.1038/80516). It acts at the TCR/antigen step, so delta UP, and it "
    "expands the autoreactive effector pool, so naive_E UP.",
    delta=ENHANCE, naive_E=ENHANCE,
)

# --------------------------------------------------------------------------- #
# Added 2026-09-21. Anchors researched by the parallel session
# (docs/TRIAL_ANCHORS.md, f2ca0d4); the dial assignments below are this file's
# own call, from pharmacology, never from the trial number. Two arms the anchor
# table offers are deliberately NOT here — see UNASSIGNED at the bottom.

UBLITUXIMAB = _p(
    "ublituximab",
    "Glycoengineered anti-CD20 monoclonal antibody; depletes B cells. THIRD arm on "
    "the same target as ocrelizumab and ofatumumab, and it inherits the same "
    "target-level magnitude for the same reason they share one: Martinez-Pasamar "
    "et al. 2013 (PMC3651362) moved K_eff 1000 -> ~850 for anti-CD20. It is added "
    "because `ke` was the most starved multi-member group in the out-of-sample "
    "ceiling (scripts/dial_ceiling.py), not because a third anti-CD20 teaches the "
    "model anything new about B cells.",
    ke=0.85,
)

IFN_BETA_1B = _p(
    "IFN-beta-1b",
    "Type I interferon, same molecular class as IFN-beta-1a; reduces T-cell "
    "activation and proliferation. Mapped to alpha_E only, for the same reason and "
    "with the same omission as IFN-beta: no Treg axis is asserted without an "
    "independent source.",
    alpha_E=SUPPRESS,
)

OZANIMOD = _p(
    "ozanimod",
    "S1P1/S1P5 receptor modulator; sequesters lymphocytes in lymph nodes by "
    "blocking egress. LUMPED onto gamma_E with fingolimod and ponesimod — no "
    "lymph-node compartment. Added with low expectations: gamma_E is already the "
    "largest group and the one no potency reaches (backtest/potency.py returns "
    "OUT OF RANGE for every arm on it).",
    gamma_E=ENHANCE,
)

RITUXIMAB = _p(
    "rituximab",
    "Chimeric anti-CD20 monoclonal antibody; depletes B cells. Same target and "
    "same target-level magnitude as the other anti-CD20 arms: Martinez-Pasamar "
    "et al. 2013 (PMC3651362), K_eff 1000 -> ~850. Its RIFUND-MS "
    "outcome is a PROPORTION of patients relapsing, not an annualised rate, so it "
    "enters the gate as direction-only — converting a proportion to an ARR would "
    "invent precision, the same rule that keeps active-comparator arms unchained.",
    ke=0.85,
)

USTEKINUMAB = _p(
    "ustekinumab",
    "Anti-IL-12/23 p40 monoclonal antibody; blocks the cytokine signals driving "
    "Th1 and Th17 effector differentiation, so it is mapped to alpha_E. **Its "
    "phase II found no significant lesion reduction at any of four doses** "
    "(PMID 18703004), which is exactly why it is valuable: alpha_E previously held "
    "three arms that all worked, and a dial whose training data is all successes "
    "cannot teach a model that the dial sometimes does nothing.",
    alpha_E=SUPPRESS,
)

ABATACEPT = _p(
    "abatacept",
    "CTLA4-Ig; binds CD80/CD86 on antigen-presenting cells and blocks the CD28 "
    "costimulation that T-cell activation requires. Mapped to `delta`, the "
    "activation step, NOT to alpha_E — it does not slow proliferation of already "
    "activated cells, it prevents activation. **No alpha_R axis is asserted.** "
    "CTLA4-Ig plausibly impairs regulatory T cells, which depend on CD28, and if "
    "it did this arm would join glatiramer's group and turn a singleton into a "
    "scored pair — which is precisely the incentive to be careful. A search on "
    "2026-09-21 found abatacept/Treg evidence only in LRBA deficiency, "
    "transplantation and haemolytic anaemia, nothing in MS, so the axis stays off.",
    delta=SUPPRESS,
)

TOVAXIN = _p(
    "Tovaxin",
    "Autologous attenuated myelin-reactive T cells (T-cell vaccination). The "
    "induced response deletes myelin-reactive EFFECTOR clones, so the defensible "
    "dial is gamma_E — effector loss — not a regulatory axis. **That assignment "
    "is deliberately the unhelpful one.** Pairing it with daclizumab on alpha_R "
    "would have turned a singleton into a scored group and improved the "
    "out-of-sample ceiling in the same commit that reports it; gamma_E is "
    "already the largest group and buys that measurement nothing. Direction-only: "
    "TERMS (n=150, RRMS, randomised placebo-controlled) MISSED its primary "
    "endpoint, and its ARR figure is post-hoc against a placebo arm the paper "
    "itself notes was lowered by prior DMT, so no magnitude is claimed. It is "
    "here because gamma_E's six arms all WORKED, and a dial whose training data "
    "is all successes cannot teach a model that the dial sometimes does nothing.",
    gamma_E=ENHANCE,
)

# NOT ASSIGNED, and the reasons matter more than the arms.
#
# laquinimod (ALLEGRO, PMID 22417253, -23.1%): aryl hydrocarbon receptor agonism
#   shifting myeloid cells toward an anti-inflammatory phenotype. No named Vélez
#   rate corresponds. It is the weakest quantified effect available and therefore
#   exactly the spread the arm set needs, which is the whole reason forcing it
#   onto a dial would be fitting the representation to the outcome.
# secukinumab (PoC, PMID 27142710, primary endpoint missed): anti-IL-17A. IL-17
#   is not a species in this model and effector function is not one of its rates.

PROFILES: dict[str, MechanismProfile] = {p.label: p for p in (
    UNTREATED, IFN_BETA, GLATIRAMER, APL_CGP77116,
    NATALIZUMAB, FINGOLIMOD, PONESIMOD, TERIFLUNOMIDE, DIMETHYL_FUMARATE,
    OCRELIZUMAB, ALEMTUZUMAB, LENERCEPT, ATACICEPT, IFN_GAMMA,
    OFATUMUMAB, DACLIZUMAB, CLADRIBINE,
    UBLITUXIMAB, IFN_BETA_1B, OZANIMOD, RITUXIMAB, USTEKINUMAB, ABATACEPT,
    TOVAXIN,
)}


def with_magnitudes(name: str, **points: float) -> MechanismProfile:
    """Return an arm's profile with some multipliers replaced by real potency.

    The migration path off `STUB_MAGNITUDE`, one axis at a time, as MRI-derived
    numbers land. Only axes the arm already touches may be changed: turning on a
    new intervention point is a change to the arm's *mechanism*, which belongs in
    this file with a citation, not in a caller.
    """
    base = PROFILES[name]
    for point, value in points.items():
        if point not in INTERVENTION_POINTS:
            raise ValueError(f"{point!r} is not an intervention point")
        if getattr(base, point) == 1.0:
            raise ValueError(
                f"{name} does not act on {point!r}. Adding an intervention point is a "
                "mechanism claim — edit bricks/profiles.py with a source instead."
            )
        if value == 1.0:
            raise ValueError(
                f"setting {name}.{point} to 1.0 would silently remove a sourced "
                "mechanism; remove it in this file if that is the intent"
            )
    return replace(base, **points, source=base.source + " [magnitude: fitted]")


def touched_points(profile: MechanismProfile) -> tuple[str, ...]:
    """The intervention points this arm actually moves."""
    return tuple(p for p in INTERVENTION_POINTS if getattr(profile, p) != 1.0)


def degeneracy_report() -> dict[str, list[tuple[str, str]]]:
    """Which arms this representation still cannot tell apart, and why.

    Two very different failures look identical in a results table, so they are
    separated here:

      "model"  — the arms touch the same points because the MODEL lacks the
                 structure to distinguish their mechanisms. A ceiling. Potency
                 numbers will not fix it; a different model would.
      "stub"   — the arms touch DIFFERENT points, or the same points for
                 different sourced reasons, and are only identical because
                 STUB_MAGNITUDE gives every axis the same size. Dissolves as
                 soon as real magnitudes land.

    An empty "model" list would mean the vocabulary separates every arm. It does
    not, and the entries say exactly what would have to change.
    """
    out: dict[str, list[tuple[str, str]]] = {"model": [], "stub": []}
    names = [n for n in PROFILES if n != "untreated"]
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            pa, pb = PROFILES[a], PROFILES[b]
            if touched_points(pa) != touched_points(pb):
                continue
            if all(getattr(pa, p) == getattr(pb, p) for p in INTERVENTION_POINTS):
                kind = "model" if (a in LUMPED and b in LUMPED) else "stub"
                out[kind].append((a, b))
    return out


if __name__ == "__main__":
    print("Blocker (4), first half — intervention points per arm.")
    print("Magnitudes are ONE shared stub; only the pattern is claimed.\n")

    width = max(len(n) for n in PROFILES)
    for name, prof in PROFILES.items():
        pts = touched_points(prof)
        desc = ", ".join(f"{p}={getattr(prof, p):g}" for p in pts) or "-"
        lump = "  [LUMPED]" if name in LUMPED else ""
        print(f"  {name:<{width}}  {desc}{lump}")

    rep = degeneracy_report()
    print(f"\nIndistinguishable pairs — MODEL ceiling ({len(rep['model'])}):")
    for a, b in rep["model"]:
        print(f"  {a} ~ {b}")
        print(f"      {LUMPED[a]}")
    print(f"\nIndistinguishable pairs — STUB only ({len(rep['stub'])}), "
          "these separate when potency lands:")
    for a, b in rep["stub"]:
        print(f"  {a} ~ {b}")

    old_classes = 2
    print(f"\n  old scheme: {old_classes} distinguishable mechanism classes")
    print(f"  now:        {len({touched_points(p) for p in PROFILES.values()})} "
          "distinct intervention patterns")
    print("  (patterns are sourced; magnitudes are not yet — see the docstring)")
