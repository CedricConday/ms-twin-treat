"""B7 — the intervention brick: what therapy is being tested, in model terms.

This is the input end of the whole pipeline. Everything downstream asks it the
same question — "how much is the autoimmune attack suppressed?" — and gets one
number back.

Deliberately a plain dict, NOT wrapped by `spine.pipeline.standin()`. The QSP
and ABM bricks read `state["intervention"]["treat"]` directly, so a wrapper
would turn a treated run into a silently untreated one: the lookup would miss,
the default would apply, and the drug would do nothing with no error anywhere.
Unvalidated status travels as `validated: False` inside the dict instead, which
is the convention the other bricks already use.

**What is real here:** the parameter plumbing. An intervention is specified
once and every scale honours it.

**What is NOT real:** the mapping from a named drug to a `treat` value. There is
no dose-response calibration behind these numbers. `IFN_BETA.treat = 0.5` does
not mean interferon beta suppresses 50% of anything — it means "a moderate
effect" in a toy model. Calibrating that mapping against real trial outcomes is
exactly the job the backtest harness exists for, and it has not been done.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any

from bricks.grounding import IMMUNOGENIC, NEUTRAL, SUPPRESSIVE, mechanism_to_params


@dataclass(frozen=True)
class Intervention:
    """One therapy, expressed in the only terms the toy models understand.

    treat        0.0 = untreated, 1.0 = attack fully suppressed. Illustrative.
    immunogenic  0.0 = neutral; >0 = the therapy PROVOKES the autoreactive
                 response instead of calming it (e.g. an altered peptide ligand
                 that activates encephalitogenic T cells). This is the mechanism
                 by which a therapy can HARM. Set from a drug's known immunology,
                 never fit to its clinical outcome.
    cns_required does the agent have to cross into the CNS to work? Antigen-
                 specific tolerance acts peripherally (False); a remyelination
                 agent must reach the lesion (True). B6 uses this to decide
                 whether barrier penetration gates the effect.
    """

    name: str
    treat: float = 0.0
    immunogenic: float = 0.0
    dose: float = 1.0
    cns_required: bool = False
    mechanism: str = ""
    notes: str = ""
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not 0.0 <= self.treat <= 1.0:
            raise ValueError(f"treat must be in [0,1], got {self.treat}")
        if not 0.0 <= self.immunogenic <= 1.0:
            raise ValueError(f"immunogenic must be in [0,1], got {self.immunogenic}")
        if self.dose < 0:
            raise ValueError(f"dose must be >= 0, got {self.dose}")

    def as_state(self) -> dict[str, Any]:
        d = asdict(self)
        d["validated"] = False  # the plumbing is real; the numbers are not
        return d


# --------------------------------------------------------------------------- #
# The library of arms. Parameters are NOT hand-set per drug — they come from the
# mechanism-class rule in bricks/grounding.py, keyed on each drug's INDEPENDENT
# (in-vitro) mechanism, never on its clinical outcome. The trial outcomes in the
# notes are REAL and citable; they are the TARGETS the gate checks against, not
# inputs to the parameters. This is what lets the clinical gate test the rule
# rather than restate four hand-tuned numbers.
# --------------------------------------------------------------------------- #

def _from_mechanism(name: str, mechanism: str, *, cns_required: bool = False,
                    notes: str = "", strength: float | None = None,
                    disrupts_regulation: bool = False) -> Intervention:
    treat, immuno = mechanism_to_params(mechanism, strength, disrupts_regulation)
    return Intervention(name, treat=treat, immunogenic=immuno,
                        cns_required=cns_required, mechanism=mechanism, notes=notes)


UNTREATED = _from_mechanism("untreated", NEUTRAL, notes="control arm")

IFN_BETA = _from_mechanism(
    "IFN-beta", SUPPRESSIVE,
    notes="approved DMT; immunomodulatory (suppressive class). The Kang 2018 dataset "
          "in this repo is IFN-beta-stimulated PBMCs. Backtest target: WORKED.")

GLATIRAMER = _from_mechanism(
    "glatiramer acetate", SUPPRESSIVE,
    notes="approved 1996; tolerance-adjacent, immunomodulatory (suppressive class), "
          "acts peripherally. Backtest target: WORKED.")

APL_CGP77116 = _from_mechanism(
    "APL CGP77116", IMMUNOGENIC,
    notes="altered peptide ligand of MBP 83-99. IMMUNOGENIC class: encephalitogenic in "
          "T-cell assays — an INDEPENDENT, in-vitro property (Bielekova et al., Nat Med "
          "2000, doi:10.1038/80516). Phase II HALTED: 3 exacerbations, 2 drug-linked. "
          "Backtest target: HARMED. Its harm parameter comes from the immunogenic class "
          "rule, NOT from its relapse number.")

# --------------------------------------------------------------------------- #
# Expanded arm set (BUILD_PLAN §8 blocker 3). Every mechanism below is the drug's
# INDEPENDENT pharmacology -- its target, known before any trial read out -- and
# never its clinical result. The trial outcomes these arms are scored against
# live in backtest/clinical.py and docs/TRIAL_ANCHORS.md; nothing here reads them.
#
# cns_required is False throughout: each of these acts on peripheral immune cells
# or at the barrier itself, none needs to reach a lesion to work.
# --------------------------------------------------------------------------- #

NATALIZUMAB = _from_mechanism(
    "natalizumab", SUPPRESSIVE,
    notes="anti-alpha4-integrin monoclonal antibody; blocks leukocyte transit across the "
          "blood-brain barrier. Suppressive class from its target, not from AFFIRM.")

FINGOLIMOD = _from_mechanism(
    "fingolimod", SUPPRESSIVE,
    notes="sphingosine-1-phosphate receptor modulator; sequesters lymphocytes in lymph "
          "nodes by blocking egress. Suppressive class from its target.")

PONESIMOD = _from_mechanism(
    "ponesimod", SUPPRESSIVE,
    notes="selective S1P1 receptor modulator; same egress-blocking mechanism as "
          "fingolimod, which is why the class rule cannot tell the two apart.")

TERIFLUNOMIDE = _from_mechanism(
    "teriflunomide", SUPPRESSIVE,
    notes="dihydroorotate dehydrogenase inhibitor; cytostatic to proliferating "
          "lymphocytes via pyrimidine synthesis. Suppressive class from its target.")

DIMETHYL_FUMARATE = _from_mechanism(
    "dimethyl fumarate", SUPPRESSIVE,
    notes="Nrf2 pathway activator with immunomodulatory and anti-oxidative effects. "
          "Suppressive class from its mechanism.")

OCRELIZUMAB = _from_mechanism(
    "ocrelizumab", SUPPRESSIVE,
    notes="anti-CD20 monoclonal antibody; depletes B cells. Suppressive class from its "
          "target.")

ALEMTUZUMAB = _from_mechanism(
    "alemtuzumab", SUPPRESSIVE,
    notes="anti-CD52 monoclonal antibody; depletes circulating T and B lymphocytes. "
          "Suppressive class from its target.")

# The two counterexamples. Both are immunosuppressive by mechanism and both HARMED
# patients in controlled trials. The class rule predicts benefit for each, so both
# are arms the current rule is expected to get wrong -- that is why they are here.
LENERCEPT = _from_mechanism(
    "lenercept", SUPPRESSIVE, disrupts_regulation=True,
    notes="TNF receptor p55-IgG fusion protein; neutralizes TNF. Suppressive by "
          "mechanism, and ALSO regulation-disrupting: TNF-deficient mice develop severe "
          "MOG-induced EAE and TNF treatment reduces its severity (Liu, Nat Med 1998, "
          "PMID 9427610), and TNFR2 is required for oligodendrocyte progenitor "
          "proliferation and remyelination (Arnett, Nat Neurosci 2001, PMID 11600888). "
          "Both flags come from animal-model mechanism work, not from the trial. "
          "Backtest target: HARMED (PMID 10449104).")

ATACICEPT = _from_mechanism(
    "atacicept", SUPPRESSIVE,
    notes="TACI-Ig fusion protein; blocks BAFF and APRIL and depletes plasma cells. "
          "Suppressive by mechanism. Backtest target: HARMED (ATAMS halted for a raised "
          "relapse rate, PMID 24613349). NOT flagged regulation-disrupting: the nearest "
          "independent evidence is for anti-CD20 depletion removing IL-10-producing "
          "regulatory B cells before disease onset (Matsushita, J Clin Invest 2008, "
          "PMID 18802481), which is a different target and carries a timing dependence "
          "this model has no way to express. Flagging it on that basis would be reasoning "
          "backwards from the trial, so it is left unflagged and the gate keeps failing "
          "it. This is the honest gap, not an oversight.")

OFATUMUMAB = _from_mechanism(
    "ofatumumab", SUPPRESSIVE,
    notes="fully human anti-CD20 monoclonal antibody, subcutaneous; depletes B cells. "
          "Same target as ocrelizumab, which is why it shares its intervention point. "
          "Suppressive class from its target, not from ASCLEPIOS.")

DACLIZUMAB = _from_mechanism(
    "daclizumab", SUPPRESSIVE,
    notes="anti-CD25 (IL-2 receptor alpha) monoclonal antibody. Blocking CD25 cuts "
          "regulatory T-cell numbers by roughly 50% over 52 weeks while raising IL-2 "
          "bioavailability, which expands CD56bright NK cells. Withdrawn from the market "
          "in 2018 after 12 cases of severe encephalitis, 3-4 fatal. Suppressive class "
          "from its target; the Treg reduction is an INDEPENDENT pharmacodynamic "
          "observation, not a reading of DECIDE.")

CLADRIBINE = _from_mechanism(
    "cladribine", SUPPRESSIVE,
    notes="purine nucleoside analogue; selectively cytotoxic to lymphocytes, producing "
          "sustained depletion followed by reconstitution. Suppressive class from its "
          "mechanism, not from CLARITY.")

IFN_GAMMA = _from_mechanism(
    "IFN-gamma", IMMUNOGENIC,
    notes="recombinant type II interferon; pro-inflammatory, activates macrophages and "
          "upregulates MHC class II. IMMUNOGENIC class from its immunology (Panitch "
          "1987, PMID 2882294). Backtest target: HARMED.")

# Added 2026-09-21 alongside bricks/profiles.py. Class comes from the target,
# never from the trial -- ustekinumab and abatacept are SUPPRESSIVE by mechanism
# and both failed their trials, which is the point of including them.
UBLITUXIMAB = _from_mechanism(
    "ublituximab", SUPPRESSIVE,
    notes="glycoengineered anti-CD20 monoclonal antibody; depletes B cells. Third "
          "arm on the same target as ocrelizumab and ofatumumab, sharing their "
          "intervention point and target-level magnitude.")

IFN_BETA_1B = _from_mechanism(
    "IFN-beta-1b", SUPPRESSIVE,
    notes="type I interferon, same molecular class as IFN-beta-1a. Suppressive from "
          "the class, not from the 1993 trial.")

OZANIMOD = _from_mechanism(
    "ozanimod", SUPPRESSIVE,
    notes="S1P1/S1P5 receptor modulator; blocks lymphocyte egress from lymph nodes. "
          "Lumped onto gamma_E with fingolimod and ponesimod.")

RITUXIMAB = _from_mechanism(
    "rituximab", SUPPRESSIVE,
    notes="chimeric anti-CD20 monoclonal antibody; depletes B cells. Enters the gate "
          "direction-only: RIFUND-MS reports the proportion of patients relapsing, "
          "not an annualised rate.")

USTEKINUMAB = _from_mechanism(
    "ustekinumab", SUPPRESSIVE,
    notes="anti-IL-12/23 p40; blocks the cytokine signals driving Th1 and Th17 "
          "effector differentiation. Suppressive by mechanism and its phase II found "
          "NO significant lesion reduction at any dose -- the class rule predicts "
          "benefit and the trial says otherwise, which is why the arm is here.")

ABATACEPT = _from_mechanism(
    "abatacept", SUPPRESSIVE,
    notes="CTLA4-Ig; blocks the CD28 costimulation T-cell activation requires. "
          "Suppressive by mechanism; ACCLAIM closed early for futility. No "
          "regulation-disrupting flag: CTLA4-Ig plausibly impairs Tregs but no MS "
          "source was found for it, and flagging it to explain the outcome would be "
          "reasoning backwards from the trial.")

LIBRARY = {i.name: i for i in (
    UNTREATED, IFN_BETA, GLATIRAMER, APL_CGP77116,
    NATALIZUMAB, FINGOLIMOD, PONESIMOD, TERIFLUNOMIDE, DIMETHYL_FUMARATE,
    OCRELIZUMAB, ALEMTUZUMAB, LENERCEPT, ATACICEPT, IFN_GAMMA,
    OFATUMUMAB, DACLIZUMAB, CLADRIBINE,
    UBLITUXIMAB, IFN_BETA_1B, OZANIMOD, RITUXIMAB, USTEKINUMAB, ABATACEPT,
)}


class InterventionStage:
    """Stage: write the chosen intervention into the state.

    Per-patient override: if the incoming state carries `intervention_name`,
    that arm is selected. Lets a virtual cohort be split across arms without
    rebuilding the pipeline.
    """

    name = "B7 intervention"
    requires: tuple[str, ...] = ()

    def __init__(self, intervention: Intervention | str = IFN_BETA) -> None:
        if isinstance(intervention, str):
            if intervention not in LIBRARY:
                raise KeyError(f"unknown arm {intervention!r}; have {sorted(LIBRARY)}")
            intervention = LIBRARY[intervention]
        self.intervention = intervention

    def run(self, state: dict) -> dict:
        chosen = self.intervention
        override = state.get("intervention_name")
        if override:
            if override not in LIBRARY:
                raise KeyError(f"unknown arm {override!r}; have {sorted(LIBRARY)}")
            chosen = LIBRARY[override]
        state["intervention"] = chosen.as_state()
        return state

    __call__ = run
