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
                    notes: str = "", strength: float | None = None) -> Intervention:
    treat, immuno = mechanism_to_params(mechanism, strength)
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

LIBRARY = {i.name: i for i in (UNTREATED, IFN_BETA, GLATIRAMER, APL_CGP77116)}


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
