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


@dataclass(frozen=True)
class Intervention:
    """One therapy, expressed in the only terms the toy models understand.

    treat        0.0 = untreated, 1.0 = attack fully suppressed. Illustrative.
    cns_required does the agent have to cross into the CNS to work? Antigen-
                 specific tolerance acts peripherally (False); a remyelination
                 agent must reach the lesion (True). B6 uses this to decide
                 whether barrier penetration gates the effect.
    """

    name: str
    treat: float = 0.0
    dose: float = 1.0
    cns_required: bool = False
    notes: str = ""
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not 0.0 <= self.treat <= 1.0:
            raise ValueError(f"treat must be in [0,1], got {self.treat}")
        if self.dose < 0:
            raise ValueError(f"dose must be >= 0, got {self.dose}")

    def as_state(self) -> dict[str, Any]:
        d = asdict(self)
        d["validated"] = False  # the plumbing is real; the numbers are not
        return d


# --------------------------------------------------------------------------- #
# A small library of arms, so a backtest can be pointed at a known outcome.
#
# The trial outcomes in the comments are REAL and citable. The `treat` values
# are NOT derived from them — they are placeholders awaiting calibration. The
# point of naming real arms is that this pipeline has somewhere to aim: a model
# earns trust by separating the ones that worked from the ones that harmed.
# --------------------------------------------------------------------------- #

UNTREATED = Intervention(
    "untreated", treat=0.0,
    notes="control arm",
)

IFN_BETA = Intervention(
    "IFN-beta", treat=0.5, cns_required=False,
    notes="approved DMT; the Kang 2018 dataset in this repo is IFN-beta stimulated PBMCs",
)

GLATIRAMER = Intervention(
    "glatiramer acetate", treat=0.45, cns_required=False,
    notes="approved 1996; tolerance-adjacent, acts peripherally. Backtest target: WORKED",
)

APL_CGP77116 = Intervention(
    "APL CGP77116", treat=-0.0, cns_required=False,
    notes=(
        "altered peptide ligand of MBP 83-99. Phase II HALTED: 3 patients had MS "
        "exacerbations, 2 immunologically linked to the drug (Bielekova et al., "
        "Nature Medicine 2000). Backtest target: HARMED. treat=0 is a placeholder - "
        "a model that earns its keep would predict a NEGATIVE effect here, and this "
        "brick cannot yet express that."
    ),
)

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
