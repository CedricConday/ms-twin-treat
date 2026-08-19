"""Mechanism -> model-parameter rule — grounding the intervention brick.

The clinical gate is only a real test if the parameters it grades were NOT set
from the answers. Hand-setting `treat`/`immunogenic` per drug fails that: each
number could be quietly nudged until its arm came out right.

This module replaces the per-drug numbers with a GENERAL rule keyed on each
drug's **mechanism class** — an independent, in-vitro-measurable property, known
before any trial reads out:

  - immunosuppressive / tolerizing : the drug calms the autoreactive response
                                     -> treat > 0, immunogenic = 0
  - immunogenic                    : the drug PROVOKES it (e.g. an altered peptide
                                     ligand that is encephalitogenic in T-cell
                                     assays) -> treat = 0, immunogenic > 0
  - none                           : untreated control -> (0, 0)

The class strengths below are single, class-level constants set from mechanism
reasoning — NOT fit per drug. So the whole library is parameterised by just two
numbers, and the clinical gate becomes a test of THIS RULE + the models, rather
than a restatement of four hand-tuned values.

Honest bound: with only four arms this is a *coarse* rule, and a rigorous
leave-one-arm-out (fit the strengths on N-1 arms, predict the Nth) needs more
arms than we have data for. What changed is real but modest: per-arm fitting is
gone; a two-parameter mechanism rule now stands in its place.
"""

from __future__ import annotations

SUPPRESSIVE = "immunosuppressive/tolerizing"
IMMUNOGENIC = "immunogenic"
NEUTRAL = "none"

# Class-level strengths. Set from mechanism class, identical within a class,
# never tuned to a drug's own clinical outcome.
SUPPRESSIVE_STRENGTH = 0.5
IMMUNOGENIC_STRENGTH = 0.4


def mechanism_to_params(mechanism: str, strength: float | None = None) -> tuple[float, float]:
    """Map an INDEPENDENT mechanism class to (treat, immunogenic).

    No clinical outcome enters here. `strength` lets a caller override the class
    default from independent potency data; left None, the class constant is used
    (which keeps every drug in a class identical — the honest, un-fit default).
    """
    if mechanism == SUPPRESSIVE:
        return (SUPPRESSIVE_STRENGTH if strength is None else float(strength), 0.0)
    if mechanism == IMMUNOGENIC:
        return (0.0, IMMUNOGENIC_STRENGTH if strength is None else float(strength))
    if mechanism == NEUTRAL:
        return (0.0, 0.0)
    raise ValueError(f"unknown mechanism class: {mechanism!r}")
