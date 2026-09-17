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

The class strengths below are single, class-level constants — NOT fit per drug.
So the whole library is parameterised by just two numbers, and the clinical gate
becomes a test of THIS RULE + the models, rather than a restatement of four
hand-tuned values.

SUPPRESSIVE_STRENGTH is DATA-GROUNDED (2026-09-06): it is the IFN-beta effect
magnitude measured in Kang 2018 (GSE96583), expressed as a fraction of the
distance between two different immune cell identities in the same data — the
in-data yardstick for the wholesale change of immune cell state that treat=1.0
is supposed to mean. Arithmetic in `scripts/derive_suppressive_strength.py`,
guarded by `tests/test_grounding.py`. Kang is an in-vitro PBMC stimulation
experiment: no relapse rate, no trial arm, no clinical endpoint enters it, so
the clinical gate remains an out-of-sample test of this number.

IMMUNOGENIC_STRENGTH is still set from mechanism reasoning. Nothing in Kang
speaks to it — Kang has no encephalitogenic arm — and only one arm in the
library (APL CGP77116) exercises it. Inventing a derivation for it would be
worse than leaving it visibly ungrounded, so it is left alone and labelled.

Honest bounds:
  - one grounded number, one reasoned number. The rule is still coarse.
  - the derived quantity is the SIZE of IFN-beta's effect on immune cells; the
    model reads that size as "fraction of the autoreactive attack suppressed".
    That equation is an assumption, not a measurement.
  - with only four arms a rigorous leave-one-arm-out (fit the strengths on N-1
    arms, predict the Nth) still needs more arms than we have data for.
"""

from __future__ import annotations

SUPPRESSIVE = "immunosuppressive/tolerizing"
IMMUNOGENIC = "immunogenic"
NEUTRAL = "none"

# --- the derivation inputs, recorded so a drift is a test failure, not a story --
# Kang 2018 (GSE96583) via data/kang.py, log-normalized, 15706 genes, the 7 cell
# types that clear the harness reliability bar (Megakaryocytes, 63/69 cells,
# reliability 0.045, dropped from both terms).
KANG_IFNB_DELTA_NORM = 15.902     # mean over cell types of ||mean_IFNb - mean_ctrl||
KANG_IDENTITY_DISTANCE = 20.325   # mean over cell-type PAIRS of ||ctrl_i - ctrl_j||

# Class-level strengths. Identical within a class, never tuned to a drug's own
# clinical outcome.
#
# GROUNDED — Kang 2018 IFN-beta magnitude: 15.902 / 20.325 = 0.7824 -> 0.78.
# IFN-beta moves an immune cell ~78% of the way to a different cell identity, and
# the model reads that as the fraction of the attack a suppressive drug removes.
# Cited the way barrier.py cites Pardridge 2019; redo it with
# `PYTHONPATH=. python scripts/derive_suppressive_strength.py`.
SUPPRESSIVE_STRENGTH = 0.78

# NOT grounded — mechanism reasoning only. Kang has no encephalitogenic arm, and
# no independent source for it has been pinned. Left where it was, on purpose.
IMMUNOGENIC_STRENGTH = 0.4

# --- the third thing a mechanism can be (2026-09-17) -------------------------
# REGULATORY_DISRUPTION. Two classes were not enough: a drug can suppress the
# immune system AND make the disease worse, because its target also carries a
# regulatory or reparative function. Lenercept (TNF blockade) and atacicept
# (BAFF/APRIL blockade) are both immunosuppressive and both harmed patients, and
# no value of SUPPRESSIVE_STRENGTH can express that — it is the rule's shape.
#
# The flag is per drug and set from INDEPENDENT biology, never from the trial:
#
#   TNF blockade — TNF-deficient mice develop SEVERE MOG-induced EAE with high
#   mortality, and treating them with TNF reduces severity (Liu et al., Nat Med
#   1998;4(1):78-83, PMID 9427610); TNFR2 is required for oligodendrocyte
#   progenitor proliferation and remyelination (Arnett et al., Nat Neurosci
#   2001;4(11):1116-22, PMID 11600888). Removing TNF removes a brake and a repair
#   signal. Both are animal-model mechanism papers, not readings of the lenercept
#   trial.
#
# A flagged drug keeps its suppressive strength and additionally gets
# IMMUNOGENIC_STRENGTH on the harm channel. No new constant is introduced: one
# tuned until lenercept came out harmful would be fitted to the outcome the gate
# is supposed to test. Whether the existing two numbers produce net harm is then
# a PREDICTION of this rule, and it is allowed to be wrong.


def mechanism_to_params(mechanism: str, strength: float | None = None,
                       disrupts_regulation: bool = False) -> tuple[float, float]:
    """Map an INDEPENDENT mechanism class to (treat, immunogenic).

    No clinical outcome enters here. `strength` lets a caller override the class
    default from independent potency data; left None, the class constant is used
    (which keeps every drug in a class identical — the honest, un-fit default).

    `disrupts_regulation` marks a suppressive drug whose target ALSO carries a
    regulatory or reparative function that the drug removes (see
    REGULATORY_DISRUPTION above). Such a drug gets both channels: it suppresses
    the attack AND provokes it. It reuses IMMUNOGENIC_STRENGTH rather than
    introducing a third constant, deliberately — a new number chosen to make the
    known-harmful arms come out harmful would be fitted to the answer, which is
    the one thing this module exists not to do.
    """
    if mechanism == SUPPRESSIVE:
        treat = SUPPRESSIVE_STRENGTH if strength is None else float(strength)
        return (treat, IMMUNOGENIC_STRENGTH if disrupts_regulation else 0.0)
    if mechanism == IMMUNOGENIC:
        return (0.0, IMMUNOGENIC_STRENGTH if strength is None else float(strength))
    if mechanism == NEUTRAL:
        return (0.0, 0.0)
    raise ValueError(f"unknown mechanism class: {mechanism!r}")
