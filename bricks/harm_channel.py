"""A harm channel that is not MRI and not the trial's own relapse rate.

THE PROBLEM THIS EXISTS FOR
---------------------------
`backtest/potency.py` fits each drug's magnitude on its trial's MRI lesion
ratio, which avoids the circularity of fitting on the outcome the gate scores.
It has one documented hole that more data does not close: **lenercept reported no
significant MRI difference while its relapse rate rose** (Neurology 1999;53:457,
PMID 10449104). Through the MRI channel its potency is ~0 — "does nothing" — for
a drug that harmed people.

WHAT WAS RULED OUT FIRST
------------------------
The obvious candidate was EAE, and it is empirically dead as a ranking signal:

    Berg I et al. "Which experimental factors govern successful animal-to-human
    translation in multiple sclerosis drug development?" eBioMedicine
    2024;110:105434. PMID 39515028. 497 studies, 15 approved and 11 failed DMTs
    (failed on efficacy OR safety), ~30,000 animals:

    "There was no association between animal study outcomes ... and successful
     approval."

plus 91% of those studies were published after the first-in-MS trial. So EAE
*efficacy outcomes* do not rank candidates. (Mechanism-level EAE findings — what
a target does — are a different claim and still stand; see bricks/profiles.py.)

THE CHANNEL
-----------
A drug that removes or blocks a target expressed on REGULATORY cells removes
regulation along with whatever it was aimed at. So score a target by where it is
expressed, with no trial outcome read at all:

    treg_ratio = nTPM(T-reg) / mean nTPM(effector T subsets)

Data: Human Protein Atlas immune cell RNA (`rna_immune_cell.tsv`,
https://www.proteinatlas.org/download/tsv/rna_immune_cell.tsv.zip), consensus
nTPM per sorted immune cell type. HPA is CC BY-SA 3.0. The values below are
transcribed for the eight targets in this repo's arm set;
`scripts/derive_harm_channel.py` re-downloads and re-derives them.

THE EXPRESSION FLOOR, AND WHY THE BEST-LOOKING RESULT IS THROWN OUT
--------------------------------------------------------------------
Raw, the ranking puts all three harm cases in the top three, which looks
decisive. It is partly an artifact and the artifact is removed here.

**TNFRSF13B (atacicept's target, TACI) ranks first at 88x only because its
denominator is ~0.025 nTPM.** TACI is a B-cell gene — 52 to 131 nTPM in B cells,
2.2 in Tregs, essentially zero in effector T cells. A ratio of two numbers that
are both near zero is not a signal about regulation, it is division by noise.

So `MIN_EFFECTOR_NTPM` requires a real denominator. Atacicept then becomes
**unclassifiable by this channel** rather than its top hit, which is the honest
outcome: TACI acts on B cells and a Treg:effector T ratio cannot see it. That
loses the prettiest number in the analysis and keeps the one that means
something.

WHAT SURVIVES, AND HOW STRONG IT ACTUALLY IS
----------------------------------------------
Of the six targets an arm maps to with a usable denominator, the two harm cases
rank **first and second**:

    IL2RA     daclizumab   31.42x    withdrawn 2018, fatal encephalitis
    TNFRSF1B  lenercept     2.83x    harm; relapses rose
    CD52      alemtuzumab   1.81x    approved, worked
    S1PR1     fingolimod    0.74x    approved, worked
    ITGA4     natalizumab   0.34x    approved, worked
    MS4A1     anti-CD20     0.29x    approved, worked

Under a null where rank is arbitrary, the chance that both harm cases land in the
top two of six is 1/C(6,2) = **1/15, about 0.067**. That is the whole statistical
claim and it is deliberately unimpressive: six targets, two events, one exact
combinatorial test, and NO correction for the fact that this channel was designed
after seeing which drugs harmed. It does not reach 0.05 and it is not meant to.
**It is a hypothesis with a suggestive first look, not a validated predictor.**
Treat the number as a reason to test it prospectively, not as evidence.

(TNFRSF1A, TNF's other receptor, sits at 0.77x and is not in the arm map —
lenercept is scored on TNFRSF1B because TNFR2 is the receptor Tregs depend on.
That choice is a mechanism claim and it favours the hypothesis, so it is named
here rather than buried.)

AN OUT-OF-SAMPLE CHECK THE CHANNEL PASSES
------------------------------------------
The 0.067 above is a binary test — harm case or not — on the two drugs whose
MS trials showed harm. There is a stronger check available that was NOT used to
build the channel, and it uses the four drugs the binary test calls negatives.

The channel predicts **regulatory failure**, not harm in general. So look at what
each drug's actual toxicity IS:

    CD52    alemtuzumab  1.81x   SECONDARY AUTOIMMUNITY in 30-48% of patients.
                                 Thyroid autoimmunity alone reached 42% over six
                                 years in the pooled CARE-MS studies, plus ITP
                                 (~2.2%) and autoimmune nephropathy (~0.34%).
    S1PR1   fingolimod   0.74x   no secondary-autoimmunity signal
    ITGA4   natalizumab  0.34x   PML — an opportunistic INFECTION from
                                 immunosuppression, not autoimmunity
    MS4A1   anti-CD20    0.29x   no secondary-autoimmunity signal

Alemtuzumab is the only one of the four with a major autoimmunity signature, and
it is the only one of the four the channel puts above 1.0. Across all six
targets the ordering tracks regulatory-failure toxicity monotonically: fatal
encephalitis, relapse worsening, 30-48% secondary autoimmunity, then three drugs
with none of it.

Natalizumab is the case that sharpens the claim rather than weakening it. It is
a dangerous drug — PML kills — and the channel ranks it fifth of six. That is
correct behaviour: PML is infection from over-suppression, which is a different
failure mode, and a channel keyed on Treg expression should not and does not
flag it.

Also at target level, not counted in the statistic because it is the SAME target
as lenercept: TNF blockade as a class causes CNS demyelination. 122 published
cases between 1990 and 2016 across etanercept (34%), infliximab (29%),
adalimumab (29%) and certolizumab (8%). Four more drugs, one more confirmation
of TNFRSF1B, zero additional independent targets.

WHY IT IS WORTH HAVING ANYWAY
------------------------------
It sees two things nothing else in this repo can:

  * **Lenercept**, which the MRI channel scores as inert.
  * **Daclizumab**, which every outcome-based channel scores as a SUCCESS — it
    cut relapses 45% (DECIDE) — and which was withdrawn for fatal encephalitis.
    Neither ARR nor MRI can express "worked and then killed people". A channel
    keyed on the target rather than the outcome can.

And the separation is not driven by the drugs working better or worse: CD52 sits
at 1.81 and alemtuzumab is one of the most effective DMTs there is. The claim is
about a regulatory-cell liability, not about efficacy.

`validated=False`. This is RNA in healthy blood, not protein, not function, not
MS tissue.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb

# Effector T-cell subsets used as the denominator. B cells, NK and monocytes are
# deliberately excluded: the question is whether a target is preferentially on
# REGULATORY T cells relative to the effector T cells a DMT means to suppress.
EFFECTOR_SUBSETS = (
    "memory CD4 T-cell", "naive CD4 T-cell",
    "memory CD8 T-cell", "naive CD8 T-cell",
)

# A denominator below this is noise, and a ratio built on it is meaningless.
# 1.0 nTPM is the conventional "not expressed" floor for bulk RNA.
MIN_EFFECTOR_NTPM = 1.0

# Human Protein Atlas consensus immune-cell nTPM, transcribed 2026-09-20 from
# rna_immune_cell.tsv. Re-derive with scripts/derive_harm_channel.py.
HPA_NTPM: dict[str, dict[str, float]] = {
    "IL2RA": {"T-reg": 271.0, "memory CD4 T-cell": 29.2, "naive CD4 T-cell": 2.0,
              "memory CD8 T-cell": 3.2, "naive CD8 T-cell": 0.1,
              "naive B-cell": 8.1, "memory B-cell": 15.4, "NK-cell": 14.3},
    "TNFRSF1A": {"T-reg": 22.1, "memory CD4 T-cell": 21.3, "naive CD4 T-cell": 19.5,
                 "memory CD8 T-cell": 45.1, "naive CD8 T-cell": 29.6,
                 "naive B-cell": 0.0, "memory B-cell": 0.1, "NK-cell": 7.9},
    "TNFRSF1B": {"T-reg": 70.9, "memory CD4 T-cell": 34.4, "naive CD4 T-cell": 9.8,
                 "memory CD8 T-cell": 40.4, "naive CD8 T-cell": 15.5,
                 "naive B-cell": 8.4, "memory B-cell": 13.5, "NK-cell": 13.4},
    "MS4A1": {"T-reg": 0.8, "memory CD4 T-cell": 2.8, "naive CD4 T-cell": 0.5,
              "memory CD8 T-cell": 5.5, "naive CD8 T-cell": 2.1,
              "naive B-cell": 630.0, "memory B-cell": 627.9, "NK-cell": 6.5},
    "CD52": {"T-reg": 11722.3, "memory CD4 T-cell": 8205.6, "naive CD4 T-cell": 5514.4,
             "memory CD8 T-cell": 6533.6, "naive CD8 T-cell": 5683.6,
             "naive B-cell": 10119.7, "memory B-cell": 10618.5, "NK-cell": 3743.6},
    "ITGA4": {"T-reg": 13.4, "memory CD4 T-cell": 47.0, "naive CD4 T-cell": 26.6,
              "memory CD8 T-cell": 52.4, "naive CD8 T-cell": 33.2,
              "naive B-cell": 33.1, "memory B-cell": 37.4, "NK-cell": 32.3},
    "S1PR1": {"T-reg": 330.3, "memory CD4 T-cell": 504.1, "naive CD4 T-cell": 384.4,
              "memory CD8 T-cell": 492.1, "naive CD8 T-cell": 397.5,
              "naive B-cell": 348.2, "memory B-cell": 295.7, "NK-cell": 398.5},
    "TNFRSF13B": {"T-reg": 2.2, "memory CD4 T-cell": 0.0, "naive CD4 T-cell": 0.0,
                  "memory CD8 T-cell": 0.0, "naive CD8 T-cell": 0.1,
                  "naive B-cell": 52.1, "memory B-cell": 131.4, "NK-cell": 0.0},
}

# Arm -> the target this channel scores it on. An arm whose drug hits several
# targets is scored on the one whose loss would remove regulation.
ARM_TARGET: dict[str, str] = {
    "daclizumab": "IL2RA",
    "lenercept": "TNFRSF1B",       # TNFR2, the receptor Tregs depend on
    "alemtuzumab": "CD52",
    "natalizumab": "ITGA4",
    "fingolimod": "S1PR1",
    "ponesimod": "S1PR1",
    "ocrelizumab": "MS4A1",
    "ofatumumab": "MS4A1",
    "atacicept": "TNFRSF13B",      # excluded by the floor; see the docstring
}


@dataclass(frozen=True)
class HarmScore:
    """A target's regulatory-cell liability, with what is needed to distrust it."""

    arm: str
    target: str
    treg_ntpm: float
    effector_ntpm: float
    treg_ratio: float | None   # None when the denominator is below the floor
    classifiable: bool
    reason: str = ""

    @property
    def validated(self) -> bool:
        return False


def score(arm: str) -> HarmScore:
    """Regulatory-cell liability for one arm. Reads no trial outcome."""
    target = ARM_TARGET.get(arm)
    if target is None:
        raise KeyError(f"no target mapped for {arm!r}")
    v = HPA_NTPM[target]
    treg = v["T-reg"]
    eff = [v[c] for c in EFFECTOR_SUBSETS if c in v]
    eff_mean = sum(eff) / len(eff) if eff else 0.0

    if eff_mean < MIN_EFFECTOR_NTPM:
        return HarmScore(
            arm=arm, target=target, treg_ntpm=treg, effector_ntpm=eff_mean,
            treg_ratio=None, classifiable=False,
            reason=(f"effector-T expression {eff_mean:.3f} nTPM is below the "
                    f"{MIN_EFFECTOR_NTPM} floor — a ratio here divides by noise. "
                    "This target is not an effector-T gene, so a Treg:effector-T "
                    "ratio cannot speak to it."),
        )
    return HarmScore(arm=arm, target=target, treg_ntpm=treg,
                     effector_ntpm=eff_mean, treg_ratio=treg / eff_mean,
                     classifiable=True)


def ranking() -> list[HarmScore]:
    """Classifiable arms, most Treg-biased target first. One row per target."""
    seen, out = set(), []
    for arm in ARM_TARGET:
        s = score(arm)
        if s.classifiable and s.target not in seen:
            seen.add(s.target)
            out.append(s)
    return sorted(out, key=lambda s: s.treg_ratio, reverse=True)


def top_k_probability(n_items: int, n_harm: int) -> float:
    """Chance that all `n_harm` harm cases land in the top `n_harm` of `n_items`.

    An exact combinatorial statement under a null of arbitrary ranking, and the
    whole statistical claim this module makes. It does NOT correct for the
    channel having been designed after seeing which drugs harmed.
    """
    if not 0 < n_harm <= n_items:
        raise ValueError("need 0 < n_harm <= n_items")
    return 1.0 / comb(n_items, n_harm)


if __name__ == "__main__":
    print("HARM CHANNEL — target expression on regulatory cells")
    print("  Human Protein Atlas immune-cell nTPM. No trial outcome is read.\n")

    # Outcomes shown for READING ONLY; nothing above consumes them.
    outcome = {"daclizumab": "withdrawn, fatal encephalitis",
               "lenercept": "HARMED, relapses rose",
               "alemtuzumab": "approved", "natalizumab": "approved",
               "fingolimod": "approved", "ocrelizumab": "approved"}

    rows = ranking()
    print(f"{'target':<11} {'arm':<13} {'T-reg':>9} {'effT':>8} {'ratio':>7}  outcome")
    print("-" * 78)
    for s in rows:
        print(f"{s.target:<11} {s.arm:<13} {s.treg_ntpm:>9.1f} {s.effector_ntpm:>8.1f} "
              f"{s.treg_ratio:>7.2f}  {outcome.get(s.arm, '')}")

    excluded = [score(a) for a in ARM_TARGET if not score(a).classifiable]
    for s in {e.target: e for e in excluded}.values():
        print(f"\n  EXCLUDED {s.target} ({s.arm}): {s.reason}")

    p = top_k_probability(len(rows), 2)
    print(f"\n  Both harm cases rank 1st and 2nd of {len(rows)}. Under arbitrary "
          f"ranking that is 1/{len(rows) * (len(rows) - 1) // 2} = {p:.3f}.")
    print("  n=7, two events, and the channel was designed after seeing which drugs")
    print("  harmed. A reason to test it prospectively, not evidence. validated=False.")
    print("\n  What it sees that nothing else here does: lenercept, which the MRI")
    print("  channel scores as inert, and daclizumab, which every outcome-based")
    print("  channel scores as a success because it cut relapses 45%.")
