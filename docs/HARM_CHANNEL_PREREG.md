# Harm channel — pre-registered enlargement of the target set

**Written and committed 2026-09-20, BEFORE any expression value for the new genes
was fetched.** That ordering is the only thing that makes the number at the end
mean anything, so it is a commit of its own.

## Why

`bricks/harm_channel.py` scores a drug target by how far its expression leans
toward regulatory T cells rather than the effector T cells a DMT means to
suppress. On the six targets in this repo's arm set with a usable denominator,
the two MS harm cases rank first and second, and the entire statistical claim is

    P(both harm cases in the top 2 of 6 | rank is arbitrary) = 1 / C(6,2) = 0.067

which does not reach 0.05, and which the channel was *designed* after seeing.
Enlarging the scored set is the cheapest available way to make that test mean
something — and it only works if the set and the labels are fixed before the
numbers are looked at.

## Inclusion rule

A gene enters the scored set if and only if:

1. it is the direct molecular target of a drug given in a **randomized
   controlled trial in multiple sclerosis**;
2. that trial's outcome is **published**, and the citation was verified against
   the Europe PMC REST API in this session — PMIDs below, none recalled;
3. the gene is not already scored, and the drug's primary target is not already
   scored (so ofatumumab adds nothing to MS4A1, ozanimod and siponimod add
   nothing to S1PR1).

**There is no expression filter in this rule.** Whether a Treg:effector ratio is
computable at all is decided afterwards, mechanically, by the existing
`MIN_EFFECTOR_NTPM = 1.0` floor. A gene the floor rejects is reported as
unclassifiable by this channel — it is not quietly dropped, and it does not get
re-included because its number looked good.

## Label rule

Assigned from the trial and post-marketing record, before scoring:

| label | definition |
|---|---|
| **REGULATORY-FAILURE HARM** | disease worsening, or new autoimmunity, attributable to the drug |
| **OTHER HARM** | serious toxicity of a different mechanism — opportunistic infection from over-suppression, hepatotoxicity, cardiac. **Counted as a negative for this channel, deliberately.** It claims to predict lost regulation; a channel that flags every dangerous drug claims nothing. Natalizumab/PML is the existing precedent (ranked fifth of six, correctly). |
| **NO HARM** | neither |

## The candidates, verified 2026-09-20

| gene | drug | trial | PMID | outcome as reported | label |
|---|---|---|---|---|---|
| IL12B | ustekinumab | phase II RRMS | 18703004 | no effect on lesion counts | NO HARM |
| CD80 | abatacept | ACCLAIM | 27481207 | stopped for futility | NO HARM |
| CD86 | abatacept | ACCLAIM | 27481207 | same trial; both genes scored, since CTLA4-Ig binds both | NO HARM |
| CD19 | inebilizumab | phase 1 RRMS | 29143550 | safe and tolerated | NO HARM |
| IL17A | secukinumab | randomized proof-of-concept RRMS | 27142710 | lesion reduction, PoC only | NO HARM |
| LINGO1 | opicinumab | AFFINITY part 1 | 41454463 | no benefit | NO HARM |
| BTK | evobrutinib | evolutionRMS1/2 | 39307151 | did not beat teriflunomide; liver-enzyme elevations | OTHER HARM |
| BTK | tolebrutinib | GEMINI 1/2, HERCULES | 40202623, 40202696 | hepatotoxicity signal | OTHER HARM |
| DHODH | teriflunomide | TEMSO | 21991951 | −31.5% ARR | NO HARM |
| HCAR2 | dimethyl fumarate | DEFINE | 22992073 | −53% ARR | NO HARM |

BTK appears once as a gene and twice as a drug; it is one target and is scored
once.

## Pre-declared exclusions

Listed so nobody can think the set was pruned after the ratios were seen.

- **vatelizumab (ITGA2)** — EMPIRE was terminated and no trial outcome paper was
  found in Europe PMC on 2026-09-20. Fails rule 2. (A mechanistic paper exists,
  PMID 30783682, and it reports that VLA-2 blockade *induces* Tregs — which is
  exactly the kind of thing that must not be allowed in through the back door
  once the labels are fixed.)
- **tabalumab / belimumab (TNFSF13B)** — no MS randomized trial. Fails rule 1.
- **IFN-γ administration, and anti-IFN-γ** — the harm case here is giving a
  cytokine, not removing a target. The channel's claim is about *removing*
  regulation, so this is out of scope rather than a missed hit.
- **IFN-β, glatiramer acetate, mitoxantrone, cladribine** — no single-gene
  molecular target to score.
- **rituximab, ublituximab, ofatumumab (MS4A1); ozanimod, siponimod (S1PR1)** —
  target already scored. Fails rule 3.

## The statistic, fixed before the numbers

The same exact combinatorial test, over however many targets clear the floor:

    P(both REGULATORY-FAILURE HARM targets in the top 2 of N) = 1 / C(N,2)

`N` is whatever the floor leaves. The test is not re-chosen afterwards, and no
alternative statistic is computed if this one disappoints.

- At N = 12, both harm cases top-2 gives p = 1/66 = **0.015**.
- **Falsifier, stated now:** any NO-HARM target scoring above TNFRSF1B (2.83x)
  breaks the top-2 and the claim weakens accordingly. IL2RA sits at 31.42x and
  is unlikely to be displaced; TNFRSF1B at 2.83x is the vulnerable one.

## What this still cannot buy

**The enlargement adds negatives only.** No new regulatory-failure harm case
with an independent immune target exists in the published MS record — the harm
set is two drugs, daclizumab and lenercept (atacicept is a third but its target
is B-cell and the floor excludes it). That is a permanent ceiling on this
channel's power, not a gap that more work closes.

And this is **not** a prospective test. The channel's design came after knowing
which drugs harmed; only the negatives here are new. A smaller p-value makes the
hypothesis more worth testing prospectively. It does not make it validated.
`validated=False` stands.

## Result

*(appended after scoring — nothing above this line changes)*
