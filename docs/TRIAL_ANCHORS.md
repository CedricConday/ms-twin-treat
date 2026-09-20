# Trial anchors — the exam the clinical gate sits

Every arm the clinical gate scores against needs three things, kept separate on
purpose:

1. **a clinical outcome** — the number the trial reported, with its comparator;
2. **a mechanism class** — assigned from the drug's *independent, in-vitro*
   pharmacology, never from the outcome in column 1;
3. **a citation for each**, so the assignment can be checked rather than trusted.

Mixing (1) into (2) is what makes a gate circular, and it is the failure mode
`GROUNDING.md` exists to prevent. The mechanism column below was written from
each drug's target, not from whether it worked.

## Comparators are not interchangeable

Three of these trials are **active-comparator** studies: their headline number is
against interferon beta-1a or teriflunomide, not placebo. Converting those to a
"vs placebo" figure by chaining two trials' effects would invent precision that
does not exist. They are stored with their real comparator and scored the way the
trial ran them — simulated drug arm against the simulated *comparator* arm.

## The anchors

Numbers verified against PubMed on 2026-09-17. Where the trial reported both arms'
annualized relapse rates (ARR), the change is computed from those; otherwise the
paper's own reported reduction is used.

Three arms carry a direction but **no scored magnitude**: lenercept and IFN-gamma
because neither paper reports a relapse-rate ratio, and atacicept because its three
dose arms are non-monotonic (+126% / +108% / +158%), the confidence intervals are
wide, and the trial was halted early. A single number picked out of that would be a
choice dressed as a measurement.

| arm | trial | comparator | ARR (drug vs comp.) | change | direction | PMID |
|---|---|---|---|---|---|---|
| untreated | — | — | — | 0% | neutral | control arm |
| IFN-beta | PRISMS 1998 | placebo | 1.73–1.82 vs 2.56 relapses | −30% | improves | 9820297 |
| glatiramer acetate | CONFIRM 2012 | placebo | 0.29 vs 0.40 | −29% | improves | 22992072 |
| natalizumab | AFFIRM 2006 | placebo | rate at 1 yr, −68% reported | −68% | improves | 16510744 |
| fingolimod 0.5 mg | FREEDOMS 2010 | placebo | 0.18 vs 0.40 | −55% | improves | 20089952 |
| teriflunomide 14 mg | TEMSO 2011 | placebo | 0.37 vs 0.54 | −31.5% | improves | 21991951 |
| dimethyl fumarate BID | DEFINE 2012 | placebo | 0.17 vs 0.36 | −53% | improves | 22992073 |
| ocrelizumab | OPERA I/II 2017 | IFN beta-1a | 0.16 vs 0.29 | −46% | improves | 28002679 |
| alemtuzumab | CARE-MS I 2012 | IFN beta-1a | rate ratio 0.45 | −55% | improves | 23122652 |
| ponesimod | OPTIMUM 2021 | teriflunomide | 0.202 vs 0.290 | −30.5% | improves | 33779698 |
| lenercept | Lenercept MS Study Group 1999 | placebo | more patients with exacerbations, occurring earlier (p=0.006) | unquantified | **harms** | 10449104 |
| atacicept | ATAMS 2014 | placebo | 0.86 / 0.79 / 0.98 (25 / 75 / 150 mg) vs 0.38 | **raised** (+126% / +108% / +158%) | **harms** | 24613349 |
| IFN-gamma | Panitch 1987 | pre/post within-patient | 7 of 18 patients had exacerbations on treatment | unquantified | **harms** | 2882294 |
| APL CGP77116 | Bielekova 2000 | (phase II, halted) | 3 exacerbations, 2 drug-linked | unquantified | **harms** | 11017150 |

## Mechanism classes, assigned from pharmacology

| arm | target / mechanism | class |
|---|---|---|
| IFN-beta | type I interferon; immunomodulatory | suppressive |
| glatiramer acetate | MBP-like random copolymer; peripheral immunomodulation | suppressive |
| natalizumab | anti-α4-integrin; blocks leukocyte transit across the BBB | suppressive |
| fingolimod, ponesimod | S1P-receptor modulators; lymphocyte egress blocked | suppressive |
| teriflunomide | DHODH inhibitor; cytostatic to proliferating lymphocytes | suppressive |
| dimethyl fumarate | Nrf2 activation; immunomodulatory | suppressive |
| ocrelizumab | anti-CD20; B-cell depletion | suppressive |
| alemtuzumab | anti-CD52; lymphocyte depletion | suppressive |
| lenercept | TNF-receptor p55-IgG fusion; TNF blockade | suppressive **+ regulation-disrupting** (TNF-/- mice get severe EAE, PMID 9427610; TNFR2 needed for remyelination, PMID 11600888) |
| atacicept | TACI-Ig; blocks BAFF and APRIL, depletes plasma cells | suppressive |
| IFN-gamma | pro-inflammatory type II interferon | immunogenic |
| APL CGP77116 | altered peptide ligand of MBP 83-99; encephalitogenic in T-cell assays | immunogenic |

## What this table is for, and the result it forces

The lenercept and atacicept rows are the point. **Lenercept and atacicept are
immunosuppressive by mechanism and harmed patients in trials.** Any rule that maps
mechanism class to clinical direction — including the two-constant rule in
`bricks/grounding.py` — predicts benefit for both and is wrong on both.

That is not a defect in the table; it is the table doing its job. A gate scored only
on IFN-beta, glatiramer and one immunogenic peptide could not expose it, which is
exactly why `BUILD_PLAN.md` §8 blocker (3) calls a four-arm gate "too few to test
anything". Recording the two counterexamples is what turns the gate from a
restatement of the setup into a test the current rule can fail.

## MRI lesion anchors — blocker (4)'s non-circular channel

Per-drug potency must not be fitted on the arm's own relapse number; that is the
circularity recorded in BUILD_PLAN §8.4. The MRI channel avoids it, because the
gate predicts ARR and this table carries lesion outcomes.

**This table is INCOMPLETE and that is the current blocker on magnitudes.** Only
outcomes verified against the trial's own report are listed. Most of these trials
are behind paywalls (NEJM, Lancet) and their MRI numbers are not in the abstracts
Europe PMC serves, so they could not be transcribed here without access. Nothing
below was taken from a meta-analysis ranking or a secondary summary.

| arm | trial | comparator | MRI outcome | source |
|---|---|---|---|---|
| ocrelizumab | OPERA I | IFN beta-1a | **94% fewer** Gd-enhancing T1 lesions (0.02 vs 0.29 per scan) | NEJM 2017, PMID 28002679 (abstract) |
| lenercept | Lenercept MS Study Group | placebo | **no significant difference** in MRI outcomes, while relapse rate rose (p=0.006-0.007) and relapses were more severe and longer | Neurology 1999;53:457, PMID 10449104 |

The lenercept row is the important one and it is not an anchor, it is a warning:
it is the documented case where the MRI channel and the relapse outcome
disagree, and `bricks/sormani.py` flags that region as `blind_spot=True`. Any
potency fitted through MRI inherits it.

**Still to extract** (needs full-text access): PRISMS, CONFIRM, AFFIRM, FREEDOMS,
TEMSO, DEFINE, CARE-MS I, OPTIMUM. For each: new/enlarging T2 lesion count and
Gd-enhancing lesion count, treated arm and its comparator, as reported.
