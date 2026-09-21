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
| ofatumumab | ASCLEPIOS I 2020 | teriflunomide | 0.11 vs 0.22 | −50% | improves | 32757523 |
| daclizumab | DECIDE 2015 | IFN beta-1a | 0.22 vs 0.39 | −45% | improves | 26444729 |
| cladribine 3.5 mg/kg | CLARITY 2010 | placebo | 0.14 vs 0.33 (96 wk) | −57.6% | improves | 20089950 |
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
| ofatumumab | anti-CD20; B-cell depletion (same target as ocrelizumab) | suppressive |
| daclizumab | anti-CD25 (IL-2Rα); cuts Treg numbers ~50%, raises IL-2 availability, expands CD56bright NK | suppressive |
| cladribine | purine nucleoside analogue; selectively lymphotoxic | suppressive |
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

## Why these three were added (2026-09-20)

Growing the arm set only helps the leave-one-MECHANISM-out test if the new arms
add mechanism COVERAGE. Six more depleting drugs would all land on `gamma_E` and
the fold count would not move. These three were chosen for what they do to the
groups, not for how many they add:

- **ofatumumab** shares anti-CD20 with ocrelizumab, turning `ke` from a
  single-arm fold into a real one. It also creates a deliberate test: the two
  share a dial AND a magnitude, so any difference their trials show is a
  difference this model provably cannot produce.
- **daclizumab** is a NEW mechanism — anti-CD25, cutting Treg numbers ~50% — and
  it is the most discriminating arm in the set. The model sees only the Treg
  loss and none of the compensating CD56bright NK expansion it has no cell type
  for, so it should predict HARM. DECIDE reported a 45% relapse reduction. An
  arm the model is expected to fail for a stateable reason is worth more than
  one it is expected to pass.
- **cladribine** adds a placebo-controlled arm to the depleting class, which was
  otherwise carried by active-comparator trials.

Arm set: 17 arms, 12 quantified, 5 mechanism groups.

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

---

## Candidate additions — researched 2026-09-21, not yet wired

The exam is thin, and that is the binding limit on what any model can demonstrate
here: the out-of-sample dial-level headroom is **0.8pp** (`gate/headroom.py`),
because twelve arms in three multi-member dial groups leave predict-the-mean
almost nothing to lose by. Widening the arm set is the cheapest available lever —
cheaper than any model port.

These rows are **research output only**. Nothing below is wired into
`bricks/profiles.py`, `bricks/intervention.py` or `backtest/clinical.py`; the
dial assignment is a separate call made from pharmacology, and the reasoning here
is an input to it, not a decision.

What the arm set is starved of, in priority order:

1. **More arms *per dial*, not more arms.** The out-of-sample variant can only
   score groups with n ≥ 2, so a third `ke` arm is worth more than a sixth
   `gamma_E` arm.
2. **The singletons are the cheapest win.** `glatiramer` (`alpha_R|delta`) and
   `daclizumab` (`alpha_R`) are fitted exactly for free and drop out of the
   honest variant entirely. Anything landing on those dials converts a free-fit
   arm into a scored one and adds two rows, not one.
3. **Failures are worth more than successes.** The set is almost all winners, so
   the outcomes cluster and the null is hard to beat. Arms that did nothing widen
   the spread the null has to cover.

### Quantified — ARR for both arms, from the trial's own report

| candidate | trial | comparator | ARR (drug vs comp.) | change | PMID | proposed dial, from pharmacology |
|---|---|---|---|---|---|---|
| ublituximab | ULTIMATE I 2022 | teriflunomide | 0.08 vs 0.19 | **−57.9%** | 36001711 | `ke` — anti-CD20 B-cell depletion, the same target class as ocrelizumab and ofatumumab, which `bricks/profiles.py` already routes to `ke` on the Martinez-Pasamar measurement |
| interferon beta-1b 8 MIU | IFNB MS Study Group 1993 | placebo | 0.84 vs 1.27 (annual exacerbation rate) | **−33.9%** | 8469318 | `alpha_E` — type I interferon, the same pharmacology as IFN beta-1a, which is already on `alpha_E` |
| ozanimod 1.0 mg | RADIANCE 2019 | IFN beta-1a | 0.17 vs 0.28 (RR 0.62) | **−39.3%** | 31492652 | `gamma_E` — S1P receptor modulator, lymph-node egress block, the same lump as fingolimod and ponesimod |
| laquinimod 0.6 mg | ALLEGRO 2012 | placebo | 0.30 vs 0.39 | **−23.1%** | 22417253 | **no confident dial** — see below |

**Value to the exam.** Ublituximab is the highest-value single row available: it
takes `ke` from 2 arms to 3, which is the group most starved in the out-of-sample
variant. IFN beta-1b takes `alpha_E` from 3 to 4 and costs nothing to justify,
since it is the same molecule class as an arm already assigned. Ozanimod grows
`gamma_E` from 5 to 6, which is the least useful addition — that group is already
the largest and is the one no potency can reach.

**Laquinimod is listed and deliberately not assigned.** Its mechanism — aryl
hydrocarbon receptor agonism shifting myeloid cells toward an anti-inflammatory
phenotype — does not map cleanly onto any of the Vélez model's named rates. It is
the weakest quantified effect available (−23.1%), which is exactly what the arm
set needs for spread, but forcing it onto a dial to get that spread would be
fitting the representation to the outcome. Assign it only if a defensible dial can
be argued from its pharmacology; otherwise leave it out.

**SUNBEAM is the second ozanimod trial** (PMID 31492651, ARR 0.18 vs 0.35 on IFN
beta-1a, −48.6%). One arm per drug is this table's existing convention — ASCLEPIOS
I was taken over ASCLEPIOS II on the same basis — so RADIANCE is proposed and
SUNBEAM recorded rather than added. The two disagree by 9pp on the same drug,
which is a useful measure of how much of any arm's number is trial noise.

### Direction-only — real randomised results, no ARR to score

These carry a direction and no magnitude, like `lenercept` and `IFN-gamma`
already do. They are the ones that widen the *failure* side of the exam.

| candidate | trial | comparator | result | direction | PMID | proposed dial |
|---|---|---|---|---|---|---|
| ustekinumab | phase II 2008 | placebo | no significant reduction in cumulative new Gd-enhancing lesions at any of four doses | **no effect** | 18703004 | `alpha_E` — anti-IL-12/23 p40, blocking Th1/Th17 differentiation, i.e. effector activation |
| abatacept | ACCLAIM 2017 | placebo | no significant difference in new Gd+ lesions or any clinical measure; enrolment closed early at 65 of 123 | **no effect** | 27481207 | `alpha_E` *with a caveat* — CTLA4-Ig blocks costimulation of effectors but also acts on regulatory populations, so the sign on `alpha_R` is not clean |
| secukinumab | proof-of-concept 2016 | placebo | primary endpoint missed: CUAL reduced 49%, CI −10 to 77, p=0.087; Gd+ lesions −67% (p=0.003), n=73 | **ambiguous** | 27142710 | `alpha_E` — anti-IL-17A |
| rituximab | RIFUND-MS 2022 | dimethyl fumarate | 3% vs 16% of patients relapsed, risk ratio 0.19 | **improves** | 35841908 | `ke` — anti-CD20, as above |

**A caution on RIFUND-MS.** Its outcome is the *proportion of patients who
relapsed*, not an annualised relapse rate. Those are different quantities and
converting one to the other would invent precision, so it belongs in this table
and not the quantified one — the same rule that keeps the active-comparator arms
un-chained.

**Ustekinumab and abatacept are the most valuable rows here** despite carrying no
magnitude, because `alpha_E` currently contains three arms that all worked. A dial
whose training data is all successes cannot teach a model that the dial sometimes
does nothing.

### Researched and rejected

**Opicinumab (AFFINITY Part 1, PMID 41454463, *Mult Scler* 2026).** Real,
randomised, and it missed: adjusted mean difference on ODRS 0.15, CI −0.05 to
0.35, p=0.148. Excluded anyway, for two reasons that are both about this model
rather than about the trial. Its target is LINGO-1 — a remyelination agent, acting
on oligodendrocyte differentiation, which is not an immune rate the Vélez model
represents at all; there is no dial it could honestly land on. And it was given as
an add-on to a background DMT, so its comparator is placebo-plus-DMT rather than
placebo, which is not a comparator any arm in this table uses.

**Mitoxantrone (MIMS, PMID 12504397).** A *progressive* MS population, not RRMS.
Every arm in this table is relapsing-remitting; mixing populations would make the
gate's arms non-comparable.

### A note on how these were verified, which is itself a result

Five candidate PMIDs were checked from memory before being searched. **Four of the
five pointed at entirely unrelated papers** — a note on relative afferent
pupillary defect, a diabetes-diagnosis review, a childhood cholera trial, and a
selenium-supplementation study. Only the 1993 interferon beta-1b PMID was correct.

Every identifier in the tables above was resolved by live query against Europe PMC
and NCBI E-utilities, and every ARR was read out of the retrieved abstract rather
than recalled. A recalled PMID is not a citation; it is a plausible-looking number
that happens to index something else.
