# Exam H (confirmation), registered 2026-11-23 before any of its outcomes was collected

**Why.** Exams v3 to v5 were each designed after reading the previous one's misses on
the same 23 arms. Exam H scores trials none of that touched, with the pipeline frozen.

**Frozen pipeline.** Exam v5 exactly as at commit 85a32fc: Pernice port, Figure S2
readings, the anchor as magnitude, nothing fitted, comparator's observed effect from the
23-arm dataset (C1), interferons unanchored (C2), two anchors for an S1P modulator (C3).
The null for every arm is the mean interval anchor of the 23 exam arms (-28.8 % against
untreated; converted with C1 for an active comparator). Model changes made after today
are scored on Exam H as separate, later, labelled runs; the first run is this pipeline.

**Arms, fixed here by rule:** randomised trials in relapsing MS that report an annualised
relapse rate per arm, not among the 23 arms' source trials.

Stratum A, drugs the exams never used:
siponimod (BOLD, 2 mg, placebo), mitoxantrone (Millefiorini 1997, placebo), firategrast
(phase 2, highest dose, placebo), evobrutinib (evolutionRMS 1 and 2, teriflunomide),
tolebrutinib (GEMINI 1 and 2, teriflunomide), fenebrutinib (FENhance, teriflunomide, if
reported), peginterferon beta-1a (ADVANCE, every 2 weeks, placebo).

Stratum B, a drug of the exams in a different trial:
fingolimod (FREEDOMS II, placebo; TRANSFORMS, IFN beta-1a), teriflunomide (TOWER, placebo),
dimethyl fumarate (CONFIRM, placebo), ocrelizumab (OPERA II), ofatumumab (ASCLEPIOS II),
ublituximab (ULTIMATE II), alemtuzumab (CARE-MS II), ozanimod (SUNBEAM, IFN beta-1a),
daclizumab (SELECT, placebo), rituximab (HERMES, placebo).

A trial whose relapse rate cannot be sourced is dropped and named. Laquinimod is excluded:
no rate in the net can be assigned from its pharmacology.

**Profiles for stratum A**, by the conventions of `bricks/profiles_pernice.py`:
siponimod as fingolimod (both barrier-passage dials, C3 anchors); firategrast as
natalizumab (both passage dials, alpha-4 saturation at trough); mitoxantrone as
cladribine (both death rates, CD4+ reduction at 12 months, nearest time point); the BTK
inhibitors on pTeff_Activation lowered, as the other B-cell-directed drugs, with target
occupancy at trough as the anchor because they do not deplete (the abatacept rule);
peginterferon unanchored (C2). Stratum B drugs keep the anchors already on file.

**Sealing.** Outcomes go to `data/confirmation_outcomes.tsv`, collected by a process that
is given no predictions. The registrant does not open that file until the predictions are
committed in `results/exam_h_predictions.json`. Stratum A and B are reported separately
and together; S1 against the null, per-arm table, same green line (80 % of null).
One run. The registrant expects misses and names one in advance: the net has no B cell,
so the BTK inhibitors are predicted as if occupancy were efficacy.

## Result, one run, 2026-11-23 (predictions at 9835a4f, then `python -m backtest.exam_h score`)

**Not confirmed.** `results/exam_h.json`. Twenty arms scored (firategrast dropped by the
annualised-rate rule before the outcomes were opened); three unanchored and scored at the null.

| | n | S1 | null | share of null | gap CI |
|---|---|---|---|---|---|
| all arms | 20 | 26.5pp | 30.8pp | 86 % | [-20.6, +12.6] |
| stratum A, drugs never used | 9 | 36.9pp | 23.5pp | 157 % | [-16.1, +41.2] |
| stratum B, known drugs in other trials | 11 | 18.0pp | 36.8pp | 49 % | [-33.3, -4.7] |

- **Stratum B holds, with a caveat about its null.** Dimethyl fumarate in CONFIRM is
  predicted to 0.1pp, ofatumumab in ASCLEPIOS II to 1.7pp, alemtuzumab in CARE-MS II to
  5.0pp. But for a head-to-head trial the null says "no better than the comparator", which
  a drug chosen for a head-to-head trial usually beats; the 49 % flatters the model.
  Against placebo the model loses to the null on three of five (fingolimod, teriflunomide,
  rituximab).
- **Stratum A fails, where the registration said it would.** Evobrutinib and tolebrutinib
  were predicted 56 to 62 % better than teriflunomide from their target occupancy; they
  were not better at all (errors 55 to 69pp). Fenebrutinib, same class and nearly the same
  occupancy, was better by about as much as predicted (errors 3 and 11pp). The model cannot
  tell them apart: occupancy is not what separates them, and the net has no B cell.
- **The S1P overshoot repeats** on siponimod, both fingolimod trials and ozanimod
  (errors 25 to 41pp), as in exams v3 to v5.
- Weak sources, as the collector flagged: siponimod and mitoxantrone outcomes from reviews,
  fenebrutinib company-reported.

So the 81 % on the 23 exam arms does not carry to drugs of a mechanism the model has not
seen. It does carry, roughly, to new trials of drug classes whose pharmacodynamics it has.
