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
