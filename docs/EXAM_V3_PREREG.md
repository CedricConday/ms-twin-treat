# Exam v3, registered 2026-11-23 before any anchor was collected or any score computed

**Question.** Does the Pernice port predict trial outcomes out of sample when each drug's
magnitude comes from a pharmacodynamic measurement made outside the efficacy outcome,
with nothing fitted to relapse data?

**What the registrant had seen.** `backtest/ceiling_per_drug.py` (BUILD_PLAN §8.4,
2026-11-23), including the magnitude each outcome asks for. That cannot be unseen, so the
primary score below has no free parameter: there is nothing to tune towards those values.
The anchor definitions were written before any anchor value was looked up.

## Arms, truth, null, green

The 23 treated arms of exam v2, their intervals and anchors as in `backtest/exam_v2.py`,
the same dial-pattern folds, the same null (mean anchor of the arms outside the fold).
Green as before: S1 at or below 80 % of its null, the fold-gap bootstrap CI excluding
zero, and the placebo-controlled subset clearing 80 % of its null too.

## The model side

Pernice 2020 port, Figure S2 readings, two-year schedule, lesion-load ratio through
`bricks/sormani`, exactly as exam v2, with the two form changes of 2026-11-23:

- D1: an active comparator is simulated at its own magnitude, not the tested drug's.
- D2: a dial the profile raises (the death rates) is parameterised by the steady-state
  count reduction `a`, multiplier `1 / (1 - a)`. A dial the profile lowers keeps
  multiplier `1 - a`.

## Magnitude: the anchor

One number `a` in [0, 0.95] per drug, the fractional reduction of the quantity the drug's
dial stands for, at the dose of the scored trial:

| dial class | arms | anchor |
|---|---|---|
| death rates raised | alemtuzumab, cladribine | reduction from baseline of circulating CD4+ T cells at 12 months (nearest reported time point) |
| effector death raised | Tovaxin | reduction of circulating myelin-reactive T-cell frequency after the vaccination course |
| barrier passage lowered | natalizumab | alpha-4 integrin saturation at trough, 300 mg every 4 weeks |
| barrier passage lowered | fingolimod, ponesimod, ozanimod | mean steady-state reduction of absolute lymphocyte count at the trial dose |
| activation lowered | ocrelizumab, ofatumumab, ublituximab, rituximab | reduction of circulating CD19+ B cells at steady state |
| activation lowered | abatacept | CD80/CD86 saturation at trough at the trial dose |
| proliferation lowered | IFN-beta, IFN-beta-1b, teriflunomide, dimethyl fumarate, ustekinumab | mean reduction from baseline of absolute lymphocyte count at steady state |
| none defined | glatiramer acetate, lenercept, APL CGP77116, IFN-gamma, atacicept, daclizumab | unanchored |

Source order: prescribing information (FDA label or EMA SmPC), else a peer-reviewed
pharmacodynamic study. A reported range enters as its midpoint; a value above 0.95 enters
as 0.95. Each value is stored with its source in `data/pd_anchors.tsv`. No efficacy
outcome (relapse, lesion, disability) is consulted for an anchor. A drug whose anchor
cannot be sourced is unanchored and is said to be.

## Scores

- **V3-P (primary), no fitted parameter.** Magnitude = anchor. An unanchored arm is
  predicted by its fold's null and so scores the null's error. A comparator without an
  anchor makes the arm unanchored.
- **V3-K (secondary).** Magnitude = min(0.95, kappa x anchor), one kappa for all drugs,
  fitted by interval least squares on the arms outside the held fold (grid 0 to 2).
- Both are also reported on the anchored arms alone, with the null on the same arms.

Reported whatever comes out: S1, null, fold-gap CI, placebo subset, per-arm table.
A red result is a result. Nothing in this file changes after the first run; a correction
is a dated section appended below.

## Result, first and only run, 2026-11-23 (`python -m backtest.exam_v3`, anchors at f379a2b)

**Red.** Numbers in `results/exam_v3.json`.

| score | S1 | null | share of null | fold-gap CI | placebo subset | anchored arms only (n = 15) |
|---|---|---|---|---|---|---|
| V3-P, anchor as magnitude, nothing fitted | 21.7pp | 19.8pp | 109 % | [-3.5, +13.5] | 22.4 vs 21.5 | 21.5 vs 18.6 |
| V3-K, one kappa fitted out of fold | 25.1pp | 19.8pp | 126 % | [-0.8, +12.6] | 23.3 vs 21.5 | 26.6 vs 18.6 |
| exam v2, shared potency (2026-09-26) | 35.1pp | 19.8pp | 177 % | | | |

Green needs 80 %. The unfitted score is the better of the two and is the closest this
exam has been to its null, but it loses. Fifteen arms were anchored; abatacept and
IFN-beta-1b could not be sourced and six had no anchor defined, so eight arms score the
null's error by rule.

Beats its null on eight arms: dimethyl fumarate, ustekinumab, natalizumab, fingolimod
(narrowly), ofatumumab, ublituximab, alemtuzumab, cladribine. Loses on six, in three groups:

- **The anchor is not the drug's mechanism.** IFN-beta: a 1.6 % fall in lymphocyte count
  (Avonex arm) predicts -3 % against -30 % observed. Lymphocyte count was the registered
  anchor for the proliferation class; for interferon it measures nothing the drug does.
  The error propagates to ozanimod and ocrelizumab, whose comparator it is.
- **The port is too sensitive to barrier passage.** The S1P modulators' measured 55 to
  70 % lymphocyte reduction predicts -76 to -84 % where -30 to -55 % was observed; the
  outcomes ask for 0.27 to 0.31 on this dial (ceiling analysis). Natalizumab errs the same
  way by less. Teriflunomide as comparator, under-predicted, widens the ponesimod miss.
- **Pharmacodynamics that did not become efficacy.** Tovaxin: a 79 % fall in
  myelin-reactive T cells predicts -78 %; the trial found nothing. Error 63pp, the largest
  on the table.

Kappa does not transfer: the folds fit 0.48 to 0.88, and scaling the anti-CD20 anchors by
0.48 costs more than it buys anywhere else.
