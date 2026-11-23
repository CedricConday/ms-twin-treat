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
