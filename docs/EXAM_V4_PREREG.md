# Exam v4, registered 2026-11-23 before any new anchor was collected or any score computed

**What the registrant had seen.** Exam v3's per-arm table (`docs/EXAM_V3_PREREG.md`,
Result) and the ceiling analysis, including the magnitude each outcome asks for. Every
change below was chosen after reading where v3 lost. So v4 is a repair of v3 on the same
23 arms, not a fresh out-of-sample test: **a green here is a candidate, to be confirmed
on arms this exam has not touched, and is reported in those words.** The new anchor
values are collected by a process that is not told which values the outcomes ask for.

Arms, intervals, folds, null and the green line are exam v2's and v3's. The model side is
v3's (Pernice port, D1 own-magnitude comparators replaced by C1 below, D2 count-reduction
parameterisation). V3-P style only: the anchor is the magnitude, nothing is fitted.

## Three changes, each with its reason

**C1. An active-comparator trial is converted with the comparator's observed effect.**
The model predicts the tested drug against untreated. To express that against an active
comparator it is divided by the comparator's effect against placebo. v3 took that effect
from the model; v4 takes the comparator's own arm in this dataset (its observed interval
anchor) whenever that arm lies outside the held fold, and falls back to the model, as in
v3, when it lies inside. Reason: the comparator's effect is a known input at the time a
head-to-head trial is designed, and a held-out drug should not inherit another drug's
prediction error. The null already uses the same out-of-fold outcomes.

**C2. A drug is anchored only if the anchor measures the quantity its dial stands for.**
Lymphocyte count was registered in v3 for the whole proliferation class. For a cytostatic
or lymphotoxic drug (teriflunomide, dimethyl fumarate) the count is that quantity. For
the interferons it is not: they act through signalling, and the count barely moves. Both
interferons are therefore unanchored in v4 and score the null's error, like glatiramer.
This removes v3's largest avoidable miss by rule and is the change most open to the
charge of being chosen with the answer in view; it is listed here so that it can be.

**C3. The S1P modulators get two anchors, one per cell type the dials carry.**
v3 put the total lymphocyte reduction on both barrier-passage dials. S1P receptor
modulation retains CCR7-positive cells in lymph nodes; effector memory T cells are
CCR7-negative. The cells that carry the effector passage dial are therefore not the ones
the total count follows. v4 anchors, per drug at the trial dose:

| dial | anchor |
|---|---|
| pPass_BBB_teff | fractional reduction from baseline of the absolute count of circulating CD4+ effector memory T cells (CCR7-negative, CD45RA-negative) at steady state |
| pPass_BBB_treg | fractional reduction from baseline of the absolute count of circulating regulatory T cells (CD4+ CD25+ FoxP3+ or CD127-low) at steady state |

A negative reduction (an increase) enters as 0. If a drug has one of the two and not the
other, or neither, it is unanchored. Natalizumab keeps v3's single anchor on both dials.
Sources and rules as v3 (`data/pd_anchors_v4.tsv`); no efficacy outcome is consulted.

Everything else is v3: abatacept and daclizumab stay unanchored, Tovaxin keeps its anchor.

## Scores

S1 over 23 arms against the null, the fold-gap CI, the placebo-controlled subset, the
anchored arms alone, and the per-arm table. One run. Appended below, whatever it is.

## Result, one run, 2026-11-23 (`python -m backtest.exam_v4`, anchors at the commit before this one)

**Red, and under its null for the first time.** `results/exam_v4.json`.

| | S1 | null | share of null | fold-gap CI | placebo subset | anchored arms only (n = 13) |
|---|---|---|---|---|---|---|
| exam v4 | 17.6pp | 19.8pp | 89 % | [-6.6, +12.2] | 21.1 vs 21.5 (98 %) | 16.9 vs 21.0 |
| exam v3, V3-P | 21.7pp | 19.8pp | 109 % | [-3.5, +13.5] | 22.4 vs 21.5 | 21.5 vs 18.6 |

Green needs 80 %, the CI excluding zero and the placebo subset at 80 %; none of the three
holds. By change:

- **C1 helped.** With the comparator's observed effect, alemtuzumab is predicted exactly
  (-55.0 against -55), ublituximab to 3.8pp, ofatumumab to 11.7pp, ocrelizumab to 16.5pp.
- **C2 did what a rule does**: the interferons score the null's error (1.2 and 5.1pp).
- **C3 failed.** Effector memory T cells are not spared in absolute counts: fingolimod
  0.83 (regulatory 0.93), ozanimod 0.64 (0.65). The predictions moved further from the
  outcomes (fingolimod -89 against -55, error 34.1, worse than v3's 28.5 and worse than
  the null). Ponesimod could not be sourced at 20 mg and scores the null. So the S1P miss
  is not the anchor; the port turns a given fall in passage into too large a fall in
  damage, for natalizumab as well (-90 against -68).

Eight arms are unanchored and carry 163 of the 405pp total error at the null's own level;
four of them are the harm arms (35.5pp each), which the port has profiles for.
