# Exam v5, registered 2026-11-23 before any new anchor was collected or any score computed

**What the registrant had seen.** Exams v3 and v4 per arm, and the ceiling analysis, which
shows the port can reach lenercept's and the altered peptide ligand's harm intervals and
at which magnitudes. As with v4, this is a repair on the same 23 arms: **a green is a
candidate needing arms the exam has not touched.** The collector of the anchors is not
told which values the outcomes ask for.

Exam v4 unchanged, plus anchors for two of the four harm arms, which v3 left as "none
defined" although `bricks/profiles_pernice.py` gives both a mechanism:

| arm | dials (from the profile, unchanged) | anchor |
|---|---|---|
| lenercept | pTeff_Dup, pTreg_Dup, pRemyelinization lowered | fraction of TNF bioactivity (or of soluble TNF) neutralised at trough at the trial's highest dose; a statement of complete neutralisation enters as 0.95 |
| APL CGP77116 | pTeff_Activation raised, multiplier 1 + a | (frequency of T cells reactive to the peptide or to MBP 83-99 on treatment) / (frequency before) - 1, capped at 0.95 |

IFN-gamma and atacicept stay unanchored (no quantity of this kind for the first, no rate
in the net for the second). Sources and rules as v3; values in `data/pd_anchors_v5.tsv`
in the v3 format; an anchor that cannot be sourced leaves the arm unanchored.

Run: `python -m backtest.exam_v4 --extra-anchors data/pd_anchors_v5.tsv --out results/exam_v5.json`.
One run, appended below whatever it is.

## Result, one run, 2026-11-23

**Red, by 0.2pp on the first criterion and clearly on the other two.** `results/exam_v5.json`.

| | S1 | null | share of null | fold-gap CI | placebo subset | anchored arms only |
|---|---|---|---|---|---|---|
| exam v5 | 16.0pp | 19.8pp | 81 % | [-12.6, +10.4] | 18.7 vs 21.5 (87 %) | 15.7 vs 22.0 (n = 14) |
| exam v4 | 17.6pp | 19.8pp | 89 % | [-6.6, +12.2] | 21.1 vs 21.5 (98 %) | 16.9 vs 21.0 (n = 13) |

Green is 80 % (15.8pp), a fold-gap CI excluding zero, and the placebo subset at 80 %.

- **APL CGP77116** sourced (T cells specific for the peptide rose from 0.80 to 210.93 per
  10^6, capped at 0.95): the port predicts +20.4 %, inside the harm interval, error 0
  against the null's 35.5. The port called a harm from the drug's pharmacodynamics alone.
- **Lenercept** could not be sourced: the MS trial measured no TNF bioactivity. A hedged
  statement of near-complete binding from a rheumatoid arthritis study at half the dose
  exists and would enter as 0.95, where the ceiling analysis says the port predicts harm.
  It was not entered, because the registrant knows that value helps. The arm scores the null.
