# Exam C — pre-registration: the Cochrane 2024 network as a second exam

Written 2026-09-26, before `backtest/exam_cochrane.py` was first run. Nothing
below changes after the first run; a Result section is appended under it.
docs/FOUR_DAY_PLAN.md day 4 item 3 names this exam.

## Source

    Gonzalez-Lorenzo M, Ridley B, Minozzi S, Del Giovane C, Peryer G, Piggott T,
    Foschi M, Filippini G, Tramacere I, Baldin E, Nonino F. "Immunomodulators and
    immunosuppressants for relapsing-remitting multiple sclerosis: a network
    meta-analysis." Cochrane Database Syst Rev 2024;1:CD011381. PMC10765473.

The open-access XML at Europe PMC (fetched 2026-09-26). The outcome used is the
**network risk ratio, versus placebo, of people with one or more relapses over
24 months**, with its 95% CI, exactly as the review's results text states it.
The per-trial forest data behind the Cochrane Library's bot wall is NOT used;
if it is obtained later it is a different exam.

## Why this exam is different from exam v2

Exam v2 scores one trial per arm against that trial's comparator. This exam
scores a network estimate against placebo for every drug at once, on one
endpoint and one horizon, from one analysis. It has fewer arms and no
direction-only arms, so it cannot replace exam v2; it answers whether the
models' within-dial and between-dial ORDER survives a change of exam.

## Rules

**C1. Arms.** Every arm in `backtest/clinical.KNOWN_OUTCOMES` for which the
review's results text gives a network RR versus placebo for relapses over 24
months with a 95% CI. The arm names are matched by molecule; interferon
beta-1a (Avonex, Rebif) maps to the `IFN-beta` arm. Arms the review pools that
the arm set lacks are listed in the data module and not scored. Certainty
grades are recorded and not used to weight.

**C2. Outcome and interval.** Percent change = (RR - 1) x 100. The scored
interval is the 95% CI converted the same way. A prediction inside the CI has
error 0; outside, the distance to the nearer bound (exam v2's `dist`). This
is a risk ratio of people relapsing, not a rate ratio; the models predict a
rate ratio through Sormani. The mismatch is stated, not corrected: on a
24-month horizon with relapse probabilities well below one, the two ratios
move together, and no model here is close enough for the difference to
matter. If one ever is, this rule is revisited under a new registration.

**C3. Folds and null.** Leave-one-dial-out over the arms of C1, the dial
patterns being each model's own (`_pattern` in `backtest/exam_v2.py`). The
null predicts the mean of the training arms' RR point estimates as percent
change. Both scored with C2's interval error. An arm alone on its dial forms
its own fold.

**C4. Prediction.** Each model's cached exam v2 response curve
(`results/exam_v2_curve*.json`, unchanged, not rebuilt) at the potency fitted
by least squares on the training folds over exam v2's `FIT_GRID`. All arms are
placebo-controlled in the network, so no comparator adjustment.

**C5. Statistics.**

    CS1  interval-LOMO MAE vs null, per model, with the bootstrap CI on the
         paired fold gap (exam v2's `_paired_ci`)
    CS2  the oracle: best single value per dial group, vs the same null
    CS3  Kendall tau between predicted and observed percent change across
         arms, at the potency fitted on ALL arms (in sample, labelled so)

**C6. Falsifiers.** Same reading order as exam v2's R6: if CS2 does not beat
the null by the 20% margin of `gate/criterion.py`, this exam has no
resolving power on this arm set and the models are not blamed for it. If
CS2 clears and CS1 does not, the model is what fails. CS3 is descriptive.

**C7. What is not done.** No new response tables. No new arms wired. The
disability-worsening endpoint (gap G10) is not scored here.
