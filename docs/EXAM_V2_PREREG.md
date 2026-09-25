# Exam v2 — pre-registered before any number was computed

**Written and committed 2026-09-25, before `backtest/exam_v2.py` was run.** The
rules below were fixed first so that the numbers at the end cannot have shaped
them. If this file and the code disagree, the code is a bug and this file wins.

## Why a second exam

The gate this repo scores today (`backtest/lomo.py`, criterion in
`gate/criterion.py`) grades a model on the **fifteen arms that report a point
annualised-relapse-rate change** and discards the other eight, because a point
error cannot be taken against "no significant difference" or "more relapses,
unquantified".

Measured 2026-09-25 on those fifteen arms, before writing this file:

| | |
|---|---|
| outcomes | mean −45.4%, SD 12.5pp, range −68 to −29 |
| between-dial-group spread | −37.1 (alpha_E), −50.9 (gamma_E), −51.3 (ke), −45.0, −29.0 (singletons) |
| within-group SD | 10.7 (alpha_E), 13.6 (gamma_E), 6.1 (ke) |
| LOMO oracle: predict every arm by its own group's true mean | **7.3pp** |
| LOMO null: predict every arm by the mean of the *other* groups | **11.8pp** |
| headroom for a perfect mechanism-level model | **4.5pp**, 38% relative |
| noise floor of one model prediction at 128 seeds (13.7% CV on ≈ −45%) | **≈ 6pp** |

So the current gate is passable in principle, by a model within about four
points of perfect, and its headroom is below the model's own noise floor. That
is a thin instrument, and it is thin for a reason the arm set cannot fix by
growing: every approved therapy lands between −30% and −68%, because that is
what gets approved. The outcomes that differ by *mechanism* — the therapies
that did nothing and the ones that made things worse — are exactly the eight
arms the point-error gate cannot read.

Eight direction-only arms, and the dial each sits on:

| arm | direction | dial |
|---|---|---|
| lenercept | harms | alpha_E\|alpha_R |
| atacicept | harms | gamma_E |
| IFN-gamma | harms | delta\|naive_E |
| APL CGP77116 | harms | delta\|naive_E |
| rituximab | improves | ke |
| ustekinumab | neutral | alpha_E |
| abatacept | neutral | delta |
| Tovaxin | neutral | gamma_E |

Two of them sit on `gamma_E` beside six arms that all worked, and one sits on
`alpha_E` beside four that all worked. A representation in which every arm on a
dial gets the same prediction cannot be right about both, and the current exam
never asks.

## Rules, fixed now

### R1. Every treated arm is scored, as an interval

| outcome type | interval | source of the width |
|---|---|---|
| point ARR change *x* | [*x*, *x*] | the trial |
| neutral, unquantified | [−15, +15] | `MAG_TOLERANCE` in `backtest/clinical.py`, the width the repo already treats as "within tolerance" |
| harms, unquantified | [+5, +∞) | `NEUTRAL_BAND`, the smallest change the repo already reads as a direction |
| improves, unquantified | (−∞, −5] | same |

Error of a prediction *p* against [lo, hi] is its distance to the interval:
0 inside, otherwise the distance to the nearer bound. For a point outcome this
is the ordinary absolute error, so the fifteen quantified arms score exactly as
they do today.

**Not adopted, and why:** atacicept's three dose-arm numbers (+108/+126/+158%)
are not used as a lower bound. The repo already declined to pick a number out
of that trial (`docs/TRIAL_ANCHORS.md`), and this file does not overrule it;
the harm interval is the one-sided one above, like every other unquantified
harm.

### R2. The null imputes a bound

Predict-the-mean needs a number per training arm. Each arm contributes:

| outcome type | value the null averages |
|---|---|
| point *x* | *x* |
| neutral | 0 |
| harms | +5 |
| improves, unquantified | −5 |

That is the nearest finite bound, so the null is charged for the direction-only
arms at their most favourable reading and is never charged for a width it
could not know.

### R3. The fold is the dial pattern, over all 23 treated arms

Leave-one-mechanism-out, as `backtest/lomo.py` does it: the fold is the tuple
of intervention points an arm moves (`bricks/profiles.touched_points`), and the
whole pattern is held out. This gives eight folds, three of which do not exist
in the current exam (`alpha_E|alpha_R`, `delta|naive_E`, `delta`) and need
their response-curve columns built the same way as the existing five, cached
in `results/exam_v2_curve.json`, never in the transcription's own file.

One global potency *s* is fitted on the training arms by least squares on the
interval distance, over the same grid `lomo.py` uses. Nothing else is fitted.

### R4. Active-comparator arms are predicted the way their trial ran

The primary figure adjusts: predicted change = (1 + drug)/(1 + comparator) − 1,
both simulated at the same potency, comparator taken from `bricks/profiles.py`.
`backtest/lomo.py` does not adjust (it predicts every arm against untreated,
while `backtest/loo.py` does adjust); the unadjusted figure is reported second,
for comparability with the existing headline, and is not the primary.

### R5. The four statistics, and what each one is allowed to say

**S1 — interval-LOMO.** Mean interval error over 23 arms, model vs null, per
fold and pooled, paired bootstrap over folds (10,000 resamples, 95% CI) on the
mean of (error − null error), same construction as `gate/evidence.py`. This is
the gate figure. It says whether the model carries mechanism information on the
widened exam.

**S2 — interval-LOMO oracle.** For each dial group, the single value that
minimises the summed interval distance of that group's arms (in-sample), scored
against the same null as S1. This is what a *perfect* dial-level model could
score. It says whether the widened exam is passable at all, which is the number
that decides whether the ranking goal survives.

**S3 — direction, three classes.** The S1 prediction read through
`NEUTRAL_BAND` (below −5 improves, above +5 harms, else neutral) against the
arm's direction. Accuracy and balanced accuracy, against the majority-class
null and a label-permutation null (10,000 shuffles; p is the fraction of
shuffles whose balanced accuracy is at least the observed one).

**S4 — does the independent channel rank within a dial?** Kendall's tau over
pairs of arms sharing a dial pattern, between the trial's observed lesion ratio
(`backtest/potency.OBSERVED_LESION_RATIOS`, never fitted) and its observed ARR
change, with a within-group permutation null. Reported twice: all same-dial
pairs, and only pairs whose lesion metric is the same string (the repo's own
warning that Gd-T1, new-T2 and CUAL ratios are not interchangeable). This does
not involve the ODE. It says whether *any* per-drug input this repo has can
order two drugs on the same dial, which is the question a screen's ranking
would eventually turn on.

### R6. Falsifiers

- If **S2** headroom is below 20% of the null (the margin in
  `gate/criterion.py`), then no model can pass this exam on these arms either,
  and the repo's ranking goal is dropped rather than the bar moved again.
- If **S2** clears 20% and **S1** does not, the model is what fails, and the
  gap between them is the prize a replacement model is playing for. That is
  the expected outcome.
- If **S3** balanced accuracy does not beat the majority null at p < 0.05, the
  model carries no direction information on the widened exam.
- If **S4** tau is not positive with p < 0.05 on the same-metric pairs, the MRI
  channel is not a within-dial ranking input, and the potency layer needs a
  different source before any within-dial ranking is attempted.

### R7. What this file does not do

It does not change `gate/criterion.py`, which stays frozen at 2026-09-20. It
does not replace `backtest/lomo.py`; both scorers run and both are reported. A
criterion for exam v2 is proposed only after S2 is known, and adopting it is
the operator's call, recorded with a date.

## Result

*(appended after running — nothing above this line changed)*

Run 2026-09-25 at 67f8b5f, `PYTHONPATH=. python -m backtest.exam_v2`, 128 seeds,
three new response columns built into `results/exam_v2_curve.json`
(`alpha_E|alpha_R`, `delta|naive_E`, `delta`). Full rows in `results/exam_v2.json`.

| statistic | model | null | verdict |
|---|---|---|---|
| S1 interval-LOMO, comparator-adjusted | **30.7pp** | 20.0pp | loses; fold-gap CI [−13.0, +20.2] |
| S1 placebo-controlled subset | 25.2pp | 21.2pp | loses |
| S1 unadjusted (lomo.py convention) | 30.7pp | 20.0pp | identical, see below |
| **S2 oracle, perfect dial-level model** | **9.1pp** | **20.0pp** | **headroom 10.9pp, 54% of the null** |
| S3 direction, accuracy / balanced | 13% / 0.33 | 70% / 0.33 | loses; permutation p = 1.000 |
| S4 MRI channel within-dial rank, all pairs (n=10) | tau **+0.80** | permutation | **p = 0.024** |
| S4 same-metric pairs only (n=2) | tau +1.00 | permutation | p = 0.247, underpowered |

**R6 falsifiers, read in order.** S2 clears the 20% margin by a wide
distance: the widened exam is passable, and the ranking goal survives. S1 does
not clear it, so the model is what fails, and the prize on this exam is
10.9pp against a 4.5pp prize on the point exam. That is the expected outcome,
and the widened exam is the one to score a replacement model on.

**Why the model scores what it scores: the fitted potency is 0.00 in all eight
folds.** The least-squares fit switched the model off. Three dial patterns
predict the wrong sign, and any positive potency is charged for them:

| pattern | arms | trial says | model at s = 0.5 |
|---|---|---|---|
| `gamma_E` | 6 improving, 1 harm, 1 neutral | mostly benefit | **+235%** (harm) |
| `delta\|naive_E` | IFN-gamma, APL | harm | **−12%** (benefit) |
| `delta` | abatacept | nothing | **−43%** |
| `alpha_R` | daclizumab | −45% | large harm |
| `alpha_E\|alpha_R` | lenercept | harm | +845% at s = 0.2, out of regime by 0.5 |

With one potency shared across dials, the arms the model gets backwards cost
more than the arms it gets right can earn, at every s above zero. On the point
exam this never showed, because the arms that expose it are the ones that exam
discards; there the fit landed at s ≈ 0.55. The unadjusted figure is identical
because at s = 0 there is nothing to adjust.

So S3 carries no information for the same reason: a model at s = 0 predicts
neutral for all 23 arms.

**S4 is the one positive result in this repository to date.** The trial's own
lesion ratio orders drugs within a dial group correctly on 8 of 10 pairs
(p = 0.024 against a within-group permutation). `docs/RECOVERABILITY.md` found
the MRI-fitted potency biased 1.55x high as a *magnitude*; bias does not
disturb rank. So the MRI channel is a within-dial *ranking* input even though
it is not a magnitude input. The same-metric restriction leaves two pairs and
cannot confirm or refute it, which is the caveat to carry.

**What this changes.** A replacement model is scored on exam v2, not on the
point exam, and has to do three things the current one cannot: keep the sign
of depletion, keep the sign of immunogenic challenge, and predict nothing for
a costimulation block. The Pernice 2020 scope (`docs/PERNICE_PORT_SCOPE.md`)
addresses the first; the other two are new acceptance tests for it.
