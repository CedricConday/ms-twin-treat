# The decision gate: can any of this be proposed to anyone?

Short answer, as of 2026-09-26: **no. The gate is red, and it is red for
three models now, not one.** The repo can say in numbers where each one
fails, and which of the three failures is the model's.

The single fact to take away, if you read nothing else: **every out-of-sample
scorer in this repository loses to predict-the-mean** — to answering with the
average of the other arms' outcomes, ignoring which drug it was asked about.
That was true of the transcription on 2026-09-21; on 2026-09-26 it is true of
the two candidate replacements as well, scored on the widened exam that was
built to give them the fairest possible hearing, and true again on a second,
independent exam (the Cochrane 2024 network) run the same day.

This document is for the reader who is not going to run the code. Everything
below is produced by the commands at the end, and every command re-measures.

## The exam, and what a perfect model would score on it

Exam v2 (`docs/EXAM_V2_PREREG.md`, registered 2026-09-25 before it was run)
scores all 23 treated arms, including the eight direction-only arms — the
three that did nothing, the four that harmed, the one that worked without a
rate — as intervals. Its oracle is a model that predicts one number per dial
group perfectly. The oracle sets the prize; the model has to collect it.

Definition of green (`docs/FOUR_DAY_PLAN.md`): S1 at or below 80% of its
null, the fold-gap CI excluding zero, the placebo-only subset clearing too.

## The three models, measured 2026-09-26

| model | dial patterns | S1 interval-LOMO | null | green needs | S2 oracle | headroom |
|---|---|---|---|---|---|---|
| Velez 2011 transcription (`bricks/qsp_velez.py`) | 8 | **30.7pp** | 20.0pp | 16.0pp | 9.1pp | 10.9pp (54%) |
| Jenner 2026 two-equation probe (`bricks/qsp_minimal.py`) | 4 | **43.8pp** | 30.6pp | 24.5pp | 12.7pp | 17.9pp (58%) |
| Pernice 2020 port, Figure S2 readings (`bricks/qsp_pernice.py`) | 11 | **35.1pp** | 19.8pp | 15.8pp | 6.3pp | 13.6pp (68%) |
| Pernice 2020 port, memory-read variant | 11 | **31.2pp** | 19.8pp | 15.8pp | 6.3pp | 13.6pp (68%) |

Every row loses. No fold-gap CI excludes zero. The nulls differ because the
folds are each model's own dial groups, and the null is drawn from the same
arms the model is scored on, in the same folds.

## The second exam, measured 2026-09-26

Exam C (`docs/EXAM_COCHRANE_PREREG.md`, registered before it was run) scores
the same cached response curves, unrebuilt, against a different truth: the
Cochrane 2024 network risk ratio versus placebo for relapses at 24 months,
ten arms, one analysis. It cannot replace exam v2 (fewer arms, no
direction-only arms); it asks whether the order survives a change of exam.

| model | dial folds | CS1 interval-LOMO | null | fold-gap CI | CS2 oracle | headroom |
|---|---|---|---|---|---|---|
| Velez 2011 transcription | 3 | **28.8pp** | 12.4pp | [-3.3, +31.9] | 0.8pp | 11.6pp (94%) |
| Jenner 2026 two-equation probe | 1 | **23.1pp** | undefined | — | 4.4pp | undefined |
| Pernice 2020 port, Figure S2 readings | 4 | **10.0pp** | 8.6pp | [-4.0, +13.3] | 0.8pp | 7.8pp (91%) |
| Pernice 2020 port, memory-read variant | 4 | **25.8pp** | 8.6pp | [-2.8, +23.0] | 0.8pp | 7.8pp (91%) |

Every model loses here too. The probe puts all ten arms on one dial, so it
has no training fold and no null: not blamed, not credited. The port on its
Figure S2 readings loses by 1.4pp with the CI straddling zero, and two of its
four dial folds beat their null (transit block, proliferation block) while
depletion and the glatiramer singleton lose. That is the closest any
out-of-sample scorer in this repository has come to its null. It is not
green: green is 80% of the null with the CI excluding zero, and this is 116%
with the CI including it. What it fails to collect is, again, per-dial
magnitude: depletion is under-predicted at the one shared potency that suits
transit block.

Read the last two columns first. **The exam is passable**: a perfect
dial-level model beats the null by 54% to 68% depending on how finely the
model divides the arms into dials, and the port's dial map — depletion,
effector deletion, transit block, activation block and proliferation block as
five separate patterns — is the best of the three by that measure. The prize
is real and it grew as the vocabulary got finer. The models did not collect it.

## Why each one loses

**The transcription** switches itself off: the least-squares potency is 0.00
in all eight folds, because three of its dial patterns point the wrong way
(effector loss predicts harm for six arms that worked; immunogenic challenge
predicts benefit; a costimulation block predicts a 43% benefit for an arm
that did nothing), and any positive potency is charged for them. It is a form
failure, not a noise failure: rebuilt at 1024 seeds instead of 128, every
number above is identical to the decimal (BUILD_PLAN §8.4, 2026-09-26).

**The two-equation probe** has load-driven damage by construction and still
gets the signs wrong, because its time-averaged inflammation is not monotone
in the disease strength: lowering the attack moves the patient from a limit
cycle onto a disease equilibrium with *higher* mean inflammation, and raising
it destroys myelin faster and leaves the inflammation less to feed on. Two of
the three acceptance tests fail structurally.

**The port** gets the signs mostly right — depletion, effector deletion and
transit block all reduce damage at every potency, immunogenic challenge
raises it — and loses anyway, because one shared potency has to serve a
transit-block curve at -74% and an activation-block curve at -24% against
truths between -30% and -68% on both; least squares settles near zero and
nine of sixteen working drugs score neutral. The costimulation block still
predicts benefit, and it shares its dial with the anti-CD20 arms. The port's
two-year behaviour does not reproduce the paper's two-year figures under
either reading of the memory arc; its 30-day figure reproduces on 10 of 12
landmarks. Both readings are scored and both lose, on both exams; the
memory-read variant's saturated untreated arm compresses every curve on
exam C and its in-sample tau turns negative.

## The one positive result, and what it is wired to

**S4.** A trial's own observed lesion ratio orders the arms inside a dial
group by relapse outcome: tau +0.80 on 10 pairs (p = 0.024) under the
transcription's dials, tau +1.00 on 7 pairs (p = 0.013) under the port's.
`docs/RECOVERABILITY.md` had found the MRI-fitted potency biased 1.55x high as
a magnitude; bias does not disturb rank. So the MRI channel is a **within-dial
rank input** and not a magnitude input, and `screen/report.py` now prints the
within-dial order of the real arms by that channel, reading the S4 statistic
from the exam's own artifact. It is never applied to a screened candidate,
which has no MRI trial. `rank_candidates()` still raises.

## What the device says

`gate.decide(candidate)` returns KILL, ABSTAIN or PASS. PASS requires a live
certificate that the predictor beats predict-the-mean out of sample. No
predictor in the repository holds one. **No input can return PASS**, on any of
the three models.

## What would change the verdict

Per-dial magnitude from an independent source is the missing input, on every
model, and it was the missing input on 2026-09-21 too. What changed is that the
port now has dials on which the sign is right, so a magnitude channel would
have something to scale. The MRI channel ranks within a dial but does not size
(1.55x bias). No other independent channel has been found. The 2024 Cochrane
network (`docs/EXAM_COCHRANE_PREREG.md`) was run as a second exam of the same
curves, not as a source of magnitudes, and it names the same missing input.
The per-trial forest data behind the Cochrane Library's bot wall is not used;
obtained, it would be a different exam under a new registration.

What has been ruled out, with the measurement: cohort size (G2, 1024 seeds,
no change); the carrying-capacity extension (inside the noise floor of the
transcription); widening the arm set from 12 to 15 (0.8pp to 0.9pp of
out-of-sample headroom on the point exam); the two-equation form (signs wrong
structurally); and the port as transcribed (signs right, magnitude absent,
on two exams).

## What can honestly be said to a third party today

- The harness runs end to end, is tested, and is reproducible from open
  data and open-access papers, three models deep.
- It reproduces the **direction** of known trial outcomes on some arms and
  not others, and reports which, per model.
- It **kills** candidates on structural grounds, and that judgement stands.
- The MRI channel **orders drugs within a mechanism class** and the repo
  says so with a p-value and a pair count.
- On a second, independent exam the best model comes within 1.4pp of
  predict-the-mean and still does not beat it; the repo says so with the CI.
- It **cannot rank** candidates or predict an effect size, on any model.
- No candidate has passed, and the criterion for passing was fixed and dated
  before any candidate or any model was scored against it.

What may not be said: that any simulated therapy is promising, that a ranking
exists, or that anything here is evidence about multiple sclerosis.

## Reproducing this page

```
PYTHONPATH=. python -m backtest.exam_v2                          # transcription
PYTHONPATH=. python -m backtest.exam_v2 --model minimal          # two-equation probe
PYTHONPATH=. python -m backtest.exam_v2 --model pernice          # the port, Figure S2 readings
PYTHONPATH=. python -m backtest.exam_v2 --model pernice --variant memread
PYTHONPATH=. python -m backtest.exam_cochrane                    # exam C, all four rows
PYTHONPATH=. python -m gate.device                               # the verdicts
PYTHONPATH=. python -m screen.report                             # the screen, with the G4 table
PYTHONPATH=. python scripts/state_of_build.py                    # every gate, one page
```

Figures on this page were measured 2026-09-26 on the 23-arm exam and the
10-arm Cochrane exam. If a number
here disagrees with what the code prints, **the code is right and this page
is stale**. `results/STATE.md` is regenerated, never edited; start there.
