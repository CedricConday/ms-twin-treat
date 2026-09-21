# The decision gate: can any of this be proposed to anyone?

Short answer, as of 2026-09-21: **no, and the repo can now say exactly why, in
numbers rather than in caution.**

The single fact to take away, if you read nothing else: **every out-of-sample
scorer in this repository loses to predict-the-mean** — to answering with the
average of the other arms' outcomes, ignoring which drug it was asked about.

This document is for the reader who is not going to run the code — a clinician,
a collaborator, a company — and who wants to know what the simulator is entitled
to claim. Everything below is produced by `python3 -m gate.device` and
`python3 -m gate.ceiling`, and both re-measure rather than recall.

## What the device does

`gate.decide(candidate)` returns one of three verdicts:

| verdict | meaning |
|---|---|
| **KILL** | doomed for a reason needing no predicted effect size: the model runs away on it, the dial cannot help at any strength, or it is indistinguishable from an existing drug |
| **ABSTAIN** | survived the filters, and the model is not entitled to say anything further |
| **PASS** | survived, *and* the predictor that scored it has been shown out-of-sample to carry drug-specific information |

Run over the 14 single-dial candidates the model admits: **13 KILL, 1 ABSTAIN,
0 PASS.** No input can currently return PASS.

Every verdict names the model that produced it — `velez2011` today — and that is
a required field, not a note. Survival is not a model-independent fact: the
published transcription and the K = 2000 carrying-capacity extension score
45.4pp and 48.0pp on the current exam — a 2.6pp gap, inside this repo's own ~10%
noise floor for these figures, so nothing here can prefer one model — and they
disagree about 11 of 40 survivors. The label carries the *parameter*
too, because K = 50000 and K = 2000 are both "the extension" and behave
oppositely (+109% and −3% on the same dial). The device also refuses to issue a
verdict at all if its kill filters and its certificate turn out to describe
different models.

## Why PASS is unavailable

PASS requires an `EvidenceCertificate` — a live measurement that the predictor
beats the simplest possible competitor, "answer with the average of the other
arms." Both of the repo's out-of-sample scorers lose to it:

| scorer | what it holds out | MAE | predict-the-mean null | 95% CI on the paired gap |
|---|---|---|---|---|
| `backtest/loo.py` | one arm | 30.0pp | 11.3pp | [+8.9, +29.3] |
| `backtest/lomo.py` | one whole mechanism | 45.4pp | 11.8pp | [+20.3, +42.8] |

Both confidence intervals lie entirely **above** zero. This is not a near miss
that a recalibration closes: the predictor is reliably worse than the null.

A pass/fail rule built on such a predictor would be a statement about the
training arms' average, not about the therapy in front of it. The device
therefore refuses to make one, and the refusal is enforced by measurement — not
by a flag someone could flip.

## What would have to change, and what has been ruled out

The obvious repair is per-drug potency from independent data (`BUILD_PLAN`
blockers 4 and 6). `gate/ceiling.py` tested whether that repair can work, by
giving every arm its own potency fitted **to its own answer** — a ceiling no
independent source can beat. Two results:

**1. Nearly half the exam is unreachable.** Seven of fifteen quantified arms —
natalizumab, fingolimod, alemtuzumab, ponesimod, daclizumab, cladribine and
ozanimod — fit at the grid edge, meaning no potency produces a benefit at all.
That is the entire depleting / sequestering / trafficking class, which the ported
QSP model collapses onto a single loss term it cannot turn downward. Oracle error
there is 50.1pp against an 11.5pp null. **No potency data of any quality touches
this**; it is a direction failure, not a calibration one.

**2. On the reachable half, the independent channel misses.** The oracle fits
those eight arms to 0.4pp — which one free parameter per target does by
arithmetic, and establishes only that the response curve passes through each
trial's answer somewhere. Scoring the arms that also carry an MRI-fitted potency
from `backtest/potency.py`, which never sees a relapse number, gives **21.4pp
against a 10.4pp null.** It loses, and it is biased high every single time —
median 1.55x the oracle potency, range 1.08x to 2.46x.

The repo's existing cross-check on ocrelizumab looks like it agrees, and on
inspection does not: its 2.1x is in multiplier space and this table's 2.00 is in
potency space, and the two ratios share the MRI estimate as a term. Read in one
consistent parameterisation the two independent estimates are **4.0x apart**,
and the potency that reproduces the trial sits halfway between them on a log
scale — 2.00x above the EAE fit, 2.00x below the MRI fit. Neither independent
source finds it, and they miss in opposite directions.
See [RECOVERABILITY.md](RECOVERABILITY.md) for the full table and the two
distinct senses of "unreachable".

**Conclusion: blockers (4) and (6) are closed routes on this model, not deferred
ones.** More MRI-channel potency data cannot lift a gate whose independent
channel misses the potency it would need, on the half of the arm set the model
can reach at all. The remaining blocker is the model's *form* — specifically
that it has no CNS compartment and no trafficking, so depletion, sequestration
and transit blockade are one lumped dial (`bricks/profiles.py` documents this
and names the only way the lump has ever broken: someone measured what the drug
does to a parameter the model already has).

## What would change the verdict — and how little is on offer

Before budgeting a model port, it is worth knowing how much there is to win.
Predict every arm by the mean of its own dial group. That is a *perfect*
dial-level model: no fitting, no simulation, nothing to get wrong.

| scoring | MAE | null | headroom |
|---|---|---|---|
| in sample, all 15 arms | 7.3pp | 10.5pp | 3.3pp |
| singleton groups dropped | 8.4pp | 10.7pp | 2.3pp |
| **out of sample, within group** | **10.7pp** | **11.5pp** | **0.9pp** |

Each row's null is drawn from the same arms the model is scored on, and scored
the same way; pooling it over arms the model was not examined on would hand the
null free information.

Glatiramer and daclizumab are alone on their dials, so the first row fits them
exactly by construction. The last row is how every other scorer in this repo is
graded, and it is the honest bound.

**The entire prize is 0.9pp**, before anything is charged for simulation or
fitting. A replacement model does not need to be *better* than this one — it
needs to land within a percentage point or two of perfect to clear a
predict-the-mean null here. The measured mechanism-holdout error is 45.4pp, so
the model sits about 39pp from a ceiling that is itself almost level with the
floor.

That is a property of the **arms**, not of the model. Fifteen quantified arms in
a handful of dial groups is a thin exam, and widening it is far cheaper than any
model port.

**This was tried, and it barely moved.** Six arms were researched, verified and
wired on 2026-09-21, taking the set from 12 quantified arms to 15. The
out-of-sample prize went from 0.8pp to **0.9pp**. The lever is real and it is
much weaker per arm than it sounds, because the arms that exist to be added land
in the groups that are already largest.

And the groups that most need widening cannot be. The two singleton dials —
`alpha_R` and `alpha_R|delta`, holding daclizumab and glatiramer — are fitted for
free and drop out of the honest variant entirely, so a second quantified arm on
either is worth more than any number of additional anti-CD20s. There is no such
arm: no MS therapy whose primary mechanism is regulatory-T-cell restoration has a
published relapse rate against a stated comparator (`docs/TRIAL_ANCHORS.md`,
searched 2026-09-21). **Those groups are singletons in the arm set because they
are nearly singletons in the field.** The part of the exam that most needs
widening is the part that cannot be widened from the literature.

(The dial-level bound originates in `scripts/dial_ceiling.py`; the out-of-sample
variant and the singleton caveat are `gate/headroom.py`'s.)

## One thing both lanes got wrong on the same night

Two independent caches in this repository were found, two hours apart, to be
answering questions about an exam that no longer existed. The gate's
mechanism-holdout cache had been measured over twelve arms and was still being
read after fifteen were wired; the screen's survivor list had been computed
against the dials the arms occupied *before* three more arms were wired onto
them, which changes which candidates count as duplicates of an existing drug.

Neither was detected by a failing test. Both were found by a person reading, and
in both cases the stale artifact was indistinguishable from a fresh one — it
parsed, it carried a date, it carried a commit.

The rule that came out of it, now enforced in code on both sides: **a cached
measurement carries the exam it was taken on, and is refused when that exam has
changed.** It is worth stating on this page because it is the failure most likely
to recur, and the one a reader has no way to see.

## What can honestly be said to a third party today

- The harness runs end to end, is tested, and is reproducible from open data.
- It reproduces the **direction** of known trial outcomes on some arms and not
  others, and reports which.
- It **kills** candidates on structural grounds, and that judgement stands on
  its own evidence.
- It **cannot rank** candidates or predict an effect size. `screen.rank_candidates()`
  raises rather than returning a list, by design.
- No candidate has passed, and the criterion for passing was fixed and dated
  before any candidate was scored against it.

What may not be said: that any simulated therapy is promising, that a ranking
exists, or that anything here is evidence about multiple sclerosis. It is a
simulator that currently loses to an average, and it says so.

## Reproducing this page

```
PYTHONPATH=. python3 -m gate.device      # the verdicts and the certificate
PYTHONPATH=. python3 -m gate.ceiling     # the ceiling and the recoverability table
PYTHONPATH=. python3 -m gate.headroom    # which of the three failures binds
PYTHONPATH=. python3 -m backtest.loo     # arm holdout
PYTHONPATH=. python3 -m backtest.lomo    # mechanism holdout (minutes)
```

Figures on this page were measured 2026-09-21 on the 15-quantified-arm set. If a
number here disagrees with what the code prints, **the code is right and this
page is stale** — that has already happened once, when six arms were wired
underneath it.

One word on what "verified" means here, because this page is the one most likely
to be read by someone who will never run the code. Every trial number is real and
cited, and every identifier was resolved live rather than recalled. That is
*verified*, and it is not *validated*: it says the citations are sound, not that
the simulator is right about multiple sclerosis. The simulator currently loses to
an average.
