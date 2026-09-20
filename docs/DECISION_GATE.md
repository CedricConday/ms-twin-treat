# The decision gate: can any of this be proposed to anyone?

Short answer, as of 2026-09-20: **no, and the repo can now say exactly why, in
numbers rather than in caution.**

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

Run over the 14 single-dial candidates the model admits: **12 KILL, 2 ABSTAIN,
0 PASS.** No input can currently return PASS.

## Why PASS is unavailable

PASS requires an `EvidenceCertificate` — a live measurement that the predictor
beats the simplest possible competitor, "answer with the average of the other
arms." Both of the repo's out-of-sample scorers lose to it:

| scorer | what it holds out | MAE | predict-the-mean null | 95% CI on the paired gap |
|---|---|---|---|---|
| `backtest/loo.py` | one arm | 28.3pp | 11.5pp | [+5.7, +28.5] |
| `backtest/lomo.py` | one whole mechanism | 45.9pp | 12.3pp | [+20.1, +43.1] |

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

**1. Half the exam is unreachable.** Six of twelve quantified arms —
natalizumab, fingolimod, alemtuzumab, ponesimod, daclizumab, cladribine — fit at
the grid edge, meaning no potency produces a benefit at all. That is the entire
depleting / sequestering / trafficking class, which the ported QSP model
collapses onto a single loss term it cannot turn downward. Oracle error there is
51.9pp against a 12.4pp null. **No potency data of any quality touches this**;
it is a direction failure, not a calibration one.

**2. On the reachable half, the independent channel misses.** The oracle fits
those six arms to 0.5pp — which one free parameter per target does by
arithmetic, and establishes only that the response curve passes through each
trial's answer somewhere. Scoring the same arms at the MRI-fitted potencies from
`backtest/potency.py`, which never see a relapse number, gives **21.4pp against
a 10.6pp null.** It loses, and it is biased high every single time — median
1.55x the oracle potency, range 1.08x to 2.46x.

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

Figures on this page were measured 2026-09-20 at `519bfe4`. If a number here
disagrees with what the code prints, the code is right and this page is stale.
