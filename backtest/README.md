# backtest/

The validation harness — the edge of the whole project. Trust nothing until it
replays known history. Scores whether the sim reproduces a *known* published
outcome before any prediction is believed. Modeled on perturbation-prediction
benchmarks (**scArchon**-style; see `../../ms-twin/docs/RESEARCH_FINDINGS.md`,
which names it as the harness to point at MS). This is what makes a claim
credible instead of a demo that looks alive but predicts nothing.

## What runs here, easiest test first

| module | question it answers | run |
|---|---|---|
| `selftest.py` | does the ruler measure? synthetic data with a known answer | `python -m backtest.selftest` |
| `clinical.py` | does the sim get the DIRECTION right on 14 cited trial arms? | `python -m backtest.clinical` |
| `loo.py` | fit on N-1 arms, predict the Nth — **calibration** | `python -m backtest.loo` |
| `lomo.py` | hold out a whole MECHANISM, predict it blind — **generalisation** | `python -m backtest.lomo` |

## Why there are two hold-out tests, and which one matters for what

`loo.py` holds out one **arm**. The held-out drug's mechanism class is still in
the training set, so the fit has seen something that works the same way. That
grades calibration, and it is the right test for *ranking known drugs*.

`lomo.py` holds out a whole **mechanism group** — every arm that moves the same
model parameters. Nothing in training acts where the held-out drug acts. That is
the situation a screened candidate is actually in, so it is the test that grades
*generating new candidates*, and it is strictly harder: for any multi-arm group
its training set is smaller than the LOO's would be.

They also run on different stacks, deliberately:

- `loo.py` uses the ABM path the clinical gate uses (`response_curve.py`).
- `lomo.py` runs on the screening stack — `bricks/qsp_velez.py` (grounded MS QSP
  model) → `bricks/profiles.py` (each arm's intervention points, from
  pharmacology) → `bricks/sormani.py` (published trial-level map from a lesion
  ratio to a relapse ratio). Running it is also a check that those three compose.

## The rule both obey

**Every benchmark reports a null beside it.** Predict-the-mean: a model that
ignores which drug it is asked about and answers with the training arms' average.
A result that does not beat that null carries no drug-specific information,
whatever mechanism story sits behind it. Neither test currently beats its null,
and both say so in their own output rather than in a footnote.

## Cohort size is not a taste decision

The QSP model's untreated damage is extraordinarily right-tailed. Measured over
200 runs at a 730-day horizon (`scripts/measure_qsp_variance.py --n 200`):
median 2.14, **mean 74.9**, SD 612, range 0.089 to 8118 — a ~90,000-fold spread,
with the top decile of runs holding 96% of all damage. The standard error at n=4
is four times the mean, and an early version of `lomo.py` reported a headline
made entirely of that noise.

So `lomo.py` uses the **median over 128 seeds**. Bootstrapped CV of the median:
628% at n=4, 41.8% at n=16, 21.2% at n=48, **13.7% at n=128**. Read that ~14% as
the resolution limit on every number this directory produces — a smaller
difference between two arms is not a result.

One intuition that did NOT survive measurement: sharing a seed between the
treated and untreated arms looks like it should cancel the infection history,
but it narrows the ratio's IQR by only ~1.1x. Pairing is kept because it is free.
It is not what makes the comparison valid; the cohort size is.

`results/mechanism_curve.json` and `results/response_curve.json` are the cached
response tables. Both are committed so these tests are reproducible without a
multi-minute rebuild; delete one to force it.
