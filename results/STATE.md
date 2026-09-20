# STATE OF THE BUILD

Generated 2026-09-20 23:22 UTC at `87393d3 (dirty)` by `python scripts/state_of_build.py`.
**Every number here was produced by the command beside it, in this run.**
Nothing on this page is transcribed, and nothing on it is interpreted —
for what the numbers mean, read BUILD_PLAN.md §8.4, where a human signs it.

| gate | command | result |
|---|---|---|
| clinical gate (ABM path) | `python -m backtest.clinical` | direction 10/17, magnitude 1/12 |
| clinical gate (grounded stack) | `python -m backtest.clinical_velez` | direction 5/16, magnitude 2/11 |
| leave-one-ARM-out | `python -m backtest.loo` | MAE 28.3pp, null 11.5pp |
| leave-one-MECHANISM-out | `python -m backtest.lomo` | MAE 45.9pp, null 12.3pp |
| LOMO under the capacity extension | `python -m backtest.lomo_capacity` | MAE 45.6pp, null 12.3pp |
| per-drug potency (MRI channel) | `python -m backtest.potency` | fitted 11, out of range 5 |
| harm channel | `python -m bricks.harm_channel` | enlarged ranking 7, p 0.048. |

## How to read this

- A gate with a null **loses to its null** wherever MAE > null. That has
  been true of every out-of-sample scorer in this repo since they were built.
- **MISSING** means the gate ran but its output no longer matches the pattern
  this script greps for. That is a script bug, not a result — fix the pattern.
- This page is regenerated, never edited. If a number here disagrees with a
  number written in a doc, the doc is wrong.

Nothing in this repo is evidence about multiple sclerosis.
