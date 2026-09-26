# STATE OF THE BUILD

Generated 2026-09-26 10:05 UTC at `4ede3a4 (dirty)` by `python scripts/state_of_build.py`.
**Every number here was produced by the command beside it, in this run.**
Nothing on this page is transcribed, and nothing on it is interpreted —
for what the numbers mean, read BUILD_PLAN.md §8.4, where a human signs it.

| gate | command | result |
|---|---|---|
| clinical gate (ABM path) | `python -m backtest.clinical` | direction 11/24, magnitude 1/15 |
| clinical gate (grounded stack) | `python -m backtest.clinical_velez` | direction 6/23, magnitude 3/14 |
| leave-one-ARM-out | `python -m backtest.loo` | MAE 30.0pp, null 11.3pp |
| leave-one-MECHANISM-out | `python -m backtest.lomo` | MAE 45.4pp, null 11.8pp |
| LOMO + capacity (K per fold, then fixed K=2000) | `python -m backtest.lomo_capacity` | MAE 42.2pp, null 11.8pp, fixed K=2000 48.0pp |
| per-drug potency (MRI channel) | `python -m backtest.potency` | fitted 11, out of range 5 |
| dial-level ceiling (out of sample) | `python -m scripts.dial_ceiling` | MAE 10.7pp, null 11.5pp, headroom 0.9pp |
| harm channel | `python -m bricks.harm_channel` | enlarged ranking 7, p 0.048 |
| exam v2, Velez transcription | `python -m backtest.exam_v2 --model velez` | S1 MAE 30.7pp, null 20.0pp, S2 headroom 10.9pp |
| exam v2, Jenner 2026 two-equation probe | `python -m backtest.exam_v2 --model minimal` | S1 MAE 43.8pp, null 30.6pp, S2 headroom 17.9pp |
| exam v2, Pernice 2020 port (Figure S2 readings) | `python -m backtest.exam_v2 --model pernice` | S1 MAE 35.1pp, null 19.8pp, S2 headroom 13.6pp |
| exam v2, Pernice 2020 port (memory-read variant) | `python -m backtest.exam_v2 --model pernice --variant memread` | S1 MAE 31.2pp, null 19.8pp, S2 headroom 13.6pp |
| exam C (Cochrane 2024 network), Pernice port | `python -m backtest.exam_cochrane` | CS1 MAE 10.0pp, null 8.6pp, CS2 headroom 7.8pp |

## How to read this

- A gate with a null **loses to its null** wherever MAE > null. That has
  been true of every out-of-sample scorer in this repo since they were built.
- **MISSING** means the gate ran but its output no longer matches the pattern
  this script greps for. That is a script bug, not a result — fix the pattern.
- This page is regenerated, never edited. If a number here disagrees with a
  number written in a doc, the doc is wrong.

Nothing in this repo is evidence about multiple sclerosis.
