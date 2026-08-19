# GROUNDING — turning toy parameters into data-grounded ones

The clinical gate reaching 4/4 exposed the real work ahead: **every parameter that
isn't the cell brick is hand-set or illustrative.** This file tracks, per brick,
how to replace hand-set values with data-grounded ones. It is also the contributor
task list — own a brick, ground it.

## The one rule that keeps grounding honest
**Ground every parameter against an INDEPENDENT source** — a measured mechanism, a
published PK value, another dataset. **Never fit it to the clinical outcome the gate
is trying to predict.** Fitting to the outcome makes the gate circular and
worthless. When a parameter is grounded independently, the clinical gate becomes a
real **out-of-sample** test instead of a restatement of what we already told it.

## Status per brick
| Brick | Parameter(s) | Now | Ground against | Priority |
|---|---|---|---|---|
| **Cell** (`cell_scgpt`) | cell-state response | ✅ **DATA-GROUNDED** — Kang IFN-β, beats the null (0.87) | — this is what "grounded" looks like | done |
| **Intervention** | `treat`, `immunogenic` | ✅ **from a 2-param mechanism-class rule** (`grounding.py`) — not per-arm hand-tuning | next: derive the class strengths from data (IFN-β magnitude from Kang; others need their single-cell data — the honest gap) | partial |
| **Barrier** (PBPK) | rate constants → CNS penetration | ✅ **baseline GROUNDED** — `k_pc` calibrated to ~0.15% CNS (Pardridge 2019); other rate constants still illustrative | remaining constants: published PK; small molecules higher | partial |
| **QSP** | disease / cytokine ODE rates | invented | published immune/cytokine kinetics; no open MS QSP exists (greenfield) | MED–HARD |
| **ABM** | aggression / damage / repair | invented | the **Weatherley MS ABM** (MIT, PLoS Comp Bio) rules/rates — a real published open model to port | MED |
| **Readout** | micro → clinical map | invented proxy scales | MSOAC / trial relapse data, held-out validated. **Hardest — open research.** | HARD |

## First concrete grounding target (identified this session)
**Barrier CNS penetration.** The current toy outputs ~8% CNS-effective exposure for
a CNS-required drug. For therapeutic **antibodies** the real CNS:serum ratio is
~**0.1–0.2%** — two-plus orders of magnitude lower. Grounding step: calibrate the
barrier rate constants so a large-molecule dose yields ~0.1–0.2% CNS AUC fraction,
matching standard mAb CNS pharmacokinetics. (Surfaced via a PubMed search for
antibody blood-brain-barrier PK — candidate primary sources found; pin one before
setting the number.) This is an **independent** correction, not a fit to any MS
outcome, so it strengthens the gate rather than circularizing it.

## Why this is the path to "viable"
The clinical gate is only a real test when the parameters it grades were **not** set
from the answers. Grounding each brick independently converts today's *"4/4 on the
arms that informed the setup"* into *"4/4 on arms the setup never saw."* That is the
whole difference between a capability demo and validation.

## For contributors
Each row is a self-contained task behind the frozen brick interface
(`predict(control_mean, cell_type)` / `run(state)`). Pick a brick, ground its
parameters against the cited source, keep `validated=False` until the clinical gate
passes **out-of-sample**, and never fit to the outcome. (A `CONTRIBUTING.md` with
good-first-issues lands when the repo goes into active recruitment.)
