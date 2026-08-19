# ms-twin-treat — Results (v0)
*Verified numbers from a single build session. Every figure below is reproducible from this repo. Read §6 before quoting anything: most of this is toy, and the honesty is the point.*

## The one honest headline
From scratch, one session: a **validated backtest harness** + a **full 10-brick multi-scale MS-intervention pipeline** that runs end to end, **recovers real immune biology** in two independent places, and is **scrupulous about what is validated vs. toy**. It de-risks candidate interventions in silico; it does **not** replace a trial, and nothing here is evidence about multiple sclerosis.

## 1. The ruler is validated (backtest harness)
The first thing built was the scorer, not the simulation. It passes its own gate on synthetic data with known structure:
- oracle (sees truth) → **1.00**, global-mean-shift null → 0.90, identity null → **0.00**. Ordering holds ⇒ the ruler discriminates.

## 2. Real-data backtest — Kang 2018 IFN-β PBMCs (8 immune cell types)
| model | aggregate delta_pearson | reading |
|---|---|---|
| identity null (predict no change) | **0.0000** | IFN-β changes a lot; predict-nothing captures nothing |
| global-mean-shift null | **0.8498** | the interferon signature is largely shared — a *hard* bar |
| **cell-transfer model (B2)** | **0.8204** | **honest: does NOT beat the bar yet** |

- Per-cell-type mean-shift ranges 0.87–0.93 **except Megakaryocytes = 0.476** — the harness auto-finds the one cell type where "everyone responds the same" breaks down.
- **Honest state:** our cell model (affine cross-cell-type transfer, leave-one-cell-type-out) scores *below* the null. A model that can't beat "everyone responds the same" hasn't earned "cell-type specific." Beating 0.85 is the job handed to the scGPT swing (2nd instance).

## 3. Recovered biology (independent sanity signals, not fitted)
- **GRN brick (B3)** top co-expression edges from raw Kang data: **IFIT1/IFIT3–ISG15** (the canonical interferon-stimulated-gene module IFN-β induces) and **HLA-DRB1–HLA-DRA** (MHC-II). The pipeline rediscovered the interferon response without being told to.
- The Megakaryocyte outlier (§2) is a second real signal the harness surfaced on its own.

## 4. The pipeline runs end to end (10 bricks)
`python spine/run_demo.py` runs a virtual cohort through intervention → cell → GRN → QSP → population → barrier → clinical readout, and labels itself: *"BUILT, NOT VALIDATED — the scales compose; that is the claim. Nothing here is evidence about MS."*
- Toy bricks are directionally coherent: QSP untreated myelin collapses to 0.09, treated preserved 0.78; ABM untreated damage 0.78, treated 0.27.
- **Barrier is consequential (B8 gating):** same drug, same sim — a CNS-required therapy that barely crosses (effective 0.08) reads ~30 lesions vs ~12 for a peripheral one (effective 1.0). Ignore delivery and the barrier brick would be decoration; it isn't.
- **The wedge (B9):** first open Python plausible-patient generator pointed at a neuroimmune model — LHS sampling + a plausibility filter that actually **rejects** the implausible region (acceptance ~0.89). MAPEL prevalence-weighting is the flagged next step.

## 5. Arm experiment — the shape of the eventual backtest
One plausible virtual population (12 patients), run through every real trial arm (`python -m results.experiment`):

| arm | treat | mean lesion_proxy | Δ vs untreated | real-world outcome |
|---|---|---|---|---|
| untreated | 0.00 | 32.7 | +0.0 | control |
| IFN-β | 0.50 | 23.7 | **−9.0** | approved / worked |
| glatiramer acetate | 0.45 | 24.8 | **−7.9** | approved / worked |
| APL CGP77116 | 0.00* | 32.7 | +0.0 | **HARMED — Phase II halted (Bielekova, *Nat Med* 2000)** |

- The pipeline **separates the therapies that worked from untreated.**
- It **cannot yet flag the harmful arm** — APL's effect is a 0.0 placeholder, so it reads like untreated. *That is the exact gap the backtest exists to close, and we are not hiding it.* A model earns trust by catching the one that harmed people; ours can't, yet.

## 6. What is real vs. toy (read this before quoting anything)
| Real / defensible | Toy / illustrative / unvalidated |
|---|---|
| The backtest harness + its gate | Every disease/PK/ABM parameter value |
| The Kang IFN-β numbers (0.00 / 0.85 / 0.82 / Megakaryocyte 0.48) | The QSP, ABM, barrier, readout *models* |
| The recovered interferon module (GRN) | The named-drug → `treat` mapping |
| The plausible-patient *method* (LHS + rejection) | The plausibility *model* it filters against |
| That the scales compose end to end | Any clinical proxy (lesion/relapse numbers) |

**In-silico evidence de-risks candidates and augments trials. It does not replace a pivotal trial or a control arm — anywhere, for any therapy. Nothing in this repo is evidence about multiple sclerosis.**

---
*Build: two instances in parallel (vps + cbk), 10 bricks + spine, one session. Author: Cedric Conday. Not pushed until go.*
