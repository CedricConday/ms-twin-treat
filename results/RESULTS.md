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
| identity null (predict no change) | 0.0000 | captures nothing |
| affine transfer (v0) | 0.8204 | below the bar |
| global-mean-shift null (**the bar**) | 0.8498 | interferon signature is largely shared — a *hard* bar |
| scGPT-blood embeddings | 0.8696 | beats the bar |
| **control-similarity transfer** | **0.8732** | **beats the bar — the winner, LOCTO** |
| noise ceiling (perfect model) | 0.9075 | the attainable limit |

- **We beat the null — 0.87 > 0.85, honestly.** A control-state-similarity-weighted transfer, leave-one-cell-type-out. Negative-control checked: scrambling which similarity pairs with which delta collapses it to 0.78 (below baseline), so the gain is the similarity, **not a leak**.
- **scGPT was run for real and did NOT earn the win.** Its embeddings clear the bar (0.8696) but do **not** beat plain Pearson correlation of the control profiles (0.8732) — paired over 8 folds the difference is undetectable (Wilcoxon p=0.55). *"scGPT beats the null"* is literally true but misleading; *"scGPT gave us the win"* is **false**. scGPT's **encoder** validated (94% leave-one-out cell-type accuracy, recovers immune lineage); its MLM **decoder** did not, so only the encoder + gene-embedding table were used.
- **Honest limits, stated not buried:** the bar is *leaky in its own favor* (it averages in the held-out cell type's true delta) — the fair leave-one-out null is 0.8232, and we beat that too. **Megakaryocytes** (63/69 cells, split-half reliability 0.06) is an unscoreable noise fold that drags every model down; on the 7 measurable cell types the picture is cleaner: bar 0.903, model 0.929, ceiling 0.990 — the model closes ~30% of the gap to what's attainable. With n=8 the margin is real but **not statistically conclusive** (p≈0.055); the model wins 7/8 cell types and loses Dendritic cells (0.915 vs 0.932).

## 3. Recovered biology (independent sanity signals, not fitted)
- **GRN brick (B3)** top co-expression edges from raw Kang data: **IFIT1/IFIT3–ISG15** (the canonical interferon-stimulated-gene module IFN-β induces) and **HLA-DRB1–HLA-DRA** (MHC-II). The pipeline rediscovered the interferon response without being told to.
- The Megakaryocyte outlier (§2) is a second real signal the harness surfaced on its own.

## 4. The pipeline runs end to end (10 bricks)
`python spine/run_demo.py` runs a virtual cohort through intervention → cell → GRN → QSP → population → barrier → clinical readout, and labels itself: *"BUILT, NOT VALIDATED — the scales compose; that is the claim. Nothing here is evidence about MS."*
- Toy bricks are directionally coherent: QSP untreated myelin collapses to 0.09, treated preserved 0.78; ABM untreated damage 0.78, treated 0.27.
- **Barrier is consequential (B8 gating):** same drug, same sim — a CNS-required therapy that barely crosses (effective 0.08) reads ~30 lesions vs ~12 for a peripheral one (effective 1.0). Ignore delivery and the barrier brick would be decoration; it isn't.
- **The wedge (B9):** first open Python virtual-population method pointed at a neuroimmune model, now **both halves**: (1) plausible-patient generation — LHS sampling + a plausibility filter that actually **rejects** the implausible region (acceptance ~0.89); (2) **MAPEL-style prevalence weighting** — the raw plausible set is skewed severe (72%), and weighting corrects it to a target prevalence (50/35/15 mild/moderate/severe), which is exactly the correction an unweighted plausible set lacks. Target illustrative, `validated=False`; full MAPEL optimizes multiple axes, this matches one.

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

## 5.1 The clinical backtest gate — what "viable" means
`python -m backtest.clinical` scores the stack against **known trial outcomes, in both directions** — the gate that must go green before any prediction is trusted. Anchors are real and cited via PubMed:

| arm | sim Δrelapse | known outcome | source | direction |
|---|---|---|---|---|
| untreated | +0% | neutral (control) | — | ✅ |
| IFN-β | −28% | **−27–33%** relapse reduction | PRISMS, *Lancet* 1998, PMID 9820297 | ✅ |
| glatiramer | −24% | ~**−29%** relapse reduction | Copolymer 1 / Johnson 1995, *Neurology*, PMID 11902590 | ✅ |
| APL CGP77116 | +0% | **HARMED** (halted; exacerbations) | Bielekova, *Nat Med* 2000, PMID 11017150, [doi](https://doi.org/10.1038/80516) | ❌ |

**DIRECTION gate: 3/4. VIABLE? Not yet.** The stack reproduces the two therapies that *worked* (direction and magnitude within 15pp) but **cannot reproduce the one that harmed** — the toy QSP/ABM only express benefit. That single red cell is the definition of the work ahead: a brick that can express immune exacerbation turns the gate green. The magnitude match on the successes is *directionally* honest but not yet *earned* (placeholder doses, not mechanism). *(via PubMed)*

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
