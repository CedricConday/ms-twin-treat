# MS-Twin — Research Findings (2026-08-19)
*Synthesized from 4 parallel research agents + direct HF-MCP / GitHub pulls. Every access claim verified. This is the "found all the pieces" deliverable — ready for adversarial review (Fable + K3).*

## TL;DR for the reviewers
A validated, open MS digital twin **does not exist** — it's mostly greenfield, but with real prior art to stand on. Every *brick* we need exists as maintained OSS. The single **first-of-kind, defensible, regulator-wanted differentiator** nobody has claimed: **rigorous virtual-population methods applied to a neuroimmune/MS model** (no Python implementation of the plausible-patient/prevalence method exists on GitHub *at all*). That is the thesis worth building around. **Honest ceiling (non-negotiable): in-silico evidence cannot replace a pivotal efficacy trial or control arm — anywhere, for any product.** This de-risks and augments; it does not replace. Everything below is scoped inside that truth.

## The bricks — verified best OSS option per scale

| Brick | Winner | License | Why | Compute |
|---|---|---|---|---|
| **B1 Cell** | **scGPT** (bowang-lab) | **MIT** ✅ commercial-safe | only permissive one; ships **brain checkpoint** (microglia/oligo/neurons) + **blood** (T/B) = exact MS cell types; 53M params | laptop/1 GPU |
| B1 alt (drug/cytokine) | **Arc State** (`ST-SE-Parse`) | NC | only model built for perturbation→response; PBMC+90-cytokine head = MS T/B arm; **no CNS cells** | 1 GPU |
| ✗ AIDO Cell (the 08-18 hype) | — | **CLOSED, waitlist** | v1.0 simulates only K-562 + Hep-G2 (cancer lines) — **MS-irrelevant.** Write off for now. | n/a |
| **B_reg Regulatory** (new scale) | **arboreto/SCENIC**, `inferelator`, `Beeline` (bench) | OSS | GRN inference = the "regulatory network" scale; models the signaling driving the autoimmune attack | 1 GPU |
| **B3 Population** | **PhysiCell** (`development` branch — releases stalled, branch active) + Studio + **PhysiBoSS/MaBoSS** | BSD | built-in **BioFVM diffusion solver** = native B3↔B4 coupling; per-agent SBML; XML not C++. Pure-Python alt: **CompuCell3D**. *(Dead — do NOT build on: Simmune, C-ImmSim.)* | CPU/multicore |
| B3 seed | **`Georgia-Weatherley/MS_ABM_Weatherley`** | **MIT** | the ONE real open MS ABM (PLoS Comp Bio, Jan 2026): T-cells×BBB, macrophages, myelin agents, oligodendrocyte repair — **fork for validated rules** (MATLAB) | — |
| **B4 Barrier** | **PK-Sim + MoBi / OSP** (whole-body PK, incl. 5-compartment brain/CSF; very actively maintained) + **Verscheijden 2019** 14-compartment brain PBPK (**R code in the supplement — most immediately runnable open CNS model**) + **ADMET-AI** (`pip`, MIT, BBB gate) | GPLv2 / MIT | no OSS covers nanoparticle/antibody CNS entry — build the natalizumab (anti-VLA-4 diapedesis) arm from published transcytosis models | modest |
| **B_mech** (QSP/ODE) | **Tellurium / COPASI** (SBML in agents), **Julia SciML**, **NVIDIA nvQSP** (GPU) | OSS | mechanistic disease/drug ODE layer, embeddable per-agent | varies |
| **Backtest harness** | **`scArchon`** (perturbation-prediction benchmark pipeline) | OSS | already scores whether perturbation models predict correctly — point it at MS | — |

## Data (B2) — RICH; processed layers open, only raw FASTQ gated
- **CNS:** Macnair 2025 (632k nuclei, 83 people, Zenodo `10.5281/zenodo.8338963`, CC-BY, patient stratification) · **Lerma-Martin GSE279183 (analysis-ready h5ad, 1.7GB — best drop-in)** · Schirmer 2019 · Meijer/Agirre (GWAS×chromatin) · Feng/Groh GSE301908 (CD8+ niches, foamy microglia) · Alsema (spatial).
- **CSF/blood:** Ban 2024 (CSF eQTL→CD8 T), Ashida GSE286068 (CITE-seq).
- **Perturbation backtest anchor:** **Kang 2018 GSE96583 — IFN-β on human PBMCs. IFN-β is a first-line MS DMT → the prime open "real-MS-drug-on-immune-cells" v0 backtest.** Plus Cui 2024 (86 cytokines), Human Cytokine Dictionary (9.7M PBMCs).
- **Clinical-layer backtest anchor:** **C-Path MSOAC placebo database** — pooled, curated MS trial control-arm data (relapse/EDSS/MRI). This is where the monoclonals get validated, and it's the exact data regulators recognize.
- **⚠️ GAP:** no open single-cell drug-response atlas for the major DMTs (ocrelizumab/fingolimod = 0). Cellular backtest = IFN-β/cytokine only; monoclonals validate at the **clinical readout** layer (relapse/MRI/EDSS from published trials + MSOAC), not cellular.

## Prior art — greenfield, but you stand on real work
- **No usable open MS twin exists.** BioModels has ZERO MS/demyelination/myelin/oligo models.
- No commercial QSP platform has an MS product (Certara/SLP consolidation; SLP going private; UISS-MS closed). Deepest autoimmune lib (SLP) = RA/lupus/IBD, **zero neuro**.
- **Regulatory reality (sharp):** EMA qualified PROCOVA (Unlearn) — the *statistics only, not the AI*; FDA declined it; real benefit 10–15%; Unlearn's MS line looks abandoned. Only EMA-qualified sim platform = Simcyp PBPK (DDI/CYP only). ICH M15 model-credibility framework hit Step 4 (Jan 2026) — a documentation standard, not a trial-replacement license.

## THE DIFFERENTIATOR (first-of-kind, unclaimed, regulator-wanted)
Rigorous **virtual-population** methods have **never been applied to any neuroinflammation or MS model**, and **no Python implementation of the plausible-patient/prevalence-selection method exists on GitHub at all.**
- The methods: **Allen–Rieger–Musante plausible-patient generation** (2016, `10.1002/psp4.12063`) + **Schmidt's MAPEL** prevalence-density selection. Reference implementation = **`BMSQSP/QSPToolbox`** — MATLAB + SimBiology (proprietary stack). That's the whole state of the art, and it's locked in a language/toolbox nobody in open science can run freely.
- EMA/ISoP explicitly say there are "no guidelines yet" for VPop generation.
- ⚠️ **Someone is building adjacent — and closing it:** Russo et al., *NX210c QSP Model of RRMS*, Int J Mol Sci 2026;27(3):1349 — a BBB module (claudin-5/TEER) + RRMS virtual populations in MS TreatSim v2.0. **Code and data explicitly withheld as proprietary.** The gap is being probed; the *open* lane is still empty. Move.

→ **A Python port of the prevalence/plausible-patient method + modern sampling (DREAM(ZS) / simulation-based inference) + the Weatherley MS ABM (or an open CNS/BBB model) = genuinely first-of-kind, defensible, and exactly what regulators are asking for.** This is the wedge, and the wedge has a clock on it.

## Recommended v0 (the hello-world, honestly scoped)
1. Fork **`MS_ABM_Weatherley`** rules → port to **PhysiCell** (BioFVM couples drug field to agents).
2. **scGPT** brain+blood checkpoints as the cell-state brick.
3. Backtest anchor: **Kang GSE96583 IFN-β response** — does the loop reproduce the known IFN-β effect on immune cells?
4. Wrap in **scArchon**-style scoring.
5. *Stretch, the real novelty:* a Python plausible-patient/prevalence VPop layer over it.
- **Build yourself (no framework gives these):** BBB diapedesis (natalizumab), microglial M1/M2 dynamics, remyelination.

## Honest risks / ceiling
- In-silico ≠ trial replacement. Ever. This de-risks candidates and augments control arms. Frame everything this way — especially to the MS community.
- Scale-bridging (the spine) is unsolved research, not glue.
- Monoclonal-DMT cellular data doesn't exist; that arm is clinical-layer only.
- Licenses: scGPT MIT is the commercial-clean core; GB.Cell/Arc-State/CellFM are non-commercial — matters for Conday Digital productization.

---
*4 research agents, ~490k tokens, verified. Ready for Fable + K3 adversarial review.*
