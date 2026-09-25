# MS-Twin — a validated in-silico testbed for MS interventions

**One honest line:** a multi-scale simulation pipeline that can be *backtested against known outcomes*, so a candidate MS therapy can be tested in silico before it burns a trial — starting narrow, built in the open, extended brick by brick.

**What this is NOT (read first, every time):** this is not "a cure," not "curing MS," not a promise to patients. It is a *tool* — a testbed that tells the good candidates from the doomed ones. The credibility of this whole project is that every claim is backed by a backtest that reproduces known history. **We demonstrate; we do not assert.** Overclaiming to a vulnerable patient community is the one thing that ends this. (See §7.)

---

## 1. The thesis (why this exists)
Nanomedicine/immunotherapy for MS stalls at Phase 2 because the field **can't predict human outcomes from animal models** — no accurate simulation of the human system. Master blocker = a prediction problem. Prediction is buildable. The value isn't proving cures — it's *separating real candidates from false ones before $200M is spent.* Source map: `docs/why-it-stalls.md` (to write from the two research passes already done).

## 2. Architecture — the spine + the bricks
**The spine (what WE build):** a multi-scale orchestration layer. State flows
`molecular → cell → cell-population → tissue/barrier → clinical readout`.
Scale-bridging is the open research problem. It is shaped like an AnyLogic multi-scale sim wired to a data pipeline — i.e. our actual skillset. The bricks are other people's models; the spine and the validation harness are ours.

> ⚠️ **This table is the 2026-08-19 morning draft and is SUPERSEDED by `docs/RESEARCH_FINDINGS.md` (same day, post-survey).** Where they disagree, the findings doc wins — it is the verified one. Kept here for provenance, not as the current plan. Current stack: B1 = **scGPT (MIT)**, B3 = **PhysiCell + `MS_ABM_Weatherley`**, B4 = **PK-Sim / Verscheijden**, plus two scales this table never had: **B_reg (GRN — arboreto/SCENIC)** and **B_mech (QSP/ODE — Tellurium/COPASI)**.

| # | Brick | Role | OSS status | Notes |
|---|---|---|---|---|
| B1 | Cell model | perturbation → predicted cell state | **superseded → scGPT** | ~~AIDO Cell (GenBio, preview)~~ — written off, see findings. OSS fallbacks became the plan: **scGPT (MIT), CellFM** (GitHub/HF). |
| B2 | Ground-truth data | real MS/immune cell states | open | Human Cell Atlas, Arc Virtual Cell Atlas, patient-derived MS single-cell datasets |
| B3 | Population model | agent-based immune assault on myelin | OSS | Mesa (Python) for v0 — **AnyLogic wheelhouse** |
| B4 | Barrier model | BBB / PBPK — what crosses into CNS | mixed | some open PBPK models exist |
| B5 | Intervention model | the therapy under test | build | tolerance nanoparticle etc. |
| B6 | Readout | map sim → clinical measures | build | MRI lesion load, relapse rate, EDSS |

## 3. v0 — the hello-world loop (first commit, tractable NOW, solo, laptop)
Do **not** simulate a human. Build the narrowest validated loop:
> **one cell type · one known MS intervention · one readout** →
> `B1 cell-model → tiny B3 agent model → B6 readout` →
> **backtest against ONE published result.**

- **Success:** it reproduces the known outcome → we have a *validated micro-loop + a spine to extend.*
- **Failure:** we learn exactly where the gap is → also the point.

## 4. The backtest discipline (non-negotiable — this is the edge)
Identical to a trading backtest: **trust no prediction until the sim replays known history.**
- Feed it a therapy that **failed** in trials → does the sim fail it?
- Feed it one that **worked** → does it work?
- Only after it reproduces the record does it get to predict the unknown.
This is precisely what the "it-worked-in-mice" era skipped. Backtest harness lives in `backtest/`.

## 5. Adversarial review (the free geniuses)
Before *and* during build, run the multi-agent adversarial pipeline (the proven method):
- **Fable** and **K3** (tested 2026-08-18) as independent reviewers, role-differentiated prompts, reviewing each other and the design.
- Cross-vendor review → local gate before anything ships. Cheap, ruthless, ours.

## 6. First-week recon (do these before building B1)
- [ ] **AIDO Cell access** — API-gated? weights downloadable? terms? If closed → start on scGPT/CellFM.
- [ ] **Datasets** — pull one open MS single-cell dataset + identify one *published intervention with a known outcome* to backtest against.
- [ ] **Population model** — stand up a trivial Mesa immune-vs-myelin toy.
- [ ] **Readout** — pick the single clinical measure for v0 (lesion count is likely simplest).
- [ ] **Reviewers** — draft the role prompts for Fable/K3.

## 7. Build-in-the-open — strategy AND guardrail
**Strategy (sound):** public repo, honest progress, real commits. The story is compelling *because it's real*. Traction and funding follow demonstrated work — same as the 40 merges.
**Guardrail (protects you):** the MS community is people hoping for their lives back. The way you win their trust and the funding is by **under-promising and over-delivering** — showing the backtest before the claim, never the reverse. Bravado is the fuel; the backtest is the currency. Do not let the bravado write a check the sim hasn't cashed. The moment this reads as "man promises to cure MS" instead of "man ships verifiable tools," it's over. Ship quietly; let the work make the noise.

## 8. Open questions / honest risks
- Scale-bridging is unsolved — this is research, not glue. (The risk *and* the reason it's worth doing.)
- Cell-model access may be gated (B1).
- Validation data for MS is real but messy.
- This is a program, not a weekend — but v0 is a weekend-to-months thing, and v0 is the whole point right now.

---
*Started 2026-08-19. Directed by Cedric Conday. Built in the open. Demonstrate, don't assert.*

---

## KEYSTONE FIND (2026-08-19) — the blueprint + the bricks are real, and pullable via HF MCP

**The blueprint paper Cedric surfaced:** *"How to build an AI-driven digital organism"* (Nature Medicine 2026) = peer-reviewed version of **arXiv 2412.06993**, by **Le Song, Eran Segal (Weizmann), Eric Xing (CMU)** — *the GenBio AI founders.* AIDO Cell (found 08-18) is one module of their **AIDO** system.
- Their thesis, verbatim: *"Manipulating biology in the physical world is extremely complex, expensive, and risky, and should be preceded by extensive computer-aided digital design, simulation, and validation as in other industrial fields such as civil, nuclear, and semiconductor engineering."* — i.e. Cedric's "simulate first, do it right once."
- Their architecture: *"a system of integrated multiscale foundation models, modular, connectable, and holistic... molecules → cells → individuals"* — i.e. "Legos, connect the bricks."
- **The named OPEN CHALLENGE (this is our lane):** *"system-wide harmonization through nested or hierarchical representation propagation"* to CONNECT the component FMs. GenBio built the bricks; **connecting them into a working, validated, disease-specific simulation is not solved — that's the spine, and that's ours.**

> ⚠️ **UNRESOLVED CONTRADICTION — reviewers, do not spend your teardown on this; it is already known.** This section (morning) says the GenBio bricks are pullable from HF. `docs/RESEARCH_FINDINGS.md` (afternoon, verified) says **AIDO Cell is closed/waitlist and its v1.0 simulates only K-562 and Hep-G2 cancer lines — MS-irrelevant.** These may both be true (public `GB.*` weights on HF ≠ access to the AIDO Cell service/checkpoints), but that has **not** been confirmed. Until it is: **B1 = scGPT.** The GenBio material below stands as the architectural blueprint argument, not as a resourced dependency.

**The bricks — HF-listed (org `genbio-ai`, ~30 models); actual downloadability UNVERIFIED, see flag above:**
| Scale | Model (HF id) | Params | License |
|---|---|---|---|
| Cell (B1) | `genbio-ai/GB.Cell-100M` | 99.9M | **other (VERIFY)** |
| Cell (small) | `genbio-ai/GB.Cell-10M`, `GB.Cell-3M` | 3–10M | other |
| Tissue | `genbio-ai/GB.Tissue-60M`, `GB.Tissue-3M` | spatial-transcriptomics | other |
| DNA | `genbio-ai/GB.DNA-300M`, `GB.DNA-7B` | | other |
| RNA | `genbio-ai/GB.RNA-650M`, `GB.RNA-1.6B` (+variants) | | other |
| Protein | `genbio-ai/GB.Protein-16B`, `GB.Protein-RAG-3B/16B` | | other |
| Structure | `genbio-ai/GB.Structure{Encoder,Decoder,Tokenizer}` | | other |
| Pathology | `genbio-ai/genbio-pathfm` (histopath) | | other |
- Framework code: **github.com/genbio-ai/AIDO** (verify license + runnability — agent on it).

**⚠️ LICENSE FLAG:** every model is tagged `license:other` — a custom GenBio license, NOT confirmed permissive. **Research/learning: almost certainly OK. Commercial (Conday Digital) product: read the actual license first.** This gates monetization, not exploration.

**COMPUTE:** GB.Cell-3M/10M/100M are small enough for a modest GPU (even CPU inference for the 3M/10M) — a v0 loop is compute-feasible on Cedric's box. The 7B–16B models need real GPU.

**STRATEGIC REFRAME (honest):** Do NOT try to out-build GenBio's digital organism — funded lab, top scientists. **Ride their open bricks; build the MS-specific spine + backtest harness on top.** The expensive part (foundation models) is done and partly open. Your slice = the neuroimmune/MS disease-twin + validation. Tractable. Mission-aligned. Yours.

---

## DATA + BACKTEST ANCHORS (2026-08-19, verified) — B2 is RICH, one honest gap

**Open MS single-cell data is abundant (processed layers open; only raw FASTQ gated):**
- **Macnair 2025** — 632k nuclei, 83 people (54 MS/29 ctrl), **fully open** processed counts (Zenodo `10.5281/zenodo.8338963`, CC-BY). Patient stratification.
- **Lerma-Martin GSE279183** — analysis-ready **h5ad** (1.7 GB) — best drop-in, no reprocessing.
- Schirmer 2019 (GSE118257), Meijer/Agirre GSE166179 (GWAS×oligodendroglial chromatin), Elkjaer, Feng/Groh GSE301908 (306k nuclei, CD8+ niches, foamy microglia), Alsema GSE208747 (spatial).
- CSF/blood: Ban 2024 (CSF eQTL→CD8 T cells, gated), Ashida GSE286068 (CITE-seq), Beltrán GSE127969 (twin CSF).

**Perturbation priors (for intervention/backtest at cellular level):**
- **Kang 2018 GSE96583 — IFN-β on human PBMCs. IFN-β is a first-line MS DMT → the single best open "real-MS-drug-on-immune-cells" dataset. Prime v0 backtest anchor.**
- Cui 2024 GSE202186 (86 cytokines), Human Cytokine Dictionary (9.7M PBMCs, 90 cytokines, CC-BY-NC), scPerturb (Zenodo `13350497`, 54 harmonized perturbation h5ads).

**⚠️ HONEST GAP (scope the backtest around this):** there is **NO open single-cell perturbation atlas for the major MS DMTs** (ocrelizumab, fingolimod = 0; alemtuzumab = 4 samples). So *cellular*-level backtesting is limited to **IFN-β / cytokine** response. Backtesting the big monoclonals must happen at the **clinical** readout layer (relapse rate / MRI / EDSS from published trials), not cellular. Two backtest layers, different data — know this before scoping.

**Minimal open v0 stack:** GSE279183 (h5ad) + Macnair Zenodo 8338963 (CNS) · GSE286068 (CSF) · Tabula Sapiens Immune (baseline) · **GSE96583 IFN-β as the perturbation backtest anchor.** Every processed layer needed is open.
