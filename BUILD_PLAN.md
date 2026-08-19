# BUILD_PLAN — ms-twin-treat (cold-boot resume doc)
*Written 2026-08-19, pre session-reset. If you are a fresh boot: read this top-to-bottom, run the COLD-BOOT checks, then execute BUILD ORDER. Goal = the 8 remaining bricks BUILT (loading + emitting output through the spine), not validated. ~4h wall-clock, install/compile-bound.*

---

## 0. STATUS (what's already real)
- **Repo:** `/home/ubuntu/repos/ms-twin-treat` — local git, 2 commits, author `Cedric Conday <cedric@condaydigital.com>`. **NOT pushed** (push held for 1300, Cedric's go only).
- **Built + validated:** brick #0 backtest harness, brick #1 Kang data loader. **2 of 10.**
- **Real result on real data:** Kang 2018 IFN-β PBMCs, 8 cell types. identity null = 0.00, global-mean-shift null = **0.85 aggregate delta_pearson** (the bar). Megakaryocytes = 0.48 (the honest outlier). No model beats the bar yet.
- **Brain repo (research):** `/home/ubuntu/repos/ms-twin/` — `docs/RESEARCH_FINDINGS.md` (all bricks, licenses, sources), `PROJECT_MS_TWIN.md`, `docs/REVIEW_CHARGE.md`.

## 1. COLD-BOOT CHECKS (run first, ~1 min)
```bash
cd /home/ubuntu/repos/ms-twin-treat
. .venv/bin/activate
python -m backtest.selftest      # must print: HARNESS VALIDATED
python -m backtest.run_kang      # must print aggregate mean-shift ~0.85
git log --oneline                # confirm 2 commits, author = condaydigital
```
If those pass, the foundation is intact — start building. If selftest FAILS, fix the ruler before anything else.

## 2. ENV FACTS / GOTCHAS (these cost real time if rediscovered)
- **CPU-only. No GPU.** No nvidia-smi. Plan every brick for CPU.
- **venv:** `.venv/` (gitignored). Pinned: numpy==2.2.6, scipy==1.14.1, pandas==2.2.3, scikit-learn==1.5.2, anndata==0.11.3, scanpy==1.10.4.
- **NUMPY ABI TRAP:** system Python's pandas is built against an older numpy and crashes on import. NEVER use system python; always the venv. Keep numpy < 2.3 or the ABI breaks again.
- **figshare downloads:** use host `ndownloader.figshare.com/files/<id>` (302→S3). `figshare.com/ndownloader/...` returns 0 bytes.
- **Datasets never committed:** `data/cache/` + `*.h5ad`/`*.parquet` are gitignored. Keep it that way.

## 3. THE BRICK CONTRACT (so bricks compose through the spine)
Two interfaces. Keep them tiny.

**(a) Harness-facing — for any brick that PREDICTS a scored perturbation** (`backtest/harness.py`):
```python
class PerturbationModel(Protocol):
    name: str
    def predict(self, control_mean: np.ndarray, cell_type: str) -> np.ndarray: ...
```

**(b) Spine-facing — for pipeline flow** (build in `spine/pipeline.py`):
```python
# MultiScaleState = a plain dict passed stage->stage. Each stage reads keys it
# needs, writes keys it produces. A brick that can't compute a key yet writes a
# clearly-labelled STANDIN value. The spine just runs stages in order.
class Stage(Protocol):
    name: str
    def run(self, state: dict) -> dict: ...   # returns the mutated/extended state
```
Rule: a brick is "BUILT" when its `run()` executes without error and writes its output key(s) to the state — even if the value is a documented stand-in. "VALIDATED" is a separate, later bar (backtest reproduces known biology).

## 4. BUILD ORDER (spine first, then bricks; parallelize installs)
Kick off heavy installs in the background FIRST (torch, pyscenic) so they finish while you build the light bricks.

- **SPINE v0** (`spine/pipeline.py`) — ~20 min. The Stage protocol + MultiScaleState + a runner that executes a list of stages and prints the state's keys at each step. This is the skeleton every brick plugs into. Do this first so bricks have a socket.

- **#2 Cell (scGPT)** — PRIMARY: scGPT small checkpoint on CPU. STAND-IN (do this first, it doubles as the first real contender vs the 0.85 bar): a **ridge cross-cell-type transfer model** in `bricks/cell_transfer.py` — hold out one cell type, predict its IFN-β delta from the others via regression on the control profile. Implements `PerturbationModel`. ~20 min. Then attempt real scGPT as time allows (torch CPU install in background).

- **#3 GRN (SCENIC/arboreto)** — PRIMARY: `pip install pyscenic` + arboreto GRNBoost2 on the Kang matrix → a gene-regulatory edge list written to state. Install is finicky; if it fights, STAND-IN: sklearn mutual-information graph over top-variable genes. ~30–45 min.

- **#4 QSP/ODE (Tellurium)** — `pip install tellurium` (libroadrunner wheel, installs clean on CPU). A toy 3-species inflammation ODE (immune→cytokine→damage) as an SBML/antimony model, integrated over time, writes a trajectory to state. ~20 min. (No open MS QSP exists — this is a labelled toy, not the real model.)

- **#5 Population ABM** — PRIMARY light stand-in (skip compiling PhysiCell now): **Mesa** (`pip install mesa`) — a tiny immune-agent vs myelin-agent grid model, writes a damage-over-time readout. ~25 min. Real PhysiCell + Weatherley port is a later, heavier job — note it, don't block on it.

- **#6 Barrier/PBPK** — STAND-IN: reimplement a minimal 2–3 compartment blood↔CSF↔CNS ODE in scipy (don't install R for Verscheijden now) → a CNS-exposure fraction written to state. ~30 min.

- **#7 Intervention model** — a parameter object that modifies the QSP/ABM inputs (dose, target) — e.g. "IFN-β on" flips the perturbation. ~15 min once #4/#5 exist.

- **#8 Clinical readout** — STAND-IN mapping function: sim damage/inflammation trajectory → a crude proxy for relapse-rate / lesion-count. Clearly labelled as an unvalidated proxy (the real micro→clinical map is open research). ~15 min.

- **#9 VPop engine (THE WEDGE)** — first Python implementation: Latin-hypercube / rejection sampling of parameter sets against plausibility bounds (Allen–Rieger flavor; MAPEL prevalence-weighting as a TODO). Generates a virtual cohort the spine can run. ~30–45 min. This is the novel piece — build a clean v1.

- **WIRE IT:** one `spine/run_demo.py` that assembles [intervention → cell → GRN → QSP → ABM → barrier → readout] over a small VPop and prints the state flowing end to end. THIS is the "built" milestone.

## 5. GUARDRAILS (non-negotiable)
- Author identity ALWAYS `cedric@condaydigital.com` (or the noreply). **NEVER** `a personal email address`.
- **Push stays held until Cedric says go (target 1300).** Local commits fine; no remote, no push.
- No datasets/weights in git. No secrets.
- **Built ≠ validated. Do not let the demo's "it runs" become "it works."** Every stand-in is labelled STANDIN in code + output. Under-promise to the MS community, always.
- Commit after each brick lands (local), so the reset/pause never costs more than one brick.

## 6. FOR THE DECK BUILDER (2nd instance — read this, don't reinvent)
**The honest narrative (use verbatim framing):**
- One engineer, from scratch, in one session: a *validated backtest harness* + a *real MS-relevant dataset* running through it, then the full multi-scale skeleton wired.
- **The discipline is the story:** the FIRST thing built was the ruler (backtest harness), not the simulation. Trust no prediction until it replays known history. That's what separates this from vaporware.
- **The real numbers** (don't inflate): identity null 0.00, mean-shift bar 0.85 on real IFN-β immune-cell data; Megakaryocytes 0.48 = the harness auto-finding where the shared-response assumption breaks.
- **The wedge:** first open Python virtual-population method ever pointed at a neuroimmune model (no impl exists on GitHub — see `ms-twin/docs/RESEARCH_FINDINGS.md`).
- **What NOT to claim:** not a cure, not a trial replacement, not a working human simulation. The bricks are BUILT (skeleton with a pulse), not validated. In-silico de-risks; it does not replace a trial. Say that on a slide — it's a strength, not a hedge.
- **Sources for facts/citations:** `ms-twin/docs/RESEARCH_FINDINGS.md` (bricks, licenses, DOIs), the two git commits, this file's §0 for numbers.

## 7. DEFINITION OF DONE (the 4h "built" milestone)
`python spine/run_demo.py` runs a small virtual cohort through all stages end-to-end without error, printing the multi-scale state accumulating keys at each brick — with every unvalidated component labelled STANDIN. Harness selftest + Kang run still green. Committed locally, brick by brick. Push still held.
