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
- Author identity ALWAYS `cedric@condaydigital.com` (or the GitHub noreply). **NEVER** a personal email address.
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

---

# 8. RESUME — 2026-09-17: the six blockers

*Sections 0–7 above describe the 2026-08-19 build session and are **historical**. §0's
"2 commits, NOT pushed" and §5's "push stays held" are both spent — the repo is public,
CI runs on every push, and the cell model now beats the bar §0 says nothing had beaten.
Read them as a record, not as instructions. This section is the live one.*

## 8.0 State, measured 2026-09-17 (not recalled — re-run before trusting)

Every gate is green and nothing is mechanically broken:

```
python -m ruff check .        # All checks passed
python -m pytest -q           # 54 passed
python -m backtest.selftest   # HARNESS VALIDATED
python spine/run_demo.py      # exit 0, 7 stages, 8/13 keys STANDIN
gh run list                   # quality + CodeQL green on the last 3 pushes
```

What is *not* green is the science. `python -m backtest.clinical`:

| arm | sim Δrelapse | real | mean adjusted damage |
|---|---|---|---|
| untreated | 0% | 0% | 0.0981 |
| IFN-beta | **−77%** | −30% | 0.0225 |
| glatiramer acetate | **−77%** | −29% | 0.0225 |
| APL CGP77116 | +28% | harms | 0.1257 |

**DIRECTION 4/4. MAGNITUDE 0/2.** That is the whole problem, and the two suppressive
arms being byte-identical is half of it.

## 8.1 The six blockers

**Nothing here is externally blocked.** 1–5 are work. 6 has an external dependency that
two open routes already clear. Ordered cheapest-first, not most-important-first.

### (5) Readout rounding quantizes the gate — *30 minutes*
`bricks/readout.py:49` rounds `relapse_proxy` to 2dp before `backtest/clinical.py:65`
averages it and `:83` takes a ratio. At treated-arm damage ~0.0225 the per-patient proxy
takes **four distinct values** — `{0.02, 0.03, 0.04, 0.05}` — around a mean of 0.0342.
That is ~25% discretization error per patient feeding the gate's only quantitative number.
**Fix:** keep the unrounded value for scoring, round only for display.

### (1) B4 QSP is still an invented toy — *~2–3 days*
`bricks/qsp.py:1` — a 3-species ODE with hand-picked rates. It was a toy *by design*
(§4 above budgeted 20 minutes for it) because in August no open MS QSP model was found.
**That finding was too strong.** There are two published, open-access MS-specific ODE
models with equations and parameter tables in the papers:

- **Vélez de Mendizábal et al. 2011**, *BMC Syst Biol* 5:114 — effector/regulatory T-cell
  cross-regulation in RRMS; relapses **emerge from the dynamics** rather than being
  scheduled. doi:10.1186/1752-0509-5-114
- **Martinez-Pasamar et al. 2013**, *BMC Syst Biol* 7:34 — the successor, adding
  antigen-specific subpopulations and microglia. PMC3651362

**Fix:** port one, the same way `bricks/abm.py` ports Weatherley — transcribe, annotate
every rate with its source, tune nothing. Re-verified 2026-09-17 that no open MS QSP
model exists to `pip install`; Russo et al.'s NX210c RRMS model (2026) is code-withheld.

### (2) B8 readout is an invented scale — *~2–3 days*
`bricks/readout.py:26-27` — `MAX_LESIONS = 40`, `MAX_RELAPSE = 1.5`, both stated as
illustrative. The docstring calls the micro→clinical map "open research". **Also too
strong:** there is a published, fitted map.

- **Kotelnikova et al. 2017**, *PLOS Comput Biol* 13(10):e1005757 — brain damage →
  disability, fitted to clinical data and validated on a 120-patient prospective cohort.

**Catch, decide before starting:** Kotelnikova maps damage → **EDSS**, and the gate is
scored on **relapse-rate** change. Either add an EDSS endpoint to the gate (PRISMS and
Copolymer-1 both report EDSS progression, so the anchors exist) or keep ARR and calibrate
separately. Do not paper over the mismatch.

**Note the dead end first:** `MAX_RELAPSE` and `MAX_LESIONS` **cancel exactly** in the
gate, because `:83` scores a ratio against untreated. The 2026-09-17 sweep in
`scripts/sweep_maxrelapse.log` (0.8 / 1.5 / 3.0 → −74.5% / −76.7% / −77.4%) is measuring
blocker (5)'s rounding noise, nothing else. **No value of these constants can move the
magnitude.** Do not sweep them again.

### (4) One class strength for every suppressive drug — *~1 day, after (3)*
`bricks/grounding.py:67` — `SUPPRESSIVE_STRENGTH = 0.78` for IFN-β *and* glatiramer *and*
anything else in the class, which is exactly why both arms print −77%. `:71`
`IMMUNOGENIC_STRENGTH = 0.4` is not grounded at all and is exercised by one arm.
`mechanism_to_params()` **already accepts a per-drug `strength` override** — the hook is
built, nothing feeds it. **Fix:** feed it per-drug potency from an independent source
(see blocker 6), never from the arm's own relapse number.

### (3) Four arms is too few to test anything — *~2–3 days*
`backtest/clinical.py:46` — `KNOWN_OUTCOMES` has 4 entries, all of which informed the
setup. `bricks/grounding.py`'s own docstring names this: a leave-one-arm-out test "needs
more arms than we have data for". **Fix:** grow to ~10–20 arms with cited trial numbers
and independent mechanism classes, then make the gate a real **LOO**: fit class/drug
strengths on N−1 arms, predict the Nth. Candidates with published ARR reductions —
natalizumab (AFFIRM), fingolimod (FREEDOMS), ocrelizumab (OPERA I/II), teriflunomide
(TEMSO), dimethyl fumarate (DEFINE/CONFIRM), alemtuzumab (CARE-MS), ponesimod (OPTIMUM),
plus the *failures and harms*, which matter more: lenercept, anti-IFN-γ, atacicept.
Aggregate arm-level numbers are public — ClinicalTrials.gov results + the papers.
**This is the one that turns the gate from a restatement into a test.** Do it before (4).

### (6) Magnitude calibration needs data the repo can legally ship — *external, 2 routes*

Patient-level MS trial data **is** obtainable; the constraint is licensing, not access.
Verified 2026-09-17:

| source | what | catch |
|---|---|---|
| **MSOAC Placebo Database** (C-Path) | 2,465 individual patient records, 9 trials, RRMS/SPMS/PPMS, longitudinal EDSS, relapse events, T25FW/9HPT/PASAT/SDMT, SF-36, CDISC SDTM v3.2. **No fee**, ~4 weeks, bio sketch + research plan | **placebo arms only** — no treatment data, no comparator, no imaging |
| **Vivli** | Sanofi MS portfolio (tolebrutinib GEMINI 1/2, HERCULES, frexalimab) | enclave-only |
| **YODA** (J&J) | 482 trials incl. ponesimod OPTIMUM; 95.5% of requests approved, median 82.5 days to access | enclave-only |
| **CSDR** | Novartis (fingolimod, siponimod) | enclave-only |
| EU CTIS portal | 11,000+ trials, protocols + CSRs | **no IPD** — dead end |
| MSBase | 52,000+ patients, 33 countries | contributor-gated to practicing neurologists — dead end |

**The collision:** everything with a drug arm is enclave-only. You analyse inside their
environment and you do not get a file, so `data/` cannot ship a loader that reproduces
from it. That breaks this repo's "fetch open data at runtime, everything reproducible"
contract. **MSOAC is the only source that can be wired into `data/`.**

**Two routes that clear this without anyone's permission:**

1. **MSOAC for the untreated arm.** Real EDSS trajectories to calibrate blocker (2)'s
   damage→disability map against Kotelnikova, and real plausibility bounds for
   `bricks/vpop.py` instead of invented ones. Free. Apply early — the 4-week clock runs
   while (1)–(5) get built.
2. **Published exposure-response models for per-drug potency**, feeding blocker (4)'s
   `strength` hook. These are patient-level data already distilled into citable
   parameters, portable into an open repo the way Weatherley was. Start:
   **Valenzuela et al. 2022**, *CPT Pharmacometrics Syst Pharmacol* — ponesimod E-R
   analysis from the phase III RMS study. PMC9574745

## 8.2 Suggested order

`(5)` → apply for MSOAC → `(3)` → `(1)` → `(2)` → `(4)`.

(5) is half an hour and makes every later measurement trustworthy. The MSOAC application
goes in next because it is the only step with a queue in front of it. (3) comes before
(4) because per-drug strengths are untestable until there are arms to hold out.

## 8.3 What stays out of reach

Not a blocker — a boundary, so nobody plans around it. Even with all six cleared, this
does not become a trustworthy predictor of a **new** drug's effect size. That needs
patient-level trajectories fitted per drug, and those live in enclaves that do not let a
reproducible open repo exist around them. What is reachable: an open, honestly-labelled
harness that ranks candidates, kills doomed ones early, and shows its work. Build that.
**Nothing in this repo is evidence about multiple sclerosis** — that line does not move.
