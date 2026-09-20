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

### (5) Readout rounding quantizes the gate — ~~*30 minutes*~~ **DONE 2026-09-17 (3fcfb2c)**
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

### (3) Four arms is too few to test anything — *~2–3 days* — **ARMS DONE 2026-09-17, LOO STILL OPEN**
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


## 8.4 Progress log — append dated entries, never rewrite

### 2026-09-17 — blocker (5) closed
`bricks/readout.py` now stores the exact `relapse_proxy` and rounds only at the display
(3fcfb2c). Regression test `test_readout_relapse_proxy_is_not_quantized` fails on the old
code and passes on the new. The clinical gate's direction and magnitude are unmoved by it,
which is the expected result: the fix makes the measurement trustworthy, it does not close
a gap.

### 2026-09-17 — blocker (3), first half: 4 arms → 14
Ten arms added to `bricks/intervention.py`, each classed from its pharmacology and never
from its outcome; anchors, comparators and PMIDs in `docs/TRIAL_ANCHORS.md`, verified
against PubMed the same day. `backtest/clinical.py` now scores each arm against **the
comparator its trial actually used** — OPERA and CARE-MS I against IFN-β, OPTIMUM against
teriflunomide — instead of converting active-comparator results into invented
"vs placebo" figures. Arms with identical intervention parameters are simulated once, so
the 14-arm gate costs three pipeline runs, not fourteen.

**Measured: DIRECTION 9/14, MAGNITUDE 1/9** (was 4/4 and 0/2 on four arms; the model is
unchanged). Failing arms: ocrelizumab, alemtuzumab, ponesimod (all +0% vs their
comparator), lenercept and atacicept (−77% where the trials harmed).

What it proves, and this is the reason the arm set was grown:
- **the rule's shape is wrong, not just its calibration.** Lenercept and atacicept are
  immunosuppressive by mechanism and harmed patients. No `SUPPRESSIVE_STRENGTH` value
  fixes a rule that maps class → direction.
- **blocker (4) is now measurable.** One strength per class returns exactly 0% on every
  within-class comparison, which is the `strength` hook's absence stated as a number.
- the single magnitude hit is natalizumab (−77% vs −68% reported); the same −77% misses
  IFN-β, glatiramer and teriflunomide by 45pp or more.

**Still open in (3):** the leave-one-arm-out test. It needs per-drug strengths to fit, so
it lands with blocker (4), not before. The arm set is now large enough to hold one out.

**Order from here:** (4) with per-drug potency from independent exposure-response data,
then the LOO that (3) still owes, then (1) and (2). The MSOAC application remains the only
step with a queue in front of it and has not been started.

### 2026-09-17 — what the rounding was actually worth, measured
Blocker (5) quantified, by running the gate twice on the same code with only the 2dp
rounding switched back on: the suppressive arms move **−77.0623% exact vs −76.7045%
rounded**, and the immunogenic arms **+28.0902% vs +28.4091%** — about **0.32–0.36pp**
on the headline. The per-patient discretization is the ~25% the blocker describes, but
it largely averages out across a cohort of 12, so at a fixed `MAX_RELAPSE` the gate-level
effect is a third of a percentage point.

That refines, rather than contradicts, the note under blocker (2): the 0.8 / 1.5 / 3.0
sweep spread ~2.9pp because each `MAX_RELAPSE` value puts the per-patient proxy on a
different quantization grid, not because rounding is worth 2.9pp at any one value.
Either way the conclusion stands — the constant cancels in the ratio, and neither it nor
the rounding can move a 45pp magnitude gap.

### 2026-09-17 — the rule's shape, not its calibration
`bricks/grounding.py` gained a per-drug `disrupts_regulation` flag so a suppressive drug
whose target also carries a regulatory/reparative function gets the harm channel too.
Lenercept is flagged on TNF biology (Liu 1998 PMID 9427610, Arnett 2001 PMID 11600888);
atacicept is deliberately NOT flagged (no independent source for BAFF/APRIL; the nearest
is anti-CD20 Bregs, PMID 18802481, different target and timing-dependent).

**Measured: lenercept −77% → −35%, sign unchanged; gate still 9/14 and 1/9.** No third
constant was added — a flagged drug reuses `IMMUNOGENIC_STRENGTH`, and a test asserts it,
because a harm constant tuned until lenercept reads harmful is fitted to the exam.

Consequence for blocker (4): the suppressive/harm *ratio* is now the load-bearing number,
and both sides of it need independent sources. Per-drug potency from exposure-response
work remains the next step; `mechanism_to_params(strength=...)` is still unfed.

### 2026-09-17 — blocker (3) closed: leave-one-arm-out is running
`backtest/response_curve.py` tabulates the pipeline once (treat -> simulated relapse
change, 0.05 grid, both harm levels, cached in `results/response_curve.json`), so fitting
is a search over the curve instead of thousands of runs. Interpolation error measured at
0.8pp in the steepest region and written into the module docstring.

`backtest/loo.py` holds out each quantified arm, fits the class strength on the other
eight, predicts the held-out one, and scores it against a predict-the-mean null.

**Result: out-of-sample MAE 26.0pp against a null of 14.0pp. The rule does not beat the
null.** Placebo-controlled arms alone: 17.1pp vs 16.0pp — still behind. The fitted
strength lands at 0.42–0.50 on every fold, i.e. wherever the training arms' average sits,
and the three active-comparator arms predict exactly 0% at every strength.

That is the honest headline for the whole build: **a single class strength carries no
drug-specific information.** The gate at 9/14 direction flattered it; the LOO does not.
Per-drug potency (blocker 4) now has a scoreboard to beat, and `tests/test_loo.py` pins
the active-comparator arms at exactly 0 so that when potency lands, the test fails and
says so.

**Obstacle found while looking for that potency, stated before anyone spends a day on it:**
the route §8 suggests — published exposure-response models — is CIRCULAR where the
"response" is the annualized relapse rate, which is the outcome the gate predicts.
Pharmacodynamic potency avoids that, but no single PD axis spans these mechanisms:
natalizumab *raises* circulating lymphocyte counts (cells stay in blood), anti-CD20
depletes B cells specifically, teriflunomide barely moves ALC. A per-drug number needs a
bridging assumption across biomarkers, and that assumption is the next thing to argue
about, not to quietly pick.

### 2026-09-20 — the goal changed, and the order with it

Stated this session: the point of this repo is **candidate generation**, not a
ranking harness. Score existing drugs well enough to propose new ones, hand a
target package to someone with a lab. §8.3's boundary is unchanged and still
holds — this does not become a predictor of a new drug's effect size, and
nothing here is evidence about MS — but the build order §8.2 gives was written
for ranking and is wrong for screening.

**Blocker (1) moves from fourth to first.** Not because it got more urgent but
because it is the *vocabulary*. `mechanism_to_params(mechanism, strength,
disrupts_regulation)` was the entire representation of a drug: one of two
classes, a float, a bool. You cannot propose a novel mechanism in that; a
"new candidate" could only differ from an existing one by a number. Screening
needs a model with named intervention points, which is what the port provides.

New order: **(1) -> (4) -> (2) -> leave-one-mechanism-out -> screen.**

#### (1) DONE 2026-09-20 (c525930) — `bricks/qsp_velez.py`

Vélez de Mendizábal 2011 transcribed from the paper *and* the authors' Vensim
model file (Additional file 2, via the Europe PMC supplementary package for
PMC3155504). Nothing tuned. Three Table-1-vs-model-file discrepancies are
recorded in the docstring; one is resolved *against* the model file, because
the paper's body text says `E' = (E/a)^2` while the MDL's equation computes
`E^n/a` — the authors' own inline comment was right and their equation is a
transcription error. Both readings remain selectable.

Reproduces the paper's Figure 3 result, not just its numbers: `alpha_R` is the
health/autoimmunity axis, and raising it moves the system from autoimmune to
homeostatic monotonically with a >10x damage swing. That is a test against a
published *result*, which nothing in this repo had before.

A drug is now a `MechanismProfile` — multipliers on the model's own named
rates (`alpha_E`, `gamma_E`, `alpha_R`, `gamma_R`, `delta`, `naive_E`). Harm is
expressible mechanistically for the first time: suppress effectors *and* strip
regulation is one object, and whether it nets to harm is a prediction of the
dynamics rather than a constant anyone picked.

Runtime 0.43s for a 5-year run, ~20x faster than the grounded ABM.

**Two traps, both now tested, both would have cost a day each to rediscover:**
- The relapse detector must be scored on the **untreated arm's** baseline. Against
  its own median, a treated arm's threshold slides down with the treatment and
  reports an unchanged relapse count for a cohort whose effector load collapsed.
  Before the fix every arm read 31-37 relapses; after it, Treg support reads 11
  against untreated's 31, and the healthy configurations read 0.
- **The published model has no carrying capacity on E.** The only thing bounding
  the effector population is the Treg feedback, so `gamma_R` above ~2 gives
  `dE/dt -> (alpha_E - eta)*E` and unbounded growth. That is the model's property,
  not a porting bug — but the damage number is then meaningless and a screen
  would rank it as spectacular harm. Runs halt at the regime ceiling and report
  `in_regime=False`; consumers must branch on that flag, never on the number.

`validated=False`. Transcribing a model is not reproducing it.

#### The circularity wall from 2026-09-17 has a way around it

§8's suggested route for blocker (4) — published exposure-response models — was
found circular where the response is ARR. It is not circular on the **MRI
channel**, and that channel is fully public:

- Every quantified arm in `docs/TRIAL_ANCHORS.md` comes from a trial that also
  reports arm-level new/enlarging T2 and Gd-enhancing lesion outcomes.
- Per-drug exposure-response models on that endpoint already exist for
  **natalizumab** (log-linear on Gd-enhancing lesion count, J Pharmacokinet
  Pharmacodyn 2017, doi:10.1007/s10928-017-9514-4) and **ponesimod** (CUALs,
  Valenzuela et al. 2022, PMC9574745).

Fit `strength` on lesion effect, predict ARR. The gate never sees the arm's own
relapse number, so the circularity is gone.

**Known hole, state it before building on it:** lenercept. That trial reported
**no significant MRI difference** while relapse rate rose significantly
(p=0.007) and relapses were more severe and longer (Neurology 1999;53:457,
PMID 10449104). The MRI channel is blind to the one harm case the model most
needs to get right, so MRI-fitted potency cannot be the only input.

#### Blocker (2) has a published map on the right endpoint

**Sormani & Bruzzi 2013**, *Lancet Neurol* 12(7):669-76, PMID 23743084 — 31
trials, 18,901 RRMS patients, trial-level regression of treatment effect on MRI
lesions against treatment effect on relapses: **slope 0.52, R^2 = 0.71**.

That is a fitted micro->clinical map scored on **ARR**, which retires §8's
"decide before starting" catch entirely: no EDSS endpoint is needed, and the
Kotelnikova damage->EDSS mismatch does not arise. Before implementing, read the
paper for the exact regression form and scale (log relative-rate vs log
relative-rate is the likely form but has NOT been verified — do not assume it).

#### Blocker (6) is retired

MSOAC re-verified 2026-09-20: still v1.0, 2,465 records, 9 trials, **placebo
arms only**, no 2026 expansion. It was only ever needed for EDSS calibration
and `bricks/vpop.py` plausibility bounds. On the ARR route the EDSS half does
not arise, and the bounds half is partly covered with **zero queue** — the
CLARITY and ADVANCE synthetic placebo arms are on Figshare as open data with
Merck and Biogen approval, carrying ARR, T2/Gd lesions and confirmed disability
worsening (PMC12488035). **Do not start the MSOAC application.**

#### Still open

- **(4)** per-drug potency from the MRI channel. Needs arm-level lesion
  outcomes extracted from the 9 quantified trials into `docs/TRIAL_ANCHORS.md`,
  then `mechanism_to_params(strength=)` fed from them. The lenercept hole above
  means the harm channel still needs a non-MRI source.
- **(2)** the Sormani map in `bricks/readout.py`, replacing `MAX_RELAPSE` /
  `MAX_LESIONS`. Verify the regression form first.
- **leave-one-mechanism-out** (`backtest/lomo.py`), beside the existing LOO.
  Holding out a drug from a class still present in training does not grade a
  screen; holding out a whole mechanism does. This becomes the headline number.
- **the screen** (`screen/`), a generator over the QSP's intervention points.
  Gated on LOMO beating its null — until then it would rank noise.
- **wiring**: `qsp_traj` is written but nothing consumes it. The clinical gate
  runs through `abm_damage`, so the port does not move the gate on its own.
  Readout must start reading `qsp_damage` for any of this to reach the score.

#### Limit of the ported model, found 2026-09-20 before mapping the arms

**Vélez de Mendizábal has no CNS compartment and no trafficking.** Checked
against the full text: zero occurrences of "trafficking", one incidental
"blood-brain" (a reference title), and "migration" appears only inside the
parameter names `gamma_E` / `gamma_R` ("Teff death, anergy and migration
Rate"). Loss from the active pool is a single lumped term; there is no
lymph-node compartment, no periphery/CNS transition, and tissue damage is an
output variable with no feedback on the T-cell dynamics.

Consequence for the arm set, and it is not small: **depletion, sequestration
and transit blockade all collapse onto `gamma_E`.** Alemtuzumab (anti-CD52
lymphocyte depletion), ocrelizumab (anti-CD20 B-cell depletion), fingolimod and
ponesimod (S1P modulators, egress block from lymph nodes) and natalizumab
(anti-alpha4-integrin, BBB transit block) are mechanistically distinct drugs
that this model cannot tell apart. Mapping them all to `gamma_E` is defensible
— each removes effectors from the pool that does damage — but it is a LUMP, and
it must be recorded as one on each arm rather than discovered later.

Two things follow:
- The mapping can proceed. Flag the lumped arms explicitly; do not pretend the
  model distinguishes them.
- **It caps the screen.** A candidate whose novelty is purely in trafficking is
  invisible here. If that is the design space worth searching, the successor
  model — Martinez-Pasamar et al. 2013, *BMC Syst Biol* 7:34, PMC3651362, which
  adds antigen-specific subpopulations and microglia — is the next port, and
  that decision should be made before the screen is built, not after it returns
  candidates it cannot rank.
