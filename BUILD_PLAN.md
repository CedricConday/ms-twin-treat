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

### 2026-09-20 — the goal changed; blocker (1) closed; blocker (6) retired

**The goal.** This repo is for **candidate generation**, not ranking. Score
existing drugs well enough to propose new mechanisms, and hand a target package
to someone with a lab. §8.3's boundary is unchanged and still holds: this does
not become a predictor of a new drug's effect size, and nothing here is
evidence about MS.

**The order in §8.2 was written for ranking and is wrong for screening.**
Blocker (1) moves from fourth to first — not because it got more urgent, but
because it is the *vocabulary*. `mechanism_to_params(mechanism, strength,
disrupts_regulation)` was the entire representation of a drug: one of two
classes, a float, a bool. A "new candidate" could only differ from an existing
one by a number. Screening needs named intervention points.

New order: **(1) → (4) → (2) → leave-one-mechanism-out → screen.**

#### (1) CLOSED (c525930) — `bricks/qsp_velez.py`

Vélez de Mendizábal et al. 2011 (*BMC Syst Biol* 5:114, PMC3155504) transcribed
from the paper **and** from the authors' own Vensim model file, shipped as
Additional file 2 and fetched via the Europe PMC supplementary package for
PMC3155504. Every rate is first-hand. Nothing tuned.

**The authors' model file has an error, and the paper's text catches it.** Three
Table-1-vs-MDL discrepancies are recorded in the module docstring. Two follow
the MDL (`d2 = 0.002` not `0.02`; resting pools start at 7.5 / 2.4, not 0). The
third is resolved *against* the MDL: it computes the damage drive as `E^n/a`,
but the body text states `E' = (E/a)^2` — a factor of a = 22800. The authors'
own inline comment on that variable was right and their equation was wrong.
Both readings stay selectable; a test pins that they differ only by that factor,
which is why no arm-vs-arm ratio was ever affected.

**Reproduces Figure 3, not just the numbers.** `alpha_R` is the paper's
health/autoimmunity axis; raising it moves the system from autoimmune to
homeostatic, monotonically, with a >10x damage swing. First test in this repo
checking a published *result* rather than a published constant.

**A drug is now a `MechanismProfile`** — multipliers on the model's own named
rates (`alpha_E`, `gamma_E`, `alpha_R`, `gamma_R`, `delta`, `naive_E`). Harm is
mechanistically expressible for the first time: suppress effectors *and* strip
regulation is one object, and whether that nets to harm is a prediction of the
dynamics rather than a constant anyone picked. Runtime 0.43s for a 5-year run,
~20x faster than the grounded ABM.

**Two traps, both tested, both would cost a day to rediscover:**
- The relapse detector must be scored on the **untreated arm's** baseline.
  Against its own median, a treated arm's threshold slides down with the
  treatment: every arm read 31–37 relapses, including ones whose effector load
  had collapsed. Fixed, Treg support now reads 11 against untreated's 31.
- **No carrying capacity on E.** The Treg feedback is the only thing bounding
  the effector population, so `gamma_R` above ~2 gives `dE/dt → (alpha_E - eta)*E`
  and unbounded growth. That is the published model's property, not a porting
  bug — but the damage number is then meaningless and a screen would rank it as
  spectacular harm. Runs halt at the regime ceiling and report `in_regime=False`;
  consumers branch on that flag, never on the number.

`validated=False`. Transcribing a model is not reproducing it.

#### What the port cannot express — check this before trusting a screen

**No CNS compartment, no trafficking.** Verified against the full text: zero
occurrences of "trafficking"; one incidental "blood-brain" (a reference title);
"migration" appears only inside the parameter names `gamma_E` / `gamma_R`
("Teff death, anergy and migration Rate"). Loss from the active pool is a single
lumped term. No lymph node, no periphery/CNS transition, and tissue damage is an
output with no feedback on the T-cell dynamics.

So **depletion, sequestration and transit blockade all collapse onto `gamma_E`**.
Alemtuzumab (anti-CD52 depletion), ocrelizumab (anti-CD20 depletion), fingolimod
and ponesimod (S1P egress block) and natalizumab (anti-α4-integrin transit
block) are mechanistically distinct and this model cannot tell them apart.
Mapping them all to `gamma_E` is defensible — each removes effectors from the
pool that does damage — but it is a LUMP and each arm must say so.

Consequences: the mapping proceeds with those arms flagged; and **the screen is
capped** — a candidate whose novelty is purely in trafficking is invisible here.
If that design space matters, the successor model (Martinez-Pasamar et al. 2013,
*BMC Syst Biol* 7:34, PMC3651362 — antigen-specific subpopulations and
microglia) is the next port, and that call belongs *before* the screen is built.

#### (6) RETIRED — do not start the MSOAC application

Re-verified 2026-09-20: MSOAC is still v1.0, 2,465 records, 9 trials, **placebo
arms only**, no 2026 expansion. It was only ever needed for EDSS calibration and
`bricks/vpop.py` plausibility bounds. On the ARR route below the EDSS half does
not arise, and the bounds half is partly covered with **zero queue** — the
CLARITY and ADVANCE synthetic placebo arms are open data on Figshare with Merck
and Biogen approval, carrying ARR, T2/Gd lesions and confirmed disability
worsening (PMC12488035).

#### The 2026-09-17 circularity wall has a way around it — the MRI channel

§8's suggested route for blocker (4), published exposure-response models, is
circular where the response is ARR. It is **not** circular on the MRI channel,
and that channel is fully public:

- Every quantified arm in `docs/TRIAL_ANCHORS.md` comes from a trial that also
  reports arm-level new/enlarging T2 and Gd-enhancing lesion outcomes.
- Per-drug exposure-response models on that endpoint already exist for
  **natalizumab** (log-linear on Gd-enhancing lesion count, *J Pharmacokinet
  Pharmacodyn* 2017, doi:10.1007/s10928-017-9514-4) and **ponesimod** (CUALs,
  Valenzuela et al. 2022, PMC9574745).

Fit `strength` on lesion effect, predict ARR. The gate never sees the arm's own
relapse number, so the circularity is gone.

**Known hole — state it before building on it.** Lenercept reported **no
significant MRI difference** while relapse rate rose significantly (p=0.007) and
relapses were more severe and longer (*Neurology* 1999;53:457, PMID 10449104).
The MRI channel is blind to the one harm case the model most needs to get right,
so MRI-fitted potency cannot be the only input.

#### (2) has a published map on the right endpoint

**Sormani & Bruzzi 2013**, *Lancet Neurol* 12(7):669–76, PMID 23743084 — 31
trials, 18,901 RRMS patients. Trial-level regression of treatment effect on MRI
lesions against treatment effect on relapses: **slope 0.52, R² = 0.71**.

A fitted micro→clinical map scored on **ARR**, which retires §8's
"decide before starting" catch: no EDSS endpoint, and the Kotelnikova
damage→EDSS mismatch does not arise. **Read the paper for the exact regression
form and scale before implementing** — log relative-rate against log
relative-rate is the likely form but has NOT been verified.

#### (4) first half DONE (ac04821) — `bricks/profiles.py`, and the ceiling it measures

Each of the 14 arms now carries a `MechanismProfile` over the Vélez model's own
rates, assigned from pharmacology and never from the trial. **Sourced: which
points each drug touches and in which direction. Stubbed: how much** — every
non-unit multiplier is one shared constant, because a hand-picked number per
drug is indistinguishable from fitting. `with_magnitudes()` is the migration
path off the stub and refuses to add or remove a point, since that is a
mechanism claim.

Two of three checked assignments came back against the obvious answer:
- **alemtuzumab gets NO Treg-depleting axis.** CD52 is on Tregs, so "it depletes
  Tregs too" is tempting and has the sign backwards — repopulation is
  Treg-BIASED, regulatory cells return faster than effectors with increased
  suppressive capacity (PMC4519957, PMC8581537).
- **teriflunomide's Treg axis is left unset**, evidence contested: Klotz 2019
  (doi:10.1126/scitranslmed.aao5563) reports no change, PMID 40879143 reports
  impaired FOXP3+ function.
- glatiramer acts at antigen presentation (MHC II competition, Arnon & Aharoni
  PNAS 2004) plus Treg restoration → `delta` + `alpha_R`; the only arm whose
  primary mechanism is at the APC step.

2 classes → **6 distinct intervention patterns**. Lenercept is one object that
suppresses *and* de-regulates, so its harm is now a prediction.

**THE FINDING, and it re-orders what is left.** `degeneracy_report()` measures
the ceiling: **six of the thirteen treated arms collapse onto `gamma_E`** —
natalizumab, fingolimod, ponesimod, ocrelizumab, alemtuzumab, atacicept. 15
indistinguishable pairs, all model ceiling, versus 4 stub-only pairs that
dissolve when magnitudes land. A test pins the sharpest case: **natalizumab
(−68%, AFFIRM) and atacicept (raised relapses, ATAMS halted) are byte-identical
in this model.** No potency number and no leave-one-out can make one right
without making the other wrong.

So the Martinez-Pasamar 2013 port is **not optional and not later** — it is on
the critical path, because the arms it would separate are the high-efficacy
drugs and one of the two harm cases. Building LOMO and the screen on the
current representation would grade a model that is provably unable to pass.

**Revised order: Martinez-Pasamar port → (4) magnitudes → (2) Sormani → LOMO →
screen.**

#### (2) CLOSED (4441b0c) — `bricks/sormani.py`

Sormani & Bruzzi 2013 (*Lancet Neurol* 12(7):669-76, PMID 23743084), 31 trials,
18,901 RRMS patients, weighted linear regression of **log-transformed relative**
treatment effects:

    log(RR_relapse) = intercept + 0.52 * log(RR_lesion)      R^2 = 0.71

Scored on ARR, so §8's "decide EDSS vs ARR before starting" catch does not arise.

Two things stated rather than smoothed: the **intercept is not published where
this repo can read it** (slope is in the abstract, intercept is behind the
paywall), so it defaults to 0 because the relation's boundary condition forces
it, and `INTERCEPT_IS_ASSUMED` rides into every prediction. And the map has a
documented failure — lenercept, which reported no significant MRI difference
while relapses rose (p=0.006). Fed a null lesion effect the map says "no change",
confidently and wrongly, so `predict_relapse_ratio` returns an object carrying
`blind_spot=True` for small lesion effects rather than a bare float.

Sanity check not part of any fit here: AFFIRM reported ~90% fewer Gd-enhancing
lesions and 68% fewer relapses; the map returns -69.8%.

Applied at the ARM level only (`readout.relapse_ratio_from_arms`). Sormani is a
between-trial relation; using it per patient would read a trial-level regression
as an individual one. `MAX_RELAPSE` stays for the per-patient proxy and is
documented as superseded for ratio scoring — it always cancelled in the ratio,
which is why it was never the thing to fix.

#### THE VARIANCE FINDING — read before trusting any cohort number in this repo

Building the LOMO response table surfaced something that invalidates any
small-cohort result on the ported QSP, including the first version of that table.

Re-derive with `scripts/measure_qsp_variance.py --n 200`. Measured 2026-09-20 at
a 730-day horizon: **median 2.14, mean 74.9, SD 612, range 0.089 to 8118** — a
~90,000-fold spread across stochastic infection histories, with the top decile
of runs holding 96% of all damage and the standard error at n=4 four times the
mean. (An earlier 40-run sample suggested "266-fold"; it undersampled the tail.
Quote the script, not that.)

The first LOMO table used 4 seeds and a mean. It produced a confident headline
made entirely of its own noise, and it also produced two "findings" that did not
survive re-measurement — that `gamma_E` has the wrong sign, and that `naive_E`
does. **Neither is established; both readings came from the noise.** Do not
quote them.

What this forces:
- the statistic is the **median of paired per-seed ratios**, never the mean;
- the cohort is **128 seeds**. Bootstrapped CV of the median: 628% at n=4,
  41.8% at n=16, 21.2% at n=48, **13.7% at n=128**. That ~14% is the noise floor
  under every number from this model;
- **pairing does not rescue a small cohort.** Sharing a seed between the treated
  and untreated arms looks like it should cancel the infection history; measured,
  it narrows the ratio's IQR by ~1.1x. Three docstrings claimed otherwise before
  it was measured. Kept because it is free, but the cohort size is what makes a
  comparison valid;
- the horizon is **730 days**, matching the trials rather than the model's
  five-year default, and ~10x cheaper per run.

Anything elsewhere in this repo that compares arms on a handful of QSP runs is
suspect until re-measured this way.

#### (4) CLOSED both halves — `bricks/profiles.py` (f5141d3) + `backtest/potency.py` (12d4627)

**Directions**, from pharmacology, never from outcome. 2 classes -> 7 distinct
intervention patterns. Lenercept is one object that suppresses AND de-regulates.

**Magnitudes**, from the MRI channel, never from ARR. `sormani.invert()` reads a
lesion ratio back out of the relapse-unit response table (the regression is
monotone in log space, so it inverts exactly), and the potency that reproduces a
trial's observed lesion ratio becomes the drug's magnitude. A test monkeypatches
`KNOWN_OUTCOMES` to raise if a fit touches it, so non-circularity is enforced
rather than asserted.

**The cross-check disagrees, and that is the most useful number here.**
Ocrelizumab is the only arm with a magnitude from two independent sources:

    ke = 0.85   Martinez-Pasamar 2013, EAE mouse T-cell dynamics
    ke = 0.40   OPERA I Gd-enhancing lesion ratio, through this module

**2.1x apart.** Neither is adjusted toward the other; a test pins the gap.
`profiles.py` keeps the EAE value as default because it measures the parameter
directly rather than inverting through two models. Treat that factor of two as
the accuracy of the whole potency layer.

**Not done: eight of nine arms have no MRI number.** Paywalled, and not in the
abstracts Europe PMC serves. `PENDING_EXTRACTION` lists them and a test asserts
every quantified arm is either fitted or explicitly pending.

#### LOMO CLOSED (28b3c4c) — `backtest/lomo.py`

Holds out a whole mechanism group and predicts it blind, on the screening stack
(qsp_velez -> profiles -> sormani), which is also a check that those compose.

**Out-of-sample MAE 71.0pp against a 14.4pp predict-the-mean null.** Per fold:
alpha_E 38.2, alpha_R|delta 29.0, gamma_E 112.4, ke 46.0. The model does not
beat the null on an unseen mechanism. **That is the gate on the screen** — until
it inverts, a generated candidate's score is not information.

`_profile_at` deliberately discards a fitted magnitude where one exists:
ocrelizumab's ke=0.85 would hand the held-out mechanism a number somebody
already measured for it.

**gamma_E predicts +60.3% where its four trials report -30 to -68%**, and that
survived the 128-seed rebuild, so it is the model and not the noise. Reading the
equations: `gamma_E` is gated by the Treg Hill term, making it the
regulatory-killing channel rather than a generic death rate, so multiplying it
misrepresents an antibody that depletes independently of regulation. Recorded,
not patched — the published model has no additive loss term to move it to. This
is the single largest contributor to the LOMO failure.

#### Still open

- **MRI extraction for eight arms** — PRISMS, CONFIRM, AFFIRM, FREEDOMS, TEMSO,
  DEFINE, CARE-MS I, OPTIMUM. New/enlarging T2 and Gd-enhancing counts, treated
  arm and comparator, as reported. Needs full-text access this box does not have.
- **a non-MRI harm channel.** Lenercept fits to potency ~0 through MRI, for a
  drug that harmed people. No amount of extraction fixes that.
- **`gamma_E` has no correct home.** Five arms sit on a dial that predicts the
  wrong sign. Fixing it means an additive effector-loss term the published model
  does not have, i.e. extending the model rather than porting one.
- **the screen** (`screen/`), a generator over the QSP's intervention points.
  Gated on LOMO beating its null — until then it ranks noise.
- **wiring.** `qsp_traj` is written but nothing consumes it; the clinical gate
  runs through `abm_damage`. **The port does not move the gate on its own** —
  readout must read `qsp_damage` before any of this reaches the score.
