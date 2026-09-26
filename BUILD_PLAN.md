# BUILD_PLAN — ms-twin-treat (cold-boot resume doc)
*Written 2026-08-19, pre session-reset. If you are a fresh boot: read this top-to-bottom, run the COLD-BOOT checks, then execute BUILD ORDER. Goal = the 8 remaining bricks BUILT (loading + emitting output through the spine), not validated. ~4h wall-clock, install/compile-bound.*

---

## 0. STATUS (what's already real)
- **Repo:** `/home/ubuntu/repos/ms-twin-treat` — local git, 2 commits, author `Cedric Conday <cedric@condaydigital.com>`. **NOT pushed** (push held for 1300, Cedric's go only).
- **Built + validated:** brick #0 backtest harness, brick #1 Kang data loader. **2 of 10.**
- **Real result on real data:** Kang 2018 IFN-β PBMCs, 8 cell types. identity null = 0.00, global-mean-shift null = **0.85 aggregate delta_pearson** (the bar). Megakaryocytes = 0.48 (the honest outlier). No model beats the bar yet.
- **Brain repo (research):** ~~`/home/ubuntu/repos/ms-twin/`~~ **folded into `docs/research/` on 2026-09-25; the separate repo is retired** — `RESEARCH_FINDINGS.md` (all bricks, licenses, sources), `PROJECT_MS_TWIN.md`, `REVIEW_CHARGE.md`.

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

~~**DIRECTION 4/4. MAGNITUDE 0/2.**~~ **SUPERSEDED — re-measured 2026-09-20 on the
17-arm set (12 quantified), every number below run fresh, none recalled:**

| gate | command | result |
|---|---|---|
| ABM path | `backtest.clinical` | **10/17 direction, 1/12 magnitude** |
| grounded stack | `backtest.clinical_velez` | **5/16 direction, 2/11 magnitude** |
| leave-one-ARM-out | `backtest.loo` | **28.3pp MAE vs 11.5pp null** (placebo-only 16.1 vs 15.2) |
| leave-one-MECHANISM-out | `backtest.lomo` | **45.9pp MAE vs 12.3pp null** |
| LOMO + capacity | `backtest.lomo_capacity` | **45.6pp vs 12.3pp null** |
| per-drug potency | `backtest.potency` | 4 arms OUT OF RANGE (fingolimod, ponesimod, cladribine, daclizumab) |
| spine | `spine/run_demo.py` | **8 of 15** state keys unvalidated (§8.0 above says 8/13; the spine grew) |

**Every gate that has a null loses to it.** The §8.4 entries below record 9/14,
1/9, 5/13, 26.0/14.0 and 71.0/14.4 — all of them true when written and all of
them one arm-set behind, because `8dc2472` grew the arms afterwards. Three
places in the CODE also printed stale numbers as if measured; fixed the same
day (`clinical_velez.py` no longer prints the other gate's headline at all,
because a hardcoded measurement goes stale silently).

The two suppressive arms being byte-identical is still half the problem.

## 8.1 The six blockers

> **STATUS 2026-09-20 — ALL SIX ARE CLOSED OR RETIRED. Read this section as the
> record of what they were, not as a worklist.** (1) closed c525930, (2) closed
> 4441b0c, (3) closed f9b0c88 (LOO) + 28b3c4c (LOMO), (4) closed f5141d3 +
> 12d4627, (5) closed 3fcfb2c, (6) **retired** — do not start the MSOAC
> application, see §8.4. Closing them did not make the build work: what remains
> is in "Still open" at the end of §8.4, and it is a model property, not a
> backlog. Verified against the tree on 2026-09-20, not recalled.

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
- **Martinez-Pasamar et al. 2013**, *BMC Syst Biol* 7:34 — ~~the successor, adding
  antigen-specific subpopulations and microglia~~ **WRONG, STRUCK 2026-09-20: it is
  the SAME four-ODE system, re-parameterised for mouse EAE. No B-cell and no
  microglia state variable. See §8.4.** PMC3651362

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

> **STRUCK 2026-09-20 — every step below has been taken or retired, and the one
> instruction in it is now WRONG: do NOT apply for MSOAC (§8.4 retires it; the
> CLARITY and ADVANCE synthetic placebo arms on Figshare, PMC12488035, carry ARR
> and lesion counts with zero queue). Kept for the record.**

~~`(5)` → apply for MSOAC → `(3)` → `(1)` → `(2)` → `(4)`.~~

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

~~So the Martinez-Pasamar 2013 port is **not optional and not later** — it is on
the critical path, because the arms it would separate are the high-efficacy
drugs and one of the two harm cases.~~ **STRUCK 2026-09-20 (see §8.4, "the
critical path was a mirage"): MP2013 has the same four state variables and the
same six dials, so it separates nothing.** The rest stands: building LOMO and
the screen on the current representation grades a model that is provably unable
to pass.

~~**Revised order: Martinez-Pasamar port → (4) magnitudes → (2) Sormani → LOMO →
screen.**~~ **STRUCK 2026-09-20 — the port is worthless, see §8.4. (4), (2) and
LOMO all landed anyway; what the order was reaching for, mechanism separation,
has no published model behind it yet.**

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

#### WIRING CLOSED (2ef3e53) — and the grounded stack scores WORSE

`backtest/clinical_velez.py` scores the same arms against the same cited
outcomes with the same thresholds, through qsp_velez -> profiles -> sormani.

    ABM path (`clinical.py`)             9/14 direction   1/9 magnitude
    grounded stack (`clinical_velez.py`) 5/13 direction   2/9 magnitude

Transcribing a published model, assigning intervention points from pharmacology
and adding a meta-analysis readout bought one magnitude hit and lost four
direction calls. **Both gates are kept and both are reported.** A new number
that quietly replaces an old one is not a comparison.

Two rows before anyone quotes the headline:
- **atacicept passes for the WRONG REASON.** It sits on `gamma_E`, `gamma_E` has
  the wrong sign, and atacicept harmed patients — a broken dial and a harmful
  drug cancelling. Fixing the defect should flip it to a miss.
- **lenercept is UNDEFINED, not wrong.** `alpha_E` down with `alpha_R` down
  strips regulation until the effectors run away; all 128 histories leave the
  regime. Counted as a miss, because an unscoreable arm must not leave the
  denominator.

**`gamma_E` is now the dominant defect in the repo**: five of thirteen arms, the
largest LOMO fold error (112pp), and the reason this gate went backwards.

#### The virtual population was decorative (28a9cab)

Found while asking whether the LOO's 26.0-vs-14.0 survived a different cohort
draw. It had never seen one. Three things stacked:

1. `sample_vpop(seed=S)` assigned each patient `"seed": i` — the plain index —
   so every cohort handed the ABM the same twelve simulation seeds.
2. The sampled `qsp_params` never reach the scored path; the readout scores
   `abm_damage`.
3. `bbb_disruption` is inert because **no arm in the library sets
   `cns_required`**, so the barrier term is always zero.

Measured before: untreated proxy **0.147137 for every seed 1-8, SD exactly
0.0000**. After: 0.142712 / 0.157593 / 0.148240, about ±10%. The clinical gate
is unmoved at 9/14 and 1/9, which is the right outcome — direction scoring
should not be fragile to the draw.

`results/response_curve.json` was built under the old seeding and is stale;
rebuilt, and the LOO re-run against it is the first honest reading of that
number.

#### vpop ported onto the grounded model (15c7765)

Swapping the pipeline to the port left `vpop` filtering candidates through the
TOY, so "plausible patient" meant plausible under a model nothing runs.
`sample_vpop_velez` runs the same Allen-Rieger method against the model in use,
and the bounds are **the paper's own**: Table 1 states `alpha_E` and `alpha_R` as
intervals (`[1:2]`, `[0.25:2]`) because those are the ranges the authors swept.

It rejects 57 of 65 candidates, and the survivors land at `alpha_R` 0.27-0.37 on
their own. The filter is told only "must develop disease and stay in regime" and
it recovers the low-`alpha_R` autoimmune configuration the paper describes —
a small independent check that the port behaves as documented.

#### The LOO headline SURVIVES the cohort fix (2026-09-20)

`results/response_curve.json` rebuilt under the corrected seeding, LOO re-run:

    out-of-sample MAE 26.0pp   predict-the-mean null 14.0pp
    placebo-controlled arms only: 17.0pp vs 16.0pp

Unchanged to the decimal (placebo-only moved 17.1 -> 17.0). The number was worth
distrusting and it holds.

**Why it holds is the point.** The LOO does not fail because of noise, so a
different cohort cannot rescue it. It fails structurally: one class strength
carries no drug-specific information, and the three active-comparator arms
predict exactly 0% at every strength because both arms of those trials get the
same value. No cohort draw changes either fact.

#### The depleting class cannot come out beneficial (ebbfb57)

`gamma_E`'s wrong sign was blamed on `gamma_E`. It is not `gamma_E`, it is the
model, and it explains the LOMO failure and the 5/13 direction score at once.

Damage is driven by `(E/a)^2`, so it is set by effector PEAK EXCURSIONS rather
than effector load. Over 48 histories at 730 days:

    arm                        med E    med PEAK E    damage
    untreated                   1126         52741     1.457
    alpha_E x0.5 (damp growth)  1020         30794     0.665
    gamma_E x1.5 (kill cells)   1077        192455    10.249

Median effector load is nearly identical; the PEAK differs six-fold and damage
follows the peak. Mechanism: effectors recruit their own regulators via
`E^h/(ke^h+E^h)`, so killing them lowers recruitment, `R` falls, the
proliferation brake releases, and the population rebounds into a larger
excursion than the killing removed.

Tested rather than argued: an additive regulation-independent loss term
`-depletion_rate * E` was added as an explicit extension whose own comment said
to REMOVE it rather than tune it if it failed. It failed (median damage
1.46 -> 4.69 -> 11.3 -> 39.8, in-regime runs 64 -> 26), so it was removed and
`qsp_velez.py` is pure transcription again.

**This is a boundary, not a bug.** Natalizumab, fingolimod, ponesimod,
alemtuzumab and atacicept cannot come out beneficial in this model under any
mapping — five of the nine quantified arms. A screen over it is blind to the
whole depleting / sequestering / trafficking class, and that has to be decided
about before the screen is built.

#### IDEAS RECOVERED FROM THE RESEARCH REPO (2026-09-20)

`~/repos/ms-twin/docs/RESEARCH_FINDINGS.md` is named on line 10 of this file as
the brain repo, and it had not been read in any recent session. Checked against
the build repo, several named, concrete resources were never carried across.
Listed so they stop being lost:

- **DREAM(ZS) / simulation-based inference.** `RESEARCH_FINDINGS.md:39` defines
  the differentiator as the plausible-patient method **plus modern sampling**.
  Only the first half exists. `bricks/vpop.py` samples by Latin hypercube with a
  rejection filter — no posterior, no DREAM(ZS), no SBI. The gap matters because
  rejection sampling yields an accepted SET while the method's value is a
  prevalence-WEIGHTED population; `weight_to_prevalence()` is a single-axis
  stand-in for MAPEL, not MAPEL.
- **Verscheijden 2019**, a 14-compartment brain PBPK **with runnable R in the
  supplement** — `RESEARCH_FINDINGS.md:17` calls it "the most immediately
  runnable open CNS model". §6 of this file deliberately deferred it ("don't
  install R for Verscheijden now") and the barrier brick has been a stand-in
  ever since. That is the same shape as blocker (1), which closed by porting a
  published model, and the barrier is currently inert anyway because no arm sets
  `cns_required`.
- **ADMET-AI** (pip, MIT) — a BBB permeability gate. Directly relevant to the
  mechanism-to-molecule gap, never installed or referenced.
- **`BMSQSP/QSPToolbox`** — the MATLAB + SimBiology reference implementation of
  the very method `vpop.py` claims to be the first open port of. Never consulted.
- **Allen-Rieger-Musante 2016, doi:10.1002/psp4.12063** — the method's citation.
  `vpop.py` named the method for a month without citing it, while claiming to be
  its first open implementation. Added.
- **scArchon** — the perturbation-prediction benchmark the harness is modelled
  on. It was cited in `backtest/README.md` and **I deleted it in this session**
  while rewriting that file. Restored.

The last one is the useful warning: a rewrite that improves a document can still
lose the one line in it that was a citation. Check what a rewrite drops, not
just what it adds.

#### THE HARM CHANNEL: EAE IS EMPIRICALLY DEAD (2026-09-20)

Searched for an independent, non-MRI source that predicts harm, because the MRI
channel fits lenercept to "does nothing" for a drug that harmed people. The
obvious candidate was EAE — it is what `bricks/profiles.py` already cites for
lenercept, and TNF-deficient mice getting severe EAE is exactly the signal MRI
missed.

**The best systematic evidence says EAE does not predict.**

    Berg I, Härvelid P, Zürrer WE, Rosso M, Reich DS, Ineichen BV.
    "Which experimental factors govern successful animal-to-human translation in
     multiple sclerosis drug development? A systematic review and meta-analysis."
    eBioMedicine 2024;110:105434. PMID 39515028, PMC11582441, CC BY.

497 animal studies, **15 approved and 11 failed DMTs** — failed on efficacy *or
safety* — roughly 30,000 animals, 274 in the meta-analysis. 86% used EAE.

    "There was no association between animal study outcomes or testing DMTs
     under varied conditions (e.g., different laboratories or models) and
     successful approval."

And the literature is contaminated by hindsight: **91% of the animal studies
were published after the first-in-MS trial, and 91% after regulatory approval.**
So a "preclinical signal" for an approved drug is usually not preclinical.

**A distinction the repo must keep, because this does NOT retract its existing
citations.** Berg et al. measure whether an EAE *efficacy outcome* — did the drug
lower the EAE score — predicts approval. `profiles.py` cites EAE work for
something different: *mechanism*, that TNF-deficient mice develop severe EAE
(Liu 1998) and that TNFR2 is required for remyelination (Arnett 2001). A
mechanism-level finding about what a target does is not an efficacy screen. The
lenercept flag stands; what dies is the idea of ranking candidates by their EAE
effect size.

**Identified but NOT yet sourced — target expression on regulatory cells.** The
remaining idea with a testable prediction: harm risk tracks whether a drug's
target is expressed on regulatory or reparative cells rather than only on
effectors. It would flag daclizumab hardest, because its target IS the canonical
Treg marker (CD25/IL2RA) — and daclizumab was withdrawn for fatal encephalitis.
It would also flag lenercept via TNFR2. Two of the harm cases, from expression
alone, with no trial outcome read.

Blocked on data, not on the idea: the Human Protein Atlas serves per-immune-cell
RNA on its website but its `search_download` API returned nothing for every
immune-cell column tried, and the Kang 2018 dataset already in this repo has no
Treg label (its CD4 T cells are not subset).

**Candidate source identified, not yet retrieved** (Europe PMC was returning 503
at the time):

    "Transcriptomic profiling of human effector and regulatory T cell subsets
     identifies predictive population signatures." PMC8085975.

It profiles effector AND regulatory subsets in the same experiment, which is the
contrast the channel needs — a target's Treg:effector expression ratio, not its
absolute level. Check its GEO accession and supplementary tables first; if they
carry per-gene values the channel is buildable for the targets in this arm set
(IL2RA, TNFRSF1B, MS4A1, CD52, ITGA4, S1PR1, TNFRSF13B).

Two further notes for whoever picks this up. The claude.ai PubMed MCP server is
installed but **not authenticated**, so it exposes only its OAuth tools —
everything cited in this repo today came through the Europe PMC REST API and
Firecrawl's research index instead, neither of which needs auth. And the
prediction to test FIRST is daclizumab: its target is the canonical Treg marker,
and it was withdrawn for fatal encephalitis while REDUCING relapses 45%. A
channel that flags it is seeing something neither MRI nor ARR can.

#### CLOSED 2026-09-20 — three of the five "blockers" were not blocked

1. **MRI extraction: 11 of 12 arms.** AFFIRM's table is in the NEJM paper and
   PRISMS has a companion MRI paper (Li, *Ann Neurol* 1999). "Predates the
   registry" was a reason to look elsewhere, not to stop. Only alemtuzumab
   remains, and for a real reason: CARE-MS I reports the PROPORTION of patients
   lesion-free and a lesion VOLUME change, neither convertible to a count ratio.

2. **The harm channel passes an out-of-sample check.** Of the four drugs the
   binary test calls negatives, alemtuzumab is the only one with a major
   secondary-autoimmunity signature (30-48% of patients) and the only one the
   channel puts above 1.0. Across all six targets the ordering tracks
   regulatory-failure toxicity monotonically. Natalizumab ranks fifth of six
   and that is correct — PML is over-suppression, not lost regulation.

3. **The Sormani intercept is empirically supported.** Refitting on this repo's
   own arms, restricted to comparable contrasts (placebo-controlled, T2
   lesion-COUNT metric), six arms give **slope 0.468, intercept -0.037, R^2
   0.347** against Sormani's published slope of 0.52. The slope agrees and the
   intercept is 1.0 within noise.

   **The restriction is the finding.** Without it, across all 11 arms, the
   relationship collapses to slope 0.081, R^2 0.077. Mixing drug-vs-drug
   contrasts with drug-vs-placebo ones destroys it. Anyone refitting must
   restrict first.

#### What is genuinely left, and it is one thing

**A model whose damage is not peak-driven.** Everything still open reduces to
this. Vélez makes killing effectors worse (damage goes as `(E/a)^2`, so it is
set by peak excursions, and removing effectors releases the proliferation brake
into a larger excursion). That single property produces:

  - five arms unscreenable at any potency, confirmed independently by
    `backtest/potency.py` returning OUT OF RANGE on real trial MRI data;
  - the grounded direction gate at 5/13 against the ABM path's 9/14;
  - LOMO at 45.9pp against a 12.3pp null, dominated by the gamma_E fold;
  - a screen in which nearly every survivor is `gamma_R-`.

**THE FIX GRID IS EXHAUSTED. Six structural variants, every one fails.** Damage
change at depletion, median over 24 infection histories at 730 days:

| variant | at 1.5x | at 2x | at 3x |
|---|---|---|---|
| published (`gamma_E`, n=2, no cap) | +1564% | +5635% | +275382% |
| `gamma_E`, linear damage (n=1) | +107% | +332% | +1719% |
| `gamma_E`, carrying capacity K=50k | +109% | +198% | +362% |
| `gamma_E`, linear + capacity | **+18%** | +36% | +64% |
| additive `-dep*E`, no cap | +281% (dep 0.05) | +1444% | +4020% |
| additive `-dep*E` + capacity | +33% (dep 0.05) | +90% | +148% |

The best case still has the wrong sign. A carrying capacity does two real things
— it stops divergence, and it makes median effector burden FALL with depletion
(1088 -> 982) — and damage still rises, because the stochastic excursions get
worse as the damping weakens.

**Why, and why it is not patchable.** Effectors recruit their own regulators
through `E^h/(ke^h+E^h)`. Remove effectors and `R` falls, the proliferation
brake releases, and the system re-equilibrates with weaker damping and larger
excursions. That loop is not an implementation detail — it is the paper's entire
contribution. You cannot remove it and still have the model.

**So the criterion stated earlier was not sharp enough.** "Does increasing
effector death reduce effector BURDEN" is passed by the capacity variant and is
still not enough. The test a replacement must pass is:

    does increasing effector death reduce DAMAGE, at every potency?

**AND THE DEFECT IS FIXABLE — found by sweeping the capacity rather than testing
one value of it.** The table above used K=50000 against a baseline effector level
of ~1000, where the cap only clips excursions. Sweep it down and the sign flips:

    K = 50000   +109%  +198%  +362%      (damage at gamma_E x1.5 / x2 / x3)
    K = 10000    +30%   +44%   +82%
    K =  3000     +4%    +0%    -4%
    K =  2000     -3%    -7%   -14%      <- depletion finally helps

The cap must BIND AT THE OPERATING POINT. That is the precise property: **effector
growth bounded by something other than the regulatory population, with the bound
active where the system actually sits.** Vélez bounds effectors only by Tregs, so
attacking effectors attacks the bound.

Available as `simulate(carrying_capacity=K)`, an EXTENSION, off by default so
every published-value test still runs on the pure transcription.

**It is not a repair.** At K=2000, three-fold depletion buys 14% less damage
against trials reporting 55-68% relapse reductions. Right sign, wrong magnitude
by a factor of four. This identifies what a replacement model needs; it does not
make this one usable for the depleting class.

**And the repo's own discarded toy passes it.** `bricks/qsp.py` has logistic
growth `r_CA*C*A*(1-A)`, so effector control does not depend on effectors
recruiting their own regulators. Damage falls monotonically with `treat`:
0.906 -> 0.383 -> 0.276 -> 0.218 -> 0.181. That is an uncomfortable result and
it is recorded rather than buried: on this one axis the invented toy is right
where the transcribed published model is wrong.

What a replacement needs, precisely: **effector growth bounded by something
OTHER than the regulatory population.** Vélez bounds it only by Tregs, so
attacking effectors attacks the bound.

Candidate checked and rejected: PMC13171755 (QSP of B-cell immune response in
mouse) — 247 ODE mentions and one "damage". It models B-cell dynamics with no
tissue-damage readout, so it cannot answer the test and is not a replacement.

#### Still open

- ~~**MRI extraction for eight arms** — PRISMS, CONFIRM, AFFIRM, FREEDOMS, TEMSO,
  DEFINE, CARE-MS I, OPTIMUM.~~ **STRUCK 2026-09-20: `backtest/potency.py`'s
  `PENDING_EXTRACTION` is now a one-tuple — alemtuzumab — and for a stated
  reason (CARE-MS I reports a lesion-free PROPORTION and a volume change,
  neither convertible to a count ratio). Eleven of twelve were in the journals.**
- ~~**a non-MRI harm channel.** EAE is ruled out empirically (above). The live
  candidate is target expression on regulatory cells, blocked on Treg-resolved
  expression data rather than on the idea.~~ **STRUCK 2026-09-20: built
  (`bricks/harm_channel.py`, 53f640d) on Human Protein Atlas immune-cell nTPM,
  and then BOUNDED — a pre-registered enlargement (docs/HARM_CHANNEL_PREREG.md)
  found 8 of 9 candidate targets unscoreable by a Treg:effector-T ratio. It
  explains daclizumab and lenercept and will not become a general filter.**
- **The depleting class has no home at all** — not a mapping problem, a model
  property (above). Options, none of them cheap: accept the blind spot and scope
  the screen to proliferation/regulation mechanisms; or find a model whose damage
  is not peak-driven. An additive loss term has been tried and does not work.
- **the screen** (`screen/`) — **half built 2026-09-20 (d7a5643).**
  `screen/kill_filter.py` runs four filters that can only ever say "doomed", and
  `rank_candidates()` raises rather than returning a sorted list. Still open is
  the generate-and-report half: enumerate the intervention-point space, screen
  it, and publish the survivors WITH their kill reasons. Ranking stays gated on
  LOMO beating its null (live: 45.9pp vs a 12.3pp null, re-measured 2026-09-20 —
  NOT the 71.0/14.4 recorded above, which was a four-fold table).
- **wiring.** `qsp_traj` is written but nothing consumes it; the clinical gate
  runs through `abm_damage`. **The port does not move the gate on its own** —
  readout must read `qsp_damage` before any of this reaches the score.
  **STRUCK 2026-09-25: refuted in code on 2026-09-20 (5e9ee90) and never struck
  here. The readout clips damage to [0,1]; the port's damage is unbounded, 78%
  of untreated runs sit above the clip (median 1.54, mean 11.0, max 248 over 32
  runs at 730 d), so wiring it would peg the cohort at the ceiling. The scored
  route is the arm-level Sormani ratio, which `clinical_velez` already uses.**

#### THE HARM CHANNEL WAS ENLARGED UNDER PRE-REGISTRATION (2026-09-20)

The channel's whole statistical claim was 1/C(6,2) = 0.067 on six targets, on a
hypothesis formed after knowing which drugs harmed. Enlarging the scored set is
the only cheap way to make that test mean anything, and it is worthless unless
the set is fixed first — so `docs/HARM_CHANNEL_PREREG.md` was written and
committed (308a56f) with the inclusion rule, the outcome labels, the statistic
and its falsifier, **before any nTPM for the new genes was fetched**. Nine
targets from published MS randomized trials, every citation verified against
Europe PMC that day rather than recalled.

**Eight of the nine fell under the effector floor.** Only DHODH clears it, at
0.76x, below TNFRSF1B — so the pre-declared falsifier did not fire and the two
harm cases stay 1st and 2nd. N goes 6 -> 7 and the exact test goes to
1/21 = 0.048.

**That number needs its provenance attached every time it is quoted.** An
earlier draft of `bricks/harm_channel.py` claimed 1/21 = 0.048 by miscounting
the arm set as seven when it is six, and a test pins the arm-set ranking at six
so it cannot recur. That test is untouched; `ranking()` still returns six and
`ranking_enlarged()` returns seven. Same arithmetic, different provenance.
Nominally crossing 0.05 by adding one target to the denominator licenses
nothing: the kill filter's REGULATORY LIABILITY flag stays soft, no veto.

**The result is the exclusion, not the p-value.** Nearly every MS drug target
not already scored is a B-cell gene (CD19, BTK, CD80/86, TACI), a secreted
ligand (IL17A, IL12B), or not an immune gene at all (LINGO1, HCAR2). A
Treg:effector-T ratio cannot speak to any of them. The channel's applicability
domain is small and cannot be grown by looking harder — and on the harm side
the MS record has exactly two regulatory-failure cases with an independent
immune target. That is the ceiling on this idea, stated so nobody plans a
general safety filter around it.

All 136 transcribed HPA values, the 64 old and the 72 new, verify against the
live download: `PYTHONPATH=. python scripts/derive_harm_channel.py --check`.

#### THE CRITICAL PATH WAS A MIRAGE — MARTINEZ-PASAMAR 2013 IS THE SAME MODEL (2026-09-20)

The entry above put the Martinez-Pasamar port on the critical path, "not optional
and not later", because six of thirteen arms collapse onto `gamma_E` and the port
was believed to add antigen-specific subpopulations and microglia. **It does not.
Checked at the source before spending the days.**

Full text via the Europe PMC REST API (`PMC3651362/fullTextXML`) and the
supplement via `PMC3651362/supplementaryFiles` (13 files, 647 KB; the model
material is `1752-0509-7-34-S1.docx`). From the Methods, verbatim:

    "The model is based in 4 differential equations describing the dynamics of
     antigen specific resting Teff (1), resting Treg (2), activated Teff (3),
     and activated Treg (4)."

    (1) dEr/dt = I_E - Er*delta - Er*beta + E*eta
    (2) dRr/dt = I_R - Rr*delta - Rr*beta + R*eta
    (3) dE/dt  = Er*delta - E*eta + E*alpha_E*(kR^h/(kR^h+R^h)) - E*gamma_E*(R^h/(kR^h+R^h))
    (4) dR/dt  = Rr*delta - R*eta + R*alpha_R*(E^h/(kE^h+E^h)) - R*gamma_E

Compare `bricks/qsp_velez.py:58-61`. **Identical in form, variable for variable.**
Table S1 lists the same parameter set — delta, beta, eta, alpha_E, alpha_R,
gamma_E, gamma_R, k_E, k_R, h — in two columns, mouse and human, and the human
column is Vélez. The paper says so itself: "a model of active Teff-Treg
cross-regulation developed previously to describe the dynamics of T-cells in
humans [10], updating here the model parameters to reproduce experimental data
from EAE studies in mice."

**There is no B-cell state variable and no microglia state variable.** "Microglia"
occurs 53 times and "B-cell" 47, all of it flow cytometry in the experimental
sections. The anti-CD20 result is an EAE experiment; the modelling contribution
is that B-cell depletion can be *represented* as a change to existing T-cell
parameters, not as a compartment. What the port would add over what is already
in the repo: mouse parameter values, and a noise-driven pulse-train input for
naive cells in place of the deterministic influx.

**So it separates nothing.** Anti-CD20, S1P modulation and alpha-4 blockade
still land on the same dials, because there is no state in this model that
distinguishes them. A port would have delivered a second copy of a model the
repo already has, after two or three days.

**What the separation actually needs, stated so the next attempt starts from
here:** a model with a compartment the lumped mechanisms differ in — periphery
vs CNS with an explicit transition (natalizumab blocks the transition,
fingolimod and ponesimod block egress from lymph node, depleters remove cells
from the pool), or an explicit B-cell population (anti-CD20). Neither Vélez 2011
nor Martinez-Pasamar 2013 has one. `bricks/brain_pbpk.py` (Verscheijden 2019)
has compartments but carries drug concentration, not cell populations, so it is
not the missing piece either. **Nothing in this repo's shortlist currently
supplies it, and that is the honest state of the mechanism-separation problem.**

One thing the check leaves intact: `profiles.py` takes ocrelizumab's `ke = 0.85`
from this paper's EAE fit. That is a parameter read, not a port, and it stands.

#### THE CAPACITY EXTENSION DOES NOT LIFT THE GATE (2026-09-20)

`bricks/qsp_velez.py` found that a carrying capacity binding at the operating
point flips the sign of depletion (db1dd2c). That was measured on the `gamma_E`
dial in isolation and was never fed into the gate it was supposed to move.
`backtest/lomo_capacity.py` now does, with **K and potency both fitted inside
each fold on the training arms only** — choosing K by looking at the headline
would be fitting the test set, and the uncapped model stays in the search so the
extension has to earn its place.

    out-of-sample MAE 45.6pp   predict-the-mean null 12.3pp   (uncapped: 45.9pp)

**It does not beat the null, and it does not move.** Fixed-K diagnostics, which
are NOT out-of-sample in K and must not be quoted as the gate: K=None 45.9,
K=10000 49.5, K=3000 58.1, K=2000 50.6, K=1500 45.8, K=1000 43.4. The best of
them is still 3.5x the null.

**What the capacity bought, exactly.** In the `gamma_E` fold the sign is no
longer wrong — it predicts **-3.6% where the five trials report -30% to -68%**.
The fold error moves 53.2pp -> 49.6pp. So the defect was never only the sign:
removing effectors in this model cannot produce a LARGE damage reduction,
because damage is set by peak excursions and a capacity only damps them. This is
the factor-of-four gap already recorded in `qsp_velez.py`, now measured at the
gate instead of on one dial.

Worth naming: the `alpha_R` fold gets worse (108.8pp). With a binding cap,
lowering `alpha_R` predicts **+63.8% where DECIDE reports -45%** for daclizumab.
One arm, so read it beside the others.

**Consequence: the screen stays kill-only, and this was the cheapest route to
ranking.** `screen.rank_candidates()` now cites both numbers in its refusal. The
remaining routes are all expensive and none is on this repo's shortlist: a model
whose damage is not peak-driven, or a representation with a compartment the
lumped mechanisms differ in — and Martinez-Pasamar 2013, the candidate the plan
named for that, turns out to be the same four ODEs (see above).

Tables are cached per capacity in `results/mechanism_curve_K*.json`; the
transcription's own curve is never overwritten. Full rows in
`results/lomo_capacity.json`.

#### THE SCREEN RAN — 98 CANDIDATES, 38 SURVIVORS, NO RANKING (2026-09-20)

`screen/report.py` enumerates every combination of up to two of the model's seven
named intervention points, both directions, at potency 0.5, and runs the four
kill filters over all 98. Artifacts: `results/screen.json`, `docs/SCREEN_RESULTS.md`.

    OUT_OF_REGIME   32      the model runs away; damage undefined
    UNREACHABLE     17      no benefit at any probed potency
    DEGENERATE      11      same dials as a drug that already exists
    SURVIVED        38

**Survivors are listed alphabetically and that is enforced in a test**, with a
fixture whose alphabetical order deliberately disagrees with its damage order.
Every survivor is by construction a mechanism pattern no existing drug occupies,
which is exactly the case `backtest/lomo.py` scores at 45.9pp against a 12.3pp
null — so an ordered list here would be a sorted list of noise.

**One earlier claim corrected by the run.** §8.4 above says "a screen in which
nearly every survivor is `gamma_R-`". Measured: **12 of 38**, about a third. The
dial is over-represented, not dominant.

Half the candidate space (49 of 98) dies on OUT_OF_REGIME or UNREACHABLE, both
of which are properties of the model rather than of the candidates — a screen
whose kill reasons are mostly "my simulator broke" is describing itself. That is
the honest reading and it is why surviving is not passing.

#### THE POTENCY CROSS-CHECK WAS QUOTED IN THE WRONG UNITS (2026-09-20)

`backtest/potency.py` called ocrelizumab's two independent magnitude estimates
"a factor of two apart". That is `ke` 0.85/0.40 = 2.1x in MULTIPLIER space. The
model fits **potency `s = 1 - ke`**, where the same two points are 0.60/0.15 =
**4.0x** apart. Neither ratio is privileged; a ratio between parameterisations is
an artifact of which one gets printed. Found by the parallel session working the
decision-rule lane, verified here before changing anything.

Quoted in the units the gate actually scores, through the cached response table:

    EAE fit  (s = 0.15)   -28.5%
    MRI fit  (s = 0.60)   -73.8%          OPERA I reported -46.0%

**45.3pp apart, straddling the trial's own number.** Both independent sources
are wrong and they are wrong in opposite directions. That is the honest size of
the uncertainty in the potency layer, and it is the same class of defect as the
hardcoded "9/14, 1/9" in `clinical_velez.py`: not stale, just never restated in
comparable units.

#### A MODEL WITH THE MISSING COMPARTMENT EXISTS — AND IT STILL WOULD NOT SEPARATE THE ARMS (2026-09-20)

The entry above left the mechanism-separation problem as "nothing on this repo's
shortlist supplies a compartment the lumped mechanisms differ in". Searched
properly rather than assumed — Europe PMC's relevance ranking is useless for a
conceptual query, OpenAlex found it immediately:

    Gazola GM, de Oliveira JVC, de Paula MAM, Quintela BM, Lobosco M.
    "Fingolimod and Neuroinflammation in MS: Representing CD8+ T-Cell Dynamics
     Through Mathematical Modeling and Clinical Evidence."
    Sclerosis 2025;3(4):38. doi:10.3390/sclerosis3040038. Open access, gold.

**It has the compartment.** Twelve state variables in two coupled subsystems:
six PDEs in brain tissue (microglia, CD8+ T, oligodendrocyte damage, antibody,
immature and activated dendritic cells) and six ODEs in the lymph node (DC,
CD8+ T, CD4+ T, B cells, plasma cells, antibody). Trafficking is explicit and
bidirectional, through `gamma_T`, `gamma_D` and `gamma_At`, scaled by the
vascular and perivascular interface areas. Equations and structure are all in
the paper; `epsilon` is fitted per patient group (0.81-0.90) against Song et
al.'s 23-patient, 360-day CD8+ cohort.

**And it would still not separate natalizumab from fingolimod. Read Equation (4).**
The tissue influx term is

    (1 - epsilon) * gamma_T * theta_BV * (T_CL(t) - T)

**One transmigration coefficient, one blockade parameter, applied to the flux
between the two compartments.** Fingolimod sequesters lymphocytes in the lymph
node by blocking EGRESS; natalizumab blocks ENTRY across the endothelium. Those
are different ends of the same flux and this model gives them the same dial, so
porting it buys the repo the same degeneracy it already has, in twelve variables
instead of four.

Worth noting precisely, because it is a criticism of the paper and not of the
idea: `epsilon` appears ONLY in the tissue influx (4), not in the lymph node
efflux (13). So as written it models blockade of entry INTO brain tissue, which
is mechanistically the natalizumab story, while the paper calls it fingolimod.
The place where the two drugs could be separated exists in this structure — put
a second parameter on the lymph node efflux term — but that is a model change,
not a port, and it would be unsourced the moment it was made.

**Cost, stated so nobody starts it casually.** A 2D PDE system on a vascular
indicator grid, against this repo's budget of 0.05s per 730-day run and a
128-seed cohort per cell of a 20-potency grid. The authors themselves could not
run Monte Carlo UQ on it and fell back to a sparse polynomial chaos surrogate.
It also has **no Treg population at all**, so the cross-regulation the entire
harm channel rests on — and the lenercept and daclizumab findings with it — has
no home in it.

**Conclusion: found, assessed, NOT recommended.** The compartment alone is not
the missing piece; the missing piece is two SEPARATE trafficking rates, egress
and entry, each with its own source. No published MS model found tonight has
that. Second candidate not yet assessed, named here so it is not lost: Pernice
et al. 2020, BMC Bioinformatics, doi:10.1186/s12859-020-03823-9, the Epimod
framework — stochastic, models blood-brain barrier integrity explicitly, and
was run for **daclizumab**, which is this repo's sharpest harm case.

#### THE BEST SUCCESSOR MODEL FOUND SO FAR — PERNICE 2020, AND HOW IT MAY LEGALLY BE USED (2026-09-20)

Named as the second candidate above; assessed now, because it is much closer to
what this repo needs than either Martinez-Pasamar or Gazola.

    Pernice S, Follia L, Maglione A, Pennisi M, Pappalardo F, Novelli F,
    Clerico M, Beccuti M, Cordero F, Rolla S.
    "Computational modeling of the immune response in multiple sclerosis using
     epimod framework." BMC Bioinformatics 2020;21(Suppl 17):550.
    doi:10.1186/s12859-020-03823-9. Open access, CC BY.

**What it has that Vélez does not.** 26 places and 55 transitions as an Extended
Stochastic Symmetric Net, from which **a system of 26 ODEs with 20 calibrated
parameters is derived** (their words). Two compartments — peripheral lymph
node/blood vessel, and CNS — interacting **through an explicit BBB place**, with
transitions suffixed `_out` and `_in` accordingly. Crucially it keeps the
Teff/Treg cross-regulation this repo's whole harm story rests on
(`TregKillsTeff_out`, `TregKillsTeff_in`), and adds NK cells, IL-17, IFN-gamma,
IL-10, and oligodendrocytes with five discrete myelination levels plus
remyelination.

**What it would actually separate, and what it would not.** Transitions are
distinct for peripheral killing, duplication, and **BBB passage**
(`Teff_pass_BBB`, `Treg_pass_BBB`). So:

  - **depletion vs trafficking blockade IS separable** — ocrelizumab and
    alemtuzumab act on peripheral Teff, natalizumab acts on `Teff_pass_BBB`.
    That is two of the arms currently collapsed onto `gamma_E`.
  - **fingolimod is STILL not separable from natalizumab.** Fingolimod's
    mechanism is retention in the lymph node, and this model lumps lymph node
    and blood vessel into ONE compartment, so there is no egress step to block.
  - **daclizumab becomes expressible as the harm case it is.** The paper gives
    it two independent parameters, `DACkillTeff` and `DACkillTreg`. That is this
    repo's harm hypothesis — a drug that strips regulation while suppressing
    effectors — as two dials rather than as one flag.

**THE LICENSING, CHECKED BEFORE ANY PORT IS PLANNED.** The analysis code and the
net files are at `github.com/qBioTurin/Multiple-Sclerosis`, and that repository
**has no licence file of any kind**, which means all rights reserved — it cannot
be vendored, copied, or adapted here. The `epimod` framework
(`github.com/qBioTurin/epimod`) declares **GPL (>=2)** in its DESCRIPTION, which
is a copyleft this repo is not going to take on for a brick.

**The route that is open is the one Vélez came in by:** transcribe from the
PAPER, which is CC BY — the net description in the Results, and the parameter
table in Additional file 1, Table S1. Cite it, tune nothing, and vendor no code.
That is a real port with a real licence story, and it is the first candidate all
night that has one.

**Not started tonight, and not to be started casually:** 26 ODEs with 20
parameters that were calibrated on 16 subjects' cytokine counts, against this
repo's 4 states and 6 dials. Scope it before committing to it.

#### THE SCREEN UNDER THE CAPACITY EXTENSION — THE VERDICTS ARE NOT STABLE (2026-09-20)

The capacity extension bought the ranking gate nothing (45.6pp against 45.9pp,
above). But a third of the candidate space was dying on OUT_OF_REGIME, which is
the model admitting failure rather than judging a candidate, and stopping
divergence is exactly what a binding capacity does. So the screen was re-run
under it. `results/screen_K2000.json`, `docs/SCREEN_RESULTS_K2000.md`, written
to separate files because it is a different model.

| verdict | transcription | K = 2000 |
|---|---|---|
| OUT_OF_REGIME | 32 | **0** |
| UNREACHABLE | 17 | **60** |
| DEGENERATE | 11 | 7 |
| SURVIVED | 38 | **31** |

**Every one of the 32 divergent candidates became UNREACHABLE. Not one of them
became a survivor.** That is worth having: "the model ran away" is an admission,
"no potency on this dial helps" is a verdict, and the extension converts the
whole of the former into the latter without turning up a single hidden
candidate. It also makes the survivor list less silly — the count claiming a
bigger effect than natalizumab drops from 32 of 38 to 20 of 31.

**And here is the finding, which is about the screen and not about the capacity.**

    survivors agreeing across the two models:   29 of 40
    survivors only under the transcription:      9
    survivors only under K = 2000:               2

**Two models that the gate cannot tell apart — 45.9pp versus 45.6pp, inside the
noise floor — disagree about 11 of 40 candidates, roughly 28%.** Nine patterns
that survive the published model are killed by the extension, including five of
the `gamma_E-` pairs, and two that the published model never scored at all now
survive.

So the survivor list is not a stable object. It carries model-selection
uncertainty that **no measurement in this repo constrains**, because the one
measurement that could discriminate the two models scores them as equal. A
candidate's survival is contingent on a modelling choice made for numerical
reasons, and anyone quoting the survivor list must quote which model produced
it. Both files ship, neither replaces the other, and `screen/report.py`
banners the extension run so the two cannot be confused.

The honest summary of the screen as a whole, after tonight: it kills reliably
(60 of 98 under either model, for reasons that are properties of the model),
it does not rank, and its remainder moves when the model does.

---

## 2026-09-20 — the accept/reject device, and three failures that are simultaneous

Built on branch `gate/decision-rule`. New package `gate/`, new docs
`docs/DECISION_GATE.md` (the one page for a reader who will not run the code)
and `docs/RECOVERABILITY.md` (the table, with its regenerating command).

**The device.** `gate.decide(candidate)` returns KILL, ABSTAIN or PASS. Over the
14 single-dial candidates: 12 KILL, 2 ABSTAIN, **0 PASS**, and no input can
return PASS today. Not hard-coded — PASS requires an `EvidenceCertificate` that
re-measures the out-of-sample scorers and asks whether they beat predict-the-mean.
Both lose, with paired-bootstrap CIs entirely above zero: LOO 28.3pp vs 11.5pp,
CI [+5.7, +28.5]; LOMO 45.9pp vs 12.3pp, CI [+20.1, +43.1]. The criterion is
frozen and dated in `gate/criterion.py` and every clause restates a bar this
repo already held. `tests/test_gate.py` drives the PASS branch with a synthetic
passing certificate, so the day a real measurement inverts, the device changes
its answer with nobody editing it.

**Three failures, and they are simultaneous rather than sequential.** This is
the part that changes what to do next.

| | measurement | verdict |
|---|---|---|
| (1) reachability | 6 of 12 arms fit at the grid edge; oracle 51.9pp vs 12.4pp null | no potency produces the observed direction |
| (2) recoverability | MRI-fitted potencies score 21.4pp vs 10.6pp null, biased high every time, median 1.55x | the independent channel cannot find the potency needed |
| (3) expressiveness | shared-potency LOO restricted to the *reachable* arms: 20.3pp vs 11.7pp (13.4 vs 11.4 dropping the single-arm fold) | one global potency cannot separate same-dial drugs |

Fixing any one leaves the gate red. (3) is the one neither of us had tested: we
had both concluded the model's form was the blocker, and it is not *the* blocker
— the scorer fails on the arms the model can already reach. So **blockers (4) and
(6) are closed routes rather than deferred ones**, and the remaining work is a
model form AND a working potency source, not a choice between them.

**How much is there to win.** A *perfect* dial-level model — predict each arm by
its dial group, no fitting, no simulation — scores 6.6pp vs a 10.6pp null in
sample. Scored the way every other scorer here is scored (each arm from the
other arms in its group, with the null drawn from the same ten arms) the
headroom is **0.8pp**: 10.9pp vs 11.7pp. Against a measured 45.9pp. A replacement model does not need to be better; it needs to be
within about a point of perfect.

That is a property of the **arm set**, not the model. 12 quantified arms in 3
multi-member dial groups is a thin exam, and widening it — especially more arms
*per dial* — is far cheaper than any port and is the first thing to do before
concluding a future model has failed.

**One correction to a number this plan and `backtest/potency.py` both carry.**
The ocrelizumab cross-check is quoted at 2.1x, which is `ke` 0.85 vs 0.40 in
*multiplier* space. The model fits *potency*, `s = 1 - ke`, and in that space the
two independent estimates are **4.0x apart**. The potency that reproduces OPERA's
ARR (0.30) sits halfway between them on a log scale — 2.00x above the EAE fit,
2.00x below the MRI fit. Both independent sources miss, in opposite directions.
Full working in `docs/RECOVERABILITY.md`, which also records why that page counts
six unreachable arms where `backtest/potency.py` counts five (alemtuzumab has no
MRI number to be out of range with).

Reproduce all of it:

```
PYTHONPATH=. python3 -m gate.device      # verdicts and the certificate
PYTHONPATH=. python3 -m gate.ceiling     # oracle ceiling and recoverability
PYTHONPATH=. python3 -m gate.headroom    # which failure binds, and the 1.0pp prize
```

#### WIDENING THE EXAM BARELY MOVED THE PRIZE (2026-09-21)

`scripts/dial_ceiling.py` had just established that the out-of-sample headroom
over a predict-the-mean null was **0.8pp**, and that the headroom is a property
of the ARMS rather than of the model — twelve quantified arms across three
multi-member dial groups is a thin exam. That made widening the arm set the
cheapest lever available, far cheaper than porting a richer model. The parallel
session researched eight candidate arms (f2ca0d4); six were wired here.

    quantified:      ublituximab -> ke, IFN-beta-1b -> alpha_E, ozanimod -> gamma_E
    direction-only:  rituximab -> ke, ustekinumab -> alpha_E, abatacept -> delta

    arms 17 -> 23, quantified 12 -> 15

**Out-of-sample headroom went 0.8pp -> 0.9pp.** Measured, not estimated, and
confirmed by two independent implementations (`scripts/dial_ceiling.py` and
`gate/headroom.py`) that agree to 1e-9 on all three variants.

**Why it barely moved, which is the useful part.** Two of the three quantified
arms landed in groups that were already the largest — `gamma_E` went 5 to 6 and
`alpha_E` 3 to 4 — and the out-of-sample variant can only score groups with at
least two members. The groups it is starved of are the singletons: `alpha_R`
(daclizumab), `alpha_R|delta` (glatiramer) and now `delta` (abatacept, which
carries no magnitude). **The lever is real and the per-arm return is far smaller
than "widen the exam" implied.** What it needs is arms in the SMALL groups, and
a quantified second arm on `alpha_R` is worth more than five more anti-CD20s.

#### THREE ASSIGNMENTS DELIBERATELY NOT MADE — and one of them was tempting

Recorded in full under "Negative results from dial assignment" in
`docs/TRIAL_ANCHORS.md`, because the tempting one is the shape of mistake this
repo keeps finding.

**Abatacept was nearly given an `alpha_R` axis.** CTLA4-Ig blocks CD28
costimulation, and regulatory T cells depend on CD28, so impairing them is
plausible. Asserting it would have put abatacept on `alpha_R|delta` beside
glatiramer — **turning a singleton into a scored pair, which is directly worth
headroom in the measurement being reported in the same commit.** That is exactly
when to go and look for the source. Europe PMC on 2026-09-21 has CTLA4-Ig/Treg
evidence in LRBA deficiency, transplantation and autoimmune haemolytic anaemia,
and none in MS. No source, no axis; abatacept sits on `delta` alone.

Laquinimod (AhR agonism) and secukinumab (IL-17) got no dial at all — neither
mechanism corresponds to a rate this model has.

#### A CITATION POINTED AT THE WRONG PAPER (2026-09-21)

The parallel session reported that four of five PMIDs it recalled from memory
resolved to unrelated work, and checked its own by live query instead. Running
that check over the EXISTING anchor table found one:

    cladribine, CLARITY 2010, cited as PMID 20089950
      -> "Signaling by the high-affinity HDL receptor scavenger receptor B type I."
    the real CLARITY paper is PMID 20089960

One digit. The row's own numbers were checked against the correct abstract and
are right — ARR 0.14 vs 0.33 at 96 weeks — so the data was sound and only the
pointer was wrong. **That is worse than a wrong number, because it survives
every sanity check a reader applies to the number itself.**
`scripts/verify_anchors.py` now resolves every identifier in the table and
flags any whose title is off topic. It is not in the test suite: it needs the
network, and a test that fails when the wifi drops teaches people to skip tests.

#### THE EXAM CANNOT BE WIDENED WHERE IT MOST NEEDS IT (2026-09-21)

The widening result above said the lever needs arms in the SMALL dial groups —
`alpha_R`, `alpha_R|delta` — not more arms in the large ones. The parallel
session then went looking for a quantified arm on a regulatory dial and came
back with a null result that is sharper than "nothing found" (e2f8cb9,
`docs/TRIAL_ANCHORS.md`). Searched live 2026-09-21:

| candidate | why it does not qualify |
|---|---|
| low-dose IL-2 | randomised in SLE and Sjögren's; no MS trial with a relapse endpoint |
| Tovaxin / TERMS | randomised, placebo-controlled, n=150 RRMS — **missed its primary**; ARR only post-hoc, against a placebo arm the paper says prior DMT lowered |
| T-cell vaccination | randomised and double-blind but n=17/arm, relapsing-PROGRESSIVE, outcome is proportion relapse-free, not an ARR |
| ATX-MS-1467 | two real 2018 *Neurology* trials, both open-label, one single-arm, MRI endpoints |
| basiliximab | same target as daclizumab; no randomised MS trial found |

**No MS therapy whose primary mechanism is regulatory-T-cell restoration has a
published ARR against a stated comparator.**

**So the singletons are singletons in the arm set because they are nearly
singletons in the FIELD.** The model's regulatory dials sit on the region of MS
pharmacology with almost no successful randomised evidence. That is a boundary
on the cheapest lever this repo has: **the exam can only be widened where the
trials already are — `alpha_E`, `gamma_E`, `ke` — and those are precisely the
groups that are already the largest.** It explains why six arms moved the prize
0.8pp to 0.9pp, and it predicts the next six will do the same.

The uncomfortable corollary, worth stating because it is not about this repo:
the dials this model can least afford to leave unmeasured are the ones the
field has the least evidence on.

#### TOVAXIN WIRED — the honest dial is the unhelpful one

TERMS is exactly the kind of arm the set is starved of: a real randomised
NEGATIVE result in RRMS. Pairing it with daclizumab on `alpha_R` would have
converted a singleton into a scored group and improved the ceiling reported in
the same commit — which is the reason to distrust that reading. Tovaxin is
autologous attenuated myelin-reactive T cells and the response it induces
deletes effector clones, so its defensible dial is **`gamma_E`**, already the
largest group, which buys the out-of-sample variant nothing.

Wired anyway, direction-only, on completeness grounds: **`gamma_E`'s other six
arms all worked**, and a dial whose training data is all successes cannot teach
a model that the dial sometimes does nothing. Arms 23 -> 24, quantified
unchanged at 15.

That is the second time in two hours that the assignment which would have
improved the measurement was the one without the evidence, and the second time
the search, not the judgement, is what settled it.

#### WHAT THE 23-ARM SET DID TO EVERY GATE (2026-09-21, 00:41 run)

Regenerated by `scripts/state_of_build.py`, every figure produced by the
command beside it in one run. The 12-arm column is the 2026-09-20 measurement
recorded earlier in this section.

| gate | 12 quantified arms | 15 quantified arms |
|---|---|---|
| clinical, ABM path | 10/17 direction, 1/12 magnitude | **11/23, 1/15** |
| clinical, grounded stack | 5/16, 2/11 | **6/22, 3/14** |
| leave-one-ARM-out | 28.3pp vs 11.5pp null | **30.0pp vs 11.3pp** |
| leave-one-MECHANISM-out | 45.9pp vs 12.3pp | **45.4pp vs 11.8pp** |
| LOMO + capacity extension | 45.6pp vs 12.3pp | **42.2pp vs 11.8pp** |
| per-drug potency | 11 fitted, 5 out of range | unchanged |
| dial ceiling, out of sample | 10.9 vs 11.7, headroom 0.8pp | **10.7 vs 11.5, headroom 0.9pp** |

**Nothing inverted. Every scorer with a null still loses to it**, and the nulls
themselves fell slightly, because the added arms cluster near the existing mean
rather than widening the outcome spread — which is the same fact the headroom
measures from the other side.

Two rows deserve a note rather than a glance.

**The direction gate got proportionally WORSE, and that was the intent.** 10/17
is 59%, 11/23 is 48%. Three of the six arms added are drugs that a
mechanism-class rule predicts will work and that did NOT: ustekinumab,
abatacept and Tovaxin all suppress by mechanism and all missed. The exam got
harder on purpose; a rule that maps class to direction cannot pass it, and now
the score says so on nine failures instead of seven.

~~**The capacity extension and the transcription have pulled apart.** They scored
45.6 and 45.9 on the 12-arm set — inside the noise floor, which is what made
`screen/`'s 28% verdict disagreement so pointed. On 15 arms they are 42.2 and
45.4, a 3.2pp gap. But the two models are no longer indistinguishable on the
evidence, and if that separation holds as arms are added it eventually becomes
possible to say which one to screen over.~~

**WRONG ON BOTH COUNTS, STRUCK 2026-09-21. The parallel session re-ran it
rather than taking my numbers and found two errors; both reproduced here
before this correction was written.**

**The label was wrong.** 42.2pp is the variant with **K fitted per fold on the
training arms** — an extra free parameter fitted inside each fold, which should
beat a fixed-K model by construction. The screen ran at **fixed K = 2000**, and
that variant scores **48.0pp, WORSE than the transcription's 45.4pp**. Full
diagnostic, re-measured:

    K=None 45.4   K=10000 51.2   K=3000 54.7   K=2000 48.0   K=1500 42.5   K=1000 40.4

So the comparison as I wrote it set a one-parameter model against a
two-parameter one and credited the difference to the wrong variant.

**And neither gap clears the noise floor.** `backtest/lomo.py` states a
bootstrapped floor of roughly 10% on these figures — about 4.5pp here. The
fold-fitted gap is 3.2pp and the fixed-K gap is 2.6pp. By this repo's own rule
that a difference below the floor is not a result, **the evidence has not
started to separate the two models.**

**The conclusion that depended on this is untouched, and is now better
supported than when I overstated it.** The provenance argument is that two
models this repo CANNOT TELL APART disagree about 11 of 40 survivors, and
2.6pp inside a 4.5pp floor is still cannot-tell-apart. Nothing about the
ranking gate changes: every variant above loses to an 11.8pp null by three to
four times.

### 2026-09-25 — exam v2: the widened exam is passable, the model switches itself off

The point-ARR exam discards the eight direction-only arms. Pre-registered
(`docs/EXAM_V2_PREREG.md`, 80a333d) and measured (`backtest/exam_v2.py`,
`results/exam_v2.json`) on all 23 treated arms with interval outcomes:

    S1 interval-LOMO      30.7pp vs 20.0pp null      loses
    S2 oracle             9.1pp vs 20.0pp null       headroom 10.9pp (54%)
    S3 direction          13% / balanced 0.33        no information
    S4 MRI within-dial    tau +0.80, p = 0.024       the channel ranks

Measured first, before the exam was widened: on the 15 quantified arms a
PERFECT mechanism-level model scores 7.3pp against an 11.8pp LOMO null, a
4.5pp headroom, under the model's own ~6pp noise floor at 128 seeds. The
widened exam more than doubles the prize.

**The fitted potency is 0.00 in every fold.** `gamma_E` predicts harm for
six arms that worked, `delta|naive_E` predicts benefit for two arms that
harmed, `delta` predicts −43% for an arm that did nothing, and the shared
potency does least damage at zero. The point exam never saw this because the
arms that expose it are the ones it discards.

**S4 is the first positive out-of-sample-shaped result here.** The trial's
observed lesion ratio orders arms within a dial on 8 of 10 pairs. The
recoverability finding (biased 1.55x high) was about magnitude; rank
survives it. Two same-metric pairs only, so the metric caveat stands.

Three acceptance tests for any replacement model, from the sign table: keep
the sign of depletion, keep the sign of immunogenic challenge, predict
nothing for costimulation block. The gap list is `docs/GATE_GAP_ANALYSIS.md`.
The 1024-seed tables (gap G2) are building detached to
`results/mechanism_curve_n1024.json`; the qsp_damage wiring item above is
struck (refuted 2026-09-20, never recorded here).

Recovered from the 2026-09-21 session record: a direction-only scorer was
proposed at 05:22 as "roughly an hour's work" and never built. S3 is it.

### 2026-09-26 — gap G2 closed: the miss is form, not noise

`scripts/build_curve_n1024.py` rebuilt every response table at 1024 seeds
(1029 tried, 5 skipped for zero untreated damage: 220, 615, 775, 877, 989),
into `results/mechanism_curve_n1024.json` and `results/exam_v2_curve_n1024.json`.
The 128-seed files are untouched. Measured with the same scorers:

    LOMO            n128 45.4 vs 11.8      n1024 45.4 vs 11.8     identical, every fold
    exam v2 S1      n128 30.7 vs 20.0      n1024 30.7 vs 20.0     potency 0.00 in all 8 folds
    curve at s=0.5  alpha_E -33.8 -> -32.6   ke -66.3 -> -69.2   gamma_E +234.9 -> +208.3

The response curves move by about a point between the two cohorts, and the
gate does not move at all. The ~6pp noise floor was real for a hypothetical
near-perfect model; it is irrelevant to this one, which sits 34pp from the
null for structural reasons. G2 is closed and nothing further is bought by
cohort size. G3, the model form, is the whole of what is left.
