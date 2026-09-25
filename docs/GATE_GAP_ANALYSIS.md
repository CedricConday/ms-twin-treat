# Every gap between this build and a green gate

Written 2026-09-25 against `cdd4831`, after reading BUILD_PLAN.md end to end,
every doc in `docs/`, the scorers, and the recoverable session records. Each
gap names the file, what is wrong, and the cheapest measured fix. Ordered by
cost. Numbers marked *measured* were produced today by the command beside them.

## The state, in one table

| | |
|---|---|
| scored path | `bricks/qsp_velez` → `bricks/profiles` → `bricks/sormani` |
| exam | 15 arms with a point ARR; 8 direction-only arms **discarded** |
| LOMO | 45.4pp vs 11.8pp null (`backtest/lomo.py`) |
| LOMO oracle, perfect mechanism-level model | **7.3pp** (*measured*, `docs/EXAM_V2_PREREG.md`) |
| LOMO headroom | **4.5pp** |
| model noise floor at 128 seeds | **≈6pp** (13.7% CV, `scripts/measure_qsp_variance.py`) |
| arms the model cannot reach at any potency | 7 of 15 (`gate/ceiling.py`) |

Read together: the gate is red for three separable reasons, and only one of
them is the model.

## G1. The exam discards the arms that carry the signal — *fixed this commit*

`backtest/lomo.py` and `gate/evidence.py` score only arms with a point ARR.
The eight direction-only arms — three that did nothing, four that harmed, one
that worked without a rate — are the only outcomes that differ *by mechanism*
rather than by trial noise, and two of them share a dial with six arms that all
worked. `docs/EXAM_V2_PREREG.md` fixes interval scoring for all 23 arms and
`backtest/exam_v2.py` measures it. The oracle on that exam (S2) is the number
that says whether ranking is achievable at all.

## G2. The headroom is under the noise floor — *one detached job*

A perfect model has 4.5pp to win; one model prediction carries ≈6pp of
stochastic noise at 128 seeds. So even a near-perfect model fails the current
criterion by luck about as often as it passes. Cohort size is the lever:
bootstrapped CV of the median was 21.2% at n=48 and 13.7% at n=128, so 1024
seeds should land near 5% (≈2pp). Cost: each response table goes from ~8 min
to ~1 h on this box. Run detached, log to a file, rebuild
`results/mechanism_curve*.json` once, and re-run every gate. This also tells us
how much of the 45pp is noise and how much is form, which nothing has
separated yet.

## G3. The model's form: damage is peak-driven — *the port, 3 to 5 days*

`bricks/qsp_velez.py`: damage goes as `(E/a)^2`, effectors recruit their own
regulators, so removing effectors releases the brake and damage *rises*. Seven
of fifteen quantified arms are unreachable at any potency, and the capacity
extension does not fix it (`backtest/lomo_capacity.py`, and see G7). The fix
grid in BUILD_PLAN §8.4 is exhausted; the plan's own conclusion is a model
whose damage is load-driven and that has a compartment the lumped mechanisms
differ in. That model is scoped, licensed and specified in
`docs/PERNICE_PORT_SCOPE.md` (26 ODEs, CC BY paper, ODC level-count readout,
separate BBB passage rates). It is the single largest item and the only one
that can move the reachable count. Acceptance test, already stated there: does
increasing effector loss reduce damage at every potency, and does depletion
separate from transit block.

## G4. Within a dial, every arm gets the same prediction — *measured today (S4)*

The only per-drug input is potency, and the MRI-fitted potency is biased 1.55x
high (`docs/RECOVERABILITY.md`). Bias does not destroy *rank*. S4 in exam v2
asks whether the observed lesion ratio orders arms within a dial group at all.
If it does, potency from MRI is a ranking input even though it is not a
magnitude input, and that is enough for a screen that ranks candidates on the
same dial. If it does not, the potency layer needs another source.

## G5. Two scorers disagree on active comparators — *fixed in exam v2*

`backtest/loo.py` scores an active-comparator arm as the ratio of two simulated
arms; `backtest/lomo.py` predicts every arm against untreated. Five of the
fifteen quantified arms are active-comparator, so the two headlines are not on
the same footing. Exam v2 adjusts (primary) and reports the unadjusted figure
beside it.

## G6. Stale and hard-coded numbers — *fixed this commit where cheap*

Recovered from the 2026-09-20/21 session records against the tree:

- BUILD_PLAN.md:996 still lists "readout must read `qsp_damage`" as open
  wiring. It was refuted in code (5e9ee90): the port's damage is unbounded,
  78% of untreated runs sit above the readout's [0,1] clip, so wiring it would
  peg the cohort at the ceiling. Struck below.
- `backtest/lomo.py:380` prints a noise floor of "roughly 10%"; the measured
  figure in the same file's header is 13.7%. BUILD_PLAN §8.4 (2026-09-21)
  derives "about 4.5pp" from the 10%; at 13.7% it is ≈6pp. Conclusion
  unchanged, figure corrected in the log.
- The "11 of 40 survivors (28%)" disagreement between the transcription and
  the K=2000 screen is a literal in `gate/device.py`, `gate/provenance.py`,
  `tests/test_gate.py` and `docs/DECISION_GATE.md`, computed on the 38/31
  survivor sets. The screens were regenerated on the 23-arm library (37/30)
  and no script recomputes the overlap. Left standing, flagged: it is the
  class of typed number the repo forbids.
- `backtest/potency.py` PENDING_EXTRACTION is four arms (alemtuzumab,
  ublituximab, IFN-beta-1b, ozanimod); BUILD_PLAN.md:973 says one.
- `results/STATE.md` was generated at `f9a13a8`, a SHA that only exists in the
  reflog; regenerate.
- README.md still describes the 17-arm exam.
- `docs/TRIAL_ANCHORS.md` candidate rows propose `alpha_E` for abatacept and
  secukinumab; the wiring put abatacept on `delta` and left secukinumab
  unassigned, and only the negative-results section says so.
- The 12-arm per-fold certificate and capacity rows were overwritten by the
  15-arm run and survive only in the session records. Not restored; the 15-arm
  numbers are the live ones.

## G7. The capacity fit sits on the grid edge — *state it*

In `results/lomo_capacity.json` the fitted K is 1000, the lowest value in the
grid, in four of five folds, and the fitted potency is 0.95, the highest, in
two of five. A fit pinned to its search bounds is a fit that wants to leave
them. `gate/ceiling.py` already refuses to count an edge fit as a fit; the
42.2pp "K fitted per fold" headline should carry the same flag.

## G8. Bricks that nothing scored flows through — *the operator's call*

`spine/run_demo.py` runs seven stages and reports 8 of 15 state keys as
stand-ins. Of the bricks: `qsp.py` (toy, superseded), `barrier.py` (inert: no
arm sets `cns_required`), `grn.py`, `cell_transfer.py`, `cell_scgpt.py`,
`brain_pbpk.py` (stand-ins), `abm.py` (the ABM-path gate, 11/24 direction),
`dream.py`, `baselines.py`. None of them reaches the LOMO figure. The scored
stack is three modules. Two honest options: label every brick SCORED or
DECORATIVE in `bricks/README.md` and keep the demo as a demo, or quarantine the
decorative ones so the spine only contains what the gate reads. Either is
better than a seven-stage pipeline of which three stages count.

## G9. There is no virtual population in the gate — *after G3*

LOMO runs one parameter set over 128 seeds. `bricks/vpop.py` and
`sample_vpop_velez` exist and are not on the scored path, so "digital twin" is
a name the gate does not earn yet. `novainsilico/vpop-calibration` (MIT,
Python, pushed 2026-09-25; verified via the GitHub API today) calibrates
virtual populations with mixed-effects and GP surrogates and would replace the
rejection filter. Worth nothing until the model form is right (G3), because a
calibrated population of a model that cannot reach seven arms is a calibrated
wrong answer.

## G10. One endpoint — *later*

Everything is ARR. Confirmed disability worsening is the endpoint that
separates progression from relapse, and the CLARITY and ADVANCE synthetic
placebo arms (PMC12488035, open on Figshare) carry it. `noemimontobbio/msprog`
(R, CRAN 1.0.0, GitHub licence field unset; verified today) computes CDW and
PIRA events from EDSS visits the way trials define them. Adds a second exam;
does not fix the first.

## The pasted resource list, checked against GitHub today

| repo | licence | last push | fits gap |
|---|---|---|---|
| `noemimontobbio/msprog` | not asserted on GitHub, check CRAN | 2026-09-04 | G10 |
| `FCACollin/rpack_pira` | none (all rights reserved) | 2021-02-20 | G10, cross-check only |
| `Open-Systems-Pharmacology/PK-Sim` | not asserted on GitHub (GPL per project) | 2026-09-24 | brain exposure, not on the scored path |
| `8navid/agent_ms_PMC12844537` | MIT, 0 stars | 2026-06-26 | remyelination ABM; the ABM path is not the scored path |
| `chenlingantelope/MSscRNAseq2019` | MIT | 2021-06-17 | blood vs CSF calibration data for a compartment model (G3) |
| `blakeaw/pysb-pkpd` | BSD-2 | 2025-06-13 | occupancy PK/PD, not needed until per-drug potency has a source |
| `anikaliu/CAMDA-DILI` | GPL-3 | 2021-10-12 | out of scope for the gate |
| `novainsilico/vpop-calibration` | MIT | 2026-09-25 | G9 |
| `openPfizer/vpop-gen` | Apache-2.0, MATLAB | 2018-06-29 | reference only |

None of them moves LOMO. The three that matter for the gate are G1 (done), G2
(a job), G3 (the port).

## Order

1. **G1** — run exam v2; read S2. If the widened exam has no headroom, stop
   pursuing ranking and say so. *(today)*
2. **G2** — 1024-seed tables, detached; re-run both exams. Separates noise
   from form. *(one job, ~6 h of box time)*
3. **G3** — Pernice transcription per the scope doc, with its two acceptance
   tests, scored on exam v2. *(the work)*
4. **G4** — if S4 is positive, MRI potency becomes the within-dial rank input.
5. G8, G9, G10 in that order, each only after the gate above it is green.
