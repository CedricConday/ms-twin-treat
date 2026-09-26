# The four-day plan — 2026-09-26 to 2026-09-29

**For a fresh instance: read this file, then `docs/GATE_GAP_ANALYSIS.md`, then the
Result section at the end of `docs/EXAM_V2_PREREG.md`, then the last two entries
of BUILD_PLAN.md §8.4. Do not read §8.1's blockers as a worklist. Everything
below is dated; strike a line when it is done, never rewrite it.**

## Definition of done

The gate is green when `PYTHONPATH=. python -m backtest.exam_v2` prints an S1
interval-LOMO error at or below **16.0pp** against its 20.0pp null (the 20%
margin in `gate/criterion.py`), the fold-gap CI excludes zero, and the
placebo-only subset also clears. Then `results/STATE.md` is regenerated and
`docs/DECISION_GATE.md` is rewritten to say so, with the date.

If day 4 ends red, done is a dated entry in BUILD_PLAN §8.4 with the measured
number of each model tried and the one structural reason each failed. A red
with its reason is a result; a red without one is not.

## What is settled, so nobody re-derives it

- The point-ARR exam cannot be won by anyone: 4.5pp headroom. Exam v2 can:
  10.9pp headroom, 54% of the null. **Score every model on exam v2.**
- Cohort size is not the problem: 1024 seeds changed nothing (017cb50).
- The current model fits potency 0.00 on exam v2 because three dial patterns
  point the wrong way. Any replacement must keep the sign of depletion, keep
  the sign of immune challenge, and predict nothing for costimulation block.
- The MRI channel orders drugs within a dial (tau 0.80, p 0.024). It is the
  within-dial ranking input once a model has the right form (gap G4).
- No open trial-replaying MS simulator exists; the closed one calibrates to the
  outcome it reports (`~/.claude` memory `ms-twin-treat-cold-start`, and the
  2026-09-26 landscape check; a table of it was stopped by the bio
  classifier, so keep landscape prose short).

## Day 1 — 2026-09-26 — probe A: the two-equation model

Source vendored at `docs/research/jenner2026/sections_2_to_6.md` (open
access, equation 2.1 transcribed there). Two state variables, four rates.

1. `bricks/qsp_minimal.py`: transcribe equation (2.1) exactly, fixed-step
   integrator, month time unit, parameter ranges as the paper states them.
   Reproduction target: the paper's own result that raising the disease
   strength parameter carries the system through a Hopf bifurcation from a
   stable state into limit cycles (its section 4). A test pins that.
2. Readout: relapse proxy = number of inflammation peaks per year above the
   untreated baseline (same detector convention as `qsp_velez.relapse_events`,
   scored on the untreated arm's baseline). Lesion ratio = time-integrated
   inflammation, treated over untreated. Sormani map unchanged.
3. `bricks/profiles_minimal.py`: each arm as multipliers on the four rates,
   from pharmacology, never from outcome. This model has no compartment and
   no cell types, so many arms will share a pattern; that is the probe's
   limit and it is stated in the module docstring.
4. `backtest/exam_v2.py` gets a `--model minimal` switch that swaps the
   response-table builder and profile map. Cache to
   `results/exam_v2_curve_minimal.json`. Deterministic model: one seed.
5. Run. Record S1, S2, S3 in BUILD_PLAN §8.4 with the date. Commit and push.

Decision at end of day 1: if S1 beats the null, the rest of the week hardens
and ranks (skip to day 4). If not, the sign table says which of the three
acceptance tests it failed, and day 2 starts the port.

~~Day 1 items 1-5~~ **done 2026-09-26**: S1 43.8pp vs 30.6pp null, loses;
tests 1 and 2 of the three fail structurally (BUILD_PLAN §8.4, same date).
The port starts. **Also 2026-09-26: the day boundaries are dropped on
Cedric's instruction; the remaining items run in order, today, none skipped.**

## Day 2 — 2026-09-27 — Pernice 2020 port, structure and reproduction

~~Items 1-3~~ **done 2026-09-26**: `bricks/qsp_pernice.py`, 10 of 12 Figure S2
landmarks, readings annotated; two-year figures not reproduced (BUILD_PLAN §8.4).
Tau-leaping not built, per item 1's condition.

Sources vendored at `docs/research/pernice2020/` (CC BY paper text and
supplement; the net diagram as `fig3_net.jpg`; the deterministic 30-day
HD-vs-MS solution as `figS2_reproduction_target.png`). Scope, licence and
the parameter tables are in `docs/PERNICE_PORT_SCOPE.md`.

1. `bricks/qsp_pernice.py`: 26 places (21 named plus ODC unfolded to 5
   myelination levels), 55 transitions. Mass-action ones from the net
   diagram; the 15 general ones from Additional file 1 §S1.1, verbatim
   functions. Table S1 (MS column where it differs) and Table S2 pinned in a
   constant block with the table reference beside each value. Initial
   marking from the paper's Table 2. Antigen injection schedule from the
   paper (days 2, 67, 127, 295, 300, 303, 307, 600). Deterministic ODE first;
   tau-leaping second only if the deterministic run reproduces.
2. Reproduction target: Figure S2, healthy versus MS, 30 days, the ODC
   irreversibly-damaged curve rising to about 13 in the MS configuration
   and staying near 0 in the healthy one. A test pins the two configurations
   differing in exactly two parameters.
3. Arc ambiguities that the diagram does not settle get a `# UNRESOLVED`
   comment naming the two readings and which one was taken. Do not guess
   silently.

## Day 3 — 2026-09-28 — Pernice arms and scoring

~~Items 1-3~~ **done 2026-09-26**: S1 35.1pp vs 19.8pp (primary reading), 31.2pp
(memory-read variant); loses. Deterministic, one run per cell, 30 s.

1. `bricks/profiles_pernice.py`: arms as multipliers on the port's own named
   transition rates, from pharmacology. Peripheral killing, BBB passage,
   activation and the daclizumab-shaped pair are separate dials here, which
   is the point of the port.
2. Response tables at 128 seeds if stochastic, one run if deterministic;
   cache under `results/exam_v2_curve_pernice*.json`. Budget the runtime
   before launching; 26 states at 730 days is roughly 10x the current model.
3. Run exam v2 with `--model pernice`. Record. Commit and push.

## Day 4 — 2026-09-29 — the best model, ranked and written up

1. For whichever model scored best: wire the MRI lesion ratio as the
   within-dial rank input (gap G4), re-run `screen/report.py`, regenerate
   `results/STATE.md` with `scripts/state_of_build.py`.
2. Rewrite `docs/DECISION_GATE.md` to the measured state, green or red.
3. Second exam, if time: the 2024 Cochrane network meta-analysis
   (PMC10765473) gives a network risk ratio per drug for one or more
   relapses at 24 months across 16 to 20 treatments, four of which the exam
   lacks. Its summary-of-findings tables are in the open-access XML; the
   per-trial forest data is behind the Cochrane Library's bot wall and needs
   a browser session on Cedric's machine.

## House rules that apply to every day

- Pre-register before measuring: a new scorer or exam gets a dated rules
  file committed before its first run. Exam v2 is already registered; a new
  model does not need a new registration, it is scored on the existing one.
- Commit and push per fix; run `pre-push` first. Author is the noreply.
- Transcribe, cite, tune nothing. A number in a comment comes with the way
  to re-derive it.
- No unrequested hygiene work. Scanners, provenance, doc rewrites cost a
  night once already.
- If the bio classifier stops a reply, do not retry the content. Keep prose
  about scoring and code; keep mechanism tables out of chat.
