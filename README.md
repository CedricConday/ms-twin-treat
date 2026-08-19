# ms-twin-treat

**A validated in-silico testbed for MS *treatment* candidates.**
Simulate a therapy for *established* multiple sclerosis, and backtest it against known clinical outcomes before anyone trusts a prediction.

## What this is
A multi-scale simulation loop — cell → cell-population → tissue/barrier → clinical readout — wired to a backtest harness. You feed it a therapy with a *known* result; if the sim reproduces that result, the loop has earned the right to predict an unknown one. Nothing is trusted until it replays history.

## What this is NOT
Not a cure. Not a trial replacement. Not a promise to patients. **In-silico evidence de-risks and augments; it does not replace a pivotal trial or a control arm — anywhere, for any therapy.** This tool separates good candidates from doomed ones before $200M is spent on the wrong one. We demonstrate; we do not assert. Overclaiming to the MS community is the one thing that ends this project.

## Scope
This repo is the **treat** arm only — intervention on someone who *already has* MS. A separate **prevent** arm (pre-symptomatic, genetically-susceptible virtual patient) is scoped in the research brain and will live in its own repo when it earns one. The two arms share a spine and will converge; they do not share a repo yet.

## What runs today (v0)
The whole spine is built and runs end to end — **and it is honest about being a skeleton, not a validated model.** The backtest harness came first, on purpose: nothing is trusted until it replays known history.

- **Backtest harness** (`backtest/`) — validated on synthetic data (oracle 1.00 › mean-shift null › identity 0.00) and run on real biology: **Kang 2018 IFN-β PBMCs**, 8 immune cell types. Two nulls, reported side by side: the **canonical bar** is leave-one-out mean shift (**0.8166**) — the fair one, restricted to what a held-out model may see; the global mean shift (**0.8498**) averages in the held-out cell type's own delta, so it is leaky in the null's favour and kept as the *harder* secondary. Every score also carries its **ceiling** — how much of the measured delta is real rather than sampling noise. That is how Megakaryocytes (63/69 cells, reliability **0.045**) is exposed as a fold nothing can score, dragging every aggregate by ~0.05.
- **Cell model** (`bricks/cell_scgpt.py`) — control-state-similarity transfer, leave-one-cell-type-out. **Scores 0.8732 — it beats both bars**, wins 7/8 cell types, and is negative-control checked (scramble the similarity and it collapses to 0.78, below baseline, so the gain is not a leak). On the 7 *measurable* cell types: bar 0.9032, model 0.9293, ceiling 0.9898 — it closes ~30% of the reachable gap.
- **scGPT was the intended swing and did not earn the win.** Its embeddings clear the bar (0.8696) but are statistically indistinguishable from plain `np.corrcoef` of the control profiles (paired over 8 folds, p=0.55). Kept, reproducible, and reported anyway, because a negative result should stay falsifiable rather than be deleted. `bricks/cell_transfer.py` (affine, 0.8204) is the earlier attempt that did *not* clear the bar, kept for comparison.
- **GRN** recovers the interferon module (IFIT1/IFIT3–ISG15) + MHC-II from raw data. **QSP / ABM / barrier / readout** are toy models, directionally coherent, every one flagged `validated=False`.
- **The wedge** (`bricks/vpop.py`) — a first open Python plausible-patient generator (LHS + a rejection filter that actually rejects) pointed at a neuroimmune model. No such implementation exists on GitHub.

**Reproduce it:**
```bash
python -m pip install -r requirements.txt
python -m backtest.selftest      # the ruler validates itself
python -m backtest.run_kang      # real IFN-β backtest — both bars + the reachable ceiling
python -m bricks.cell_scgpt      # the cell model's full LOCTO scorecard vs every null
python spine/run_demo.py         # the full pipeline, end to end, self-labeling
python -m results.experiment     # virtual cohort across real trial arms
python -m backtest.clinical      # two-directional clinical gate vs REAL trial outcomes (defines "viable")
python -m pytest                 # invariant tests (needs requirements-dev.txt)
```

See **[`results/RESULTS.md`](results/RESULTS.md)** for every verified number and an explicit real-vs-toy table, and `results/figures/` for the deck figures. Full brick/data map + citations: `../ms-twin/docs/RESEARCH_FINDINGS.md`.

## Honest state of the science
Every disease/PK/ABM parameter is illustrative, not fitted. The pipeline separates the therapies that *worked* (IFN-β, glatiramer) from untreated, and now also directs the one that **harmed** patients correctly (APL CGP77116, halted in Phase II) — its harm emerges from an `immunogenic` parameter set from the drug's *documented encephalitogenic mechanism*, not fitted to its relapse number.

**Read that carefully: the clinical gate at 4/4 is a capability milestone, not validation.** All four arms are ones whose outcomes informed the setup, so direction-correctness is *necessary, not sufficient*. Real viability still needs data-grounded parameters, **out-of-sample** arms, and validated magnitudes. **Nothing in this repo is evidence about multiple sclerosis.**

## Layout
- `spine/` — the multi-scale orchestration layer (ours)
- `bricks/` — the ten bricks (cell, GRN, QSP, ABM, barrier, intervention, readout, vpop, …)
- `data/` — loaders + backtest anchors (open datasets only; no raw patient data in-repo)
- `backtest/` — the validation harness. Trust nothing until it replays known history.
- `results/` — end-to-end experiment, RESULTS.md, deck figures
- `tests/` — invariant tests (the ruler discriminates, the wedge rejects, the pipeline never self-reports as validated)

## License
This repository's own code is **Apache-2.0** (`LICENSE`, `NOTICE`). It bundles **no** model weights and **no** datasets — loaders fetch open data at runtime into a git-ignored cache, and every runtime dependency is permissive (BSD/Apache). Restricted models named in the roadmap are optional, user-supplied, and never distributed here; their licenses travel with them. Full inventory and the boundary: **[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md)**.
