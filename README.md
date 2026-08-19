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

- **Backtest harness** (`backtest/`) — validated on synthetic data (oracle 1.00 › mean-shift null › identity 0.00) and run on real biology: **Kang 2018 IFN-β PBMCs**, 8 immune cell types. Identity null 0.00; the global-mean-shift bar is **0.85**; Megakaryocytes (0.48) are the outlier the harness surfaces on its own.
- **Cell model** (`bricks/cell_transfer.py`) — cross-cell-type transfer, leave-one-cell-type-out. **Scores 0.82 — it does NOT beat the 0.85 null yet.** Said out loud, because a model that can't beat "everyone responds the same" hasn't earned "cell-type specific." (scGPT is the swing to beat it.)
- **GRN** recovers the interferon module (IFIT1/IFIT3–ISG15) + MHC-II from raw data. **QSP / ABM / barrier / readout** are toy models, directionally coherent, every one flagged `validated=False`.
- **The wedge** (`bricks/vpop.py`) — a first open Python plausible-patient generator (LHS + a rejection filter that actually rejects) pointed at a neuroimmune model. No such implementation exists on GitHub.

**Reproduce it:**
```bash
python -m pip install -r requirements.txt
python -m backtest.selftest      # the ruler validates itself
python -m backtest.run_kang      # real IFN-β backtest — the 0.85 bar
python spine/run_demo.py         # the full pipeline, end to end, self-labeling
python -m results.experiment     # virtual cohort across real trial arms
python -m backtest.clinical      # two-directional clinical gate vs REAL trial outcomes (defines "viable")
python -m pytest                 # invariant tests (needs requirements-dev.txt)
```

See **[`results/RESULTS.md`](results/RESULTS.md)** for every verified number and an explicit real-vs-toy table, and `results/figures/` for the deck figures. Full brick/data map + citations: `../ms-twin/docs/RESEARCH_FINDINGS.md`.

## Honest state of the science
Every disease/PK/ABM parameter is illustrative, not fitted. The pipeline separates the therapies that *worked* (IFN-β, glatiramer) from untreated — and **openly cannot yet flag the one that harmed patients** (APL CGP77116, halted in Phase II). That gap is exactly what the backtest exists to close, and it is not hidden. **Nothing in this repo is evidence about multiple sclerosis.**

## Layout
- `spine/` — the multi-scale orchestration layer (ours)
- `bricks/` — the ten bricks (cell, GRN, QSP, ABM, barrier, intervention, readout, vpop, …)
- `data/` — loaders + backtest anchors (open datasets only; no raw patient data in-repo)
- `backtest/` — the validation harness. Trust nothing until it replays known history.
- `results/` — end-to-end experiment, RESULTS.md, deck figures
- `tests/` — invariant tests (the ruler discriminates, the wedge rejects, the pipeline never self-reports as validated)

## License
This repository's own code is **Apache-2.0** (`LICENSE`, `NOTICE`). It bundles **no** model weights and **no** datasets — loaders fetch open data at runtime into a git-ignored cache, and every runtime dependency is permissive (BSD/Apache). Restricted models named in the roadmap are optional, user-supplied, and never distributed here; their licenses travel with them. Full inventory and the boundary: **[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md)**.
