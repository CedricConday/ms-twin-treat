# ms-twin-treat

**A validated in-silico testbed for MS *treatment* candidates.**
Simulate a therapy for *established* multiple sclerosis, and backtest it against known clinical outcomes before anyone trusts a prediction.

## What this is
A multi-scale simulation loop — cell → cell-population → tissue/barrier → clinical readout — wired to a backtest harness. You feed it a therapy with a *known* result; if the sim reproduces that result, the loop has earned the right to predict an unknown one. Nothing is trusted until it replays history.

## What this is NOT
Not a cure. Not a trial replacement. Not a promise to patients. **In-silico evidence de-risks and augments; it does not replace a pivotal trial or a control arm — anywhere, for any therapy.** This tool separates good candidates from doomed ones before $200M is spent on the wrong one. We demonstrate; we do not assert. Overclaiming to the MS community is the one thing that ends this project.

## Scope
This repo models intervention on someone who *already has* MS. It is the whole project: the research notes that used to live in a separate `ms-twin` repo are in [`docs/research/`](docs/research/) since 2026-09-25, and a *prevent* arm that was once planned as a second repo was never built.

## What runs today (v0)
The whole spine is built and runs end to end — **and it is honest about being a skeleton, not a validated model.** The backtest harness came first, on purpose: nothing is trusted until it replays known history.

- **Backtest harness** (`backtest/`) — validated on synthetic data (oracle 1.00 › mean-shift null › identity 0.00) and run on real biology: **Kang 2018 IFN-β PBMCs**, 8 immune cell types. Two nulls, reported side by side: the **canonical bar** is leave-one-out mean shift (**0.8166**) — the fair one, restricted to what a held-out model may see; the global mean shift (**0.8498**) averages in the held-out cell type's own delta, so it is leaky in the null's favour and kept as the *harder* secondary. Every score also carries its **ceiling** — how much of the measured delta is real rather than sampling noise. That is how Megakaryocytes (63/69 cells, reliability **0.045**) is exposed as a fold nothing can score, dragging every aggregate by ~0.05.
- **Cell model** (`bricks/cell_scgpt.py`) — control-state-similarity transfer, leave-one-cell-type-out. **Scores 0.8732 — it beats both bars**, wins 7/8 cell types, and is negative-control checked (scramble the similarity and it collapses to 0.78, below baseline, so the gain is not a leak). On the 7 *measurable* cell types: bar 0.9032, model 0.9293, ceiling 0.9898 — it closes ~30% of the reachable gap.
- **scGPT was the intended swing and did not earn the win.** Its embeddings clear the bar (0.8696) but are statistically indistinguishable from plain `np.corrcoef` of the control profiles (paired over 8 folds, p=0.55). Kept, reproducible, and reported anyway, because a negative result should stay falsifiable rather than be deleted. `bricks/cell_transfer.py` (affine, 0.8204) is the earlier attempt that did *not* clear the bar, kept for comparison.
- **GRN** recovers the interferon module (IFIT1/IFIT3–ISG15) + MHC-II from raw data. **QSP / ABM / barrier / readout** are toy models, directionally coherent, every one flagged `validated=False`.
- **The wedge** (`bricks/vpop.py`) — a first open Python plausible-patient generator (LHS + a rejection filter that actually rejects) pointed at a neuroimmune model. No such implementation exists on GitHub.

**Reproduce it** (in a venv — the pins matter, a system Python gives different numbers):
```bash
python -m pip install -r requirements-dev.txt   # harness, toy pipeline, tests, lint
python -m ruff check .           # lint gate
python -m pytest                 # invariant tests + scorer input guards
python -m backtest.selftest      # the ruler validates itself
python spine/run_demo.py         # the full pipeline, end to end, self-labeling

python -m pip install -r requirements-data.txt  # + scanpy, for the real Kang path only
python -m backtest.run_kang      # real IFN-β backtest — both bars + the reachable ceiling
python -m bricks.cell_scgpt      # the cell model's full LOCTO scorecard vs every null
python -m results.experiment     # virtual cohort across real trial arms
python -m backtest.clinical      # two-directional clinical gate vs REAL trial outcomes (defines "viable")
```
Full setup, pins and gates: **[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md)**. What the
repo's evidence labels claim, and what they don't: **[`docs/QUALITY.md`](docs/QUALITY.md)**.

See **[`results/RESULTS.md`](results/RESULTS.md)** for every verified number and an explicit real-vs-toy table, and `results/figures/` for the deck figures. Full brick/data map + citations: [`docs/research/RESEARCH_FINDINGS.md`](docs/research/RESEARCH_FINDINGS.md).

## Honest state of the science
Every disease/PK/ABM parameter is illustrative, not fitted. The pipeline separates the therapies that *worked* from untreated, and directs the ones that **harmed** patients correctly (APL CGP77116 and IFN-γ) — that harm emerges from an `immunogenic` parameter set from each drug's *documented mechanism*, not fitted to its relapse number.

> **Live numbers: [`results/STATE.md`](results/STATE.md)**, regenerated by `python scripts/state_of_build.py`, which runs every gate and writes its output with a timestamp and the git SHA. Three separate figures in this repo went stale in a week because they were written down by hand; if a number in any document disagrees with that file, the document is wrong.

**The clinical gate is scored on 17 arms and currently gets 10 of 17 directions and 1 of 12 magnitudes right** (re-measured 2026-09-20; `python -m backtest.clinical`, anchors in [`docs/TRIAL_ANCHORS.md`](docs/TRIAL_ANCHORS.md)). It read 9/14 and 1/9 on 2026-09-17 and 4/4 on four arms before that. Nothing about the model improved or regressed when that changed — the exam got harder, and the 4/4 version could not show either of the two things now visible:

- **lenercept and atacicept are immunosuppressive by mechanism and harmed patients in trials.** A rule that maps mechanism class to clinical direction predicts benefit for both and is wrong on both. Lenercept now carries a regulation-disruption flag earned by its TNF biology (TNF-deficient mice get *worse* EAE; TNFR2 is required for remyelination), which moves it from −77% to −35% without crossing zero — the model can represent the case, and the current constants still get the sign wrong.
- **one strength per class cannot separate two drugs inside a class.** Ocrelizumab vs interferon beta-1a (OPERA), alemtuzumab vs interferon beta-1a (CARE-MS I) and ponesimod vs teriflunomide (OPTIMUM) are all scored against an active comparator, and the simulation returns exactly 0% difference for each.

**The out-of-sample tests run too, and every one of them loses to its null** (re-measured 2026-09-21 on the 15-arm exam; live values in [`results/STATE.md`](results/STATE.md)):

| test | what is held out | result |
|---|---|---|
| `backtest.loo` | one arm | **30.0pp MAE vs an 11.3pp predict-the-mean null** |
| `backtest.lomo` | a whole mechanism | **45.4pp vs 11.8pp** |
| `backtest.lomo_capacity` | a whole mechanism, under the one structural fix that flips the depletion sign | **42.2pp vs 11.8pp** with K fitted per fold; **48.0pp at the fixed K=2000 the screen runs at** — *worse* than the transcription |

A single class strength carries no drug-specific information, and on an unseen
mechanism the model answers with roughly the training arms' average whatever it
is asked. **That is the honest headline for this build**, and it is why
`screen/` kills candidates but refuses to rank them.

Both failures are the point of the arm set, not defects in it. Real viability still needs data-grounded parameters, **out-of-sample** arms, and validated magnitudes. **Nothing in this repo is evidence about multiple sclerosis.**

## Layout
- `spine/` — the multi-scale orchestration layer (ours)
- `bricks/` — the ten bricks (cell, GRN, QSP, ABM, barrier, intervention, readout, vpop, …)
- `data/` — loaders + backtest anchors (open datasets only; no raw patient data in-repo)
- `backtest/` — the validation harness. Trust nothing until it replays known history.
- `results/` — end-to-end experiment, RESULTS.md, deck figures
- `tests/` — invariant tests (the ruler discriminates, the wedge rejects, the pipeline never self-reports as validated)

## License
This repository's own code is **Apache-2.0** (`LICENSE`, `NOTICE`). It bundles **no** model weights and **no** datasets — loaders fetch open data at runtime into a git-ignored cache, and every runtime dependency is permissive (BSD/Apache). Restricted models named in the roadmap are optional, user-supplied, and never distributed here; their licenses travel with them. Full inventory and the boundary: **[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md)**.
