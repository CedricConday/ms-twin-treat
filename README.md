# ms-twin-treat

**A validated in-silico testbed for MS *treatment* candidates.**
Simulate a therapy for *established* multiple sclerosis, and backtest it against known clinical outcomes before anyone trusts a prediction.

## What this is
A multi-scale simulation loop — cell → cell-population → tissue/barrier → clinical readout — wired to a backtest harness. You feed it a therapy with a *known* result; if the sim reproduces that result, the loop has earned the right to predict an unknown one. Nothing is trusted until it replays history.

## What this is NOT
Not a cure. Not a trial replacement. Not a promise to patients. **In-silico evidence de-risks and augments; it does not replace a pivotal trial or a control arm — anywhere, for any therapy.** This tool separates good candidates from doomed ones before $200M is spent on the wrong one. We demonstrate; we do not assert. Overclaiming to the MS community is the one thing that ends this project.

## Scope
This repo is the **treat** arm only — intervention on someone who *already has* MS. A separate **prevent** arm (pre-symptomatic, genetically-susceptible virtual patient) is scoped in the research brain and will live in its own repo when it earns one. The two arms share a spine and will converge; they do not share a repo yet.

## v0 — the smallest honest loop
> one cell type · one known MS intervention · one readout · backtested against one published result

- **Cell brick:** scGPT (MIT — brain + blood checkpoints)
- **Population:** PhysiCell, seeded by the rules in `MS_ABM_Weatherley` (MIT)
- **Backtest anchor:** Kang GSE96583 — IFN-β (a real first-line MS drug) on immune cells. Does the loop reproduce the known IFN-β effect?
- **Success** = it reproduces the known outcome → validated micro-loop + a spine to extend.
- **Failure** = we learn exactly where the gap is → also the point.

## Status
Pre-v0. Scaffolding. Research and full brick/data map live in the brain repo (`ms-twin/docs/RESEARCH_FINDINGS.md`).

## Layout
- `spine/` — the multi-scale orchestration layer (ours)
- `bricks/` — adapters to the component models (other people's models, our wrappers)
- `data/` — loaders + backtest anchors (open datasets only; no raw patient data in-repo)
- `backtest/` — the validation harness. Trust nothing until it replays known history.

## License
TBD before first public push — leaning permissive for the spine/harness (ours), with per-brick licenses respected (scGPT MIT is the commercial-clean core; some models are non-commercial — see the brain's findings doc).
