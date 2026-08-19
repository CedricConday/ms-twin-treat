# Architecture — ms-twin-treat

## The split (why this is its own repo)
`ms-twin` (the brain) holds research + the shared method. This repo holds the **treat** build. A future `ms-twin-prevent` will hold the prevent build. The shared skeleton (spine, brick adapters, virtual-population engine, backtest harness) is factored into `ms-twin-core` **only when both arms genuinely use it** — not before. Premature abstraction is a trap; so is a forked spine. Build treat, watch where prevent will actually reuse code, factor the core the day they cross.

## The spine (what WE build)
State flows: `molecular → cell → cell-population → tissue/barrier → clinical readout`.
Scale-bridging is the open research problem — this is research, not glue. The bricks are other people's models; the spine and the validation harness are ours. Shaped like a multi-scale AnyLogic sim wired to a data pipeline — the actual skillset.

## The bricks (verified — full detail in `ms-twin/docs/RESEARCH_FINDINGS.md`)
| Layer | Pick | License | Note |
|---|---|---|---|
| Cell | scGPT | MIT | brain + blood checkpoints = MS cell types; laptop-scale |
| Population | PhysiCell + `MS_ABM_Weatherley` seed | BSD / MIT | BioFVM couples drug field to agents |
| Barrier | PK-Sim + Verscheijden brain-PBPK | GPLv2 / open R | Verscheijden ships runnable R in supplement |
| Backtest data | Kang GSE96583 (IFN-β) cellular; MSOAC clinical | open | two backtest layers, different data |
| VPop (the wedge) | Python port of Allen/Rieger + MAPEL | to build | no open Python impl exists — first-of-kind |

## The backtest discipline (the edge)
Like a trading backtest: trust no prediction until the sim replays known history.
- Feed a therapy that **failed** in trials → does the sim fail it?
- Feed one that **worked** → does it work?
- Only after it reproduces the record does it get to predict the unknown.

## Adversarial review
Fable + K3, role-differentiated (scientist vs. strategist), reviewing the design and each other. Charge lives in `ms-twin/docs/REVIEW_CHARGE.md`. This repo is what they review.
