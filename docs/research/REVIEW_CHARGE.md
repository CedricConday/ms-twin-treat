# MS-Twin — Adversarial Review Charge

**How to run this (Cedric's notes, not part of the prompt):** send this document + both attached files to **Fable** and **K3 separately**. Identical prompt, no cross-contamination — neither sees the other's output in Round 1. Both answer *every* question, so each objection gets two independent votes and the disagreements are the signal. Round 2 = each reads the other's answers to the same questions and adjudicates. Then I synthesize. Nobody here is validating anything; they're here to find the reason this dies.

---

## Context

A solo engineer (strong at multi-scale simulation, data pipelines, backtesting; **not** a bench biologist) is scoping an **open, validated in-silico testbed for MS interventions** — a pipeline that simulates a candidate MS therapy and is *backtested against known clinical outcomes* before anyone trusts it. Explicit, non-negotiable framing: **this is a tool that separates good therapy candidates from doomed ones — NOT a cure, NOT a trial replacement, NOT a promise to patients.** His wife has MS; overclaiming to that community is the one failure that ends the project.

Attached: the project scaffold (`PROJECT_MS_TWIN.md`) and the research findings from ~490k tokens of OSS/literature/regulatory survey across 4 parallel agents (`RESEARCH_FINDINGS.md`). Everything in the findings doc is claimed-verified. **Assume nothing is verified until you have reasoned about whether it *could* be true.** Your job is to find the fatal flaw, not to improve the wording.

**You are one of two independent reviewers running this identical charge.** You will not see the other's answers until Round 2. Do not hedge, do not cover the waterfront, and do not soften a call because you assume someone else will make it. Where you are uncertain, say *uncertain* and say what would resolve it — an unmarked guess from you is worse than a gap, because the whole point of running two of you is that a confident wrong answer gets caught.

**Standing rule:** the honest bottom line is *"this de-risks candidates; it cannot replace a trial."* If that line is being quietly walked back anywhere in the plan — in the architecture, the framing, the funding story, the public-repo strategy — flag it hard. That line is the project's spine.

---

# PART A — Science & feasibility

Attack the *science and the build*, not the pitch.

**A1. The bricks.** For each proposed component — scGPT as the cell brick, PhysiCell + `MS_ABM_Weatherley` for population, PK-Sim/Verscheijden for the barrier, the GRN (SCENIC/arboreto) and QSP/ODE (Tellurium/COPASI) layers — is it actually fit for *this* purpose, or is it a plausible-sounding name that doesn't do what the plan needs? Where is a brick being asked to do something it was never built for?

**A2. The spine.** The plan admits "scale-bridging is the open research problem." Is it *tractably* open, or the kind of open that eats a decade? Is a solo engineer wiring foundation-model outputs into an agent-based model committing a category error — passing embeddings between scales that don't compose?

**A3. The backtest.** Is the v0 anchor (Kang GSE96583, IFN-β on PBMCs) a *real* validation of anything about MS, or a toy that proves the pipeline runs without proving it's right? What would a *genuine* falsifying backtest look like, and does the plan have the data for it? Push hard on the honest gap — no open single-cell data for the major DMTs, so the monoclonals only validate at the clinical readout layer.

**A4. Kill question.** What is the single most likely reason this produces a demo that looks alive but predicts nothing? Name it.

---

# PART B — Strategy, credibility & regulatory

**For this part only, assume Part A is solved.** Grant the science works exactly as the plan hopes, then attack everything else. If your real answer to a B question is "it doesn't matter, the science fails," that belongs in Part A — say it there, then come back here and answer B on its own terms. A strategy pass that collapses into "but the science is shaky" is a non-answer and wastes half this review.

**B1. The wedge.** The claimed differentiator: no open Python implementation of rigorous virtual-population methods (Allen/Rieger plausible-patients, Schmidt MAPEL) has ever been run on a neuroinflammation model. Is that a real, defensible moat — or a niche nobody occupies *because it isn't valuable*? A competitor (Russo 2026, NX210c QSP/RRMS) is already probing MS virtual-pops and withholding the code. Does the clock argument hold, or is second-and-open a losing position?

**B2. The credibility ceiling.** This lives or dies on "demonstrate, don't assert" to an MS patient community. Where does the plan *structurally* tempt overclaiming? "Built in the open" means every half-working commit is public — is radical transparency an asset here, or a way to broadcast failure to vulnerable people in real time?

**B3. The trap.** Is this a fundable, shippable slice — or a disguised PhD / $200M-lab problem that consumes a solo founder with no revenue and no endpoint? What is the smallest thing that produces *external* validation (a citation, a maintainer, a grant, a collaborator) inside 90 days, and does the plan aim at it or at the cathedral?

**B4. Kill question.** If this quietly fails, what does it cost *him* — reputation, the community's trust, his wife's hope? Is that risk priced in?

---

# Required output format

Answer in this exact shape. The two reviews get diffed line against line, so structure matters more than prose.

```
## PART A
A1 — <verdict in one line>  |  confidence: high / medium / low
    <argument, max 6 lines>
A2 — ...
A3 — ...
A4 — <the named failure mode, one sentence>

### Top 3 technical objections (ranked, hardest first)
1. <objection> — WHAT WOULD CHANGE MY MIND: <specific evidence>
2. ...
3. ...

## PART B
B1 — <verdict in one line>  |  confidence: high / medium / low
    <argument, max 6 lines>
B2 — ...
B3 — ...
B4 — <the cost, named>

### Top 3 strategic objections (ranked, hardest first)
1. <objection> — CHEAPEST DE-RISKING EXPERIMENT: <concrete, ≤2 weeks>
2. ...
3. ...

## VERDICT
GO / NO-GO / GO-IF: <condition>
FATAL FLAW: <the one thing most likely to kill this>
SMALLEST FIRST MOVE: <one action, doable solo, this month>
```

---

# Round 2 — Cross-examination (sent after both return)

You now have the other reviewer's answers to the **same** questions.

1. **Every item where you disagree** (different verdict, or same verdict at different confidence): name it, say who is right, and say why. Do not split the difference — if you're changing your answer, change it and say so.
2. **Compounding.** Where does one of their objections make one of yours *worse*? A technical flaw and a strategic flaw that multiply are the real finding — e.g. "their 90-day external-validation demand is unreachable *because* my A2 means there's nothing demonstrable in 90 days."
3. **Agreement is not confirmation.** Where you both said the same thing, ask once whether you're both wrong for the same reason — same training data, same obvious read, same missing domain knowledge.

Keep it to one page. Revised verdict at the bottom, in the same VERDICT block format.

---

# What I need out of this

Not a grade. A **go / no-go with the fatal flaw named**, and if go — the *one* smallest first move that survives both teardowns. Under-promise, over-deliver. If the honest version isn't fundable, better to know now than after a public failure.
