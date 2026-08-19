# SCGPT_HANDOFF — 2nd instance, one job: beat the bar

**You are a new Claude Code instance. Owner name: `sg` (scGPT).**
Work in `~/repos/ms-twin-treat/`. Read this fully, run the cold-boot check, claim
your lane, then build. Everything you need is here.

---

## 0. THE ONE JOB
The cell brick (B2) currently scores **0.82 aggregate delta_pearson** on the Kang
IFN-β backtest — it does **NOT** beat the global-mean-shift null (**0.85**, "the
bar"). Your job: produce a cell-state model — **scGPT** the primary swing — that
**beats 0.85** on the same honest leave-one-cell-type-out (LOCTO) evaluation.

If scGPT beats it: we have "a model that beats the null," the single number that
upgrades the whole story. If it doesn't: **return an honest account of why**
(zero-shot cross-cell-type transfer is genuinely hard; scGPT's perturbation
variant is GEARS-style and may not fit this setup). **A truthful "did not beat it,
here is the number and the reason" is a WIN. A faked win is the one unforgivable
outcome — it poisons an MS patient project. Do not fake it.**

## 1. COLD-BOOT CHECK (~1 min)
```bash
cd ~/repos/ms-twin-treat
. .venv/bin/activate
python -m backtest.selftest      # HARNESS VALIDATED
python -m bricks.cell_transfer   # shows the 0.82 you must beat, per cell type
git log --oneline | head -3
```

## 2. CLAIM YOUR LANE (atomic, do first)
```bash
cd ~/repos/ms-twin-treat
mkdir claims/B2-SCGPT 2>/dev/null && {
  echo "owner: sg" > claims/B2-SCGPT/OWNER
  echo "started: $(date -u +%FT%TZ)" >> claims/B2-SCGPT/OWNER
  echo "claimed B2-SCGPT"
} || echo "ALREADY CLAIMED: $(cat claims/B2-SCGPT/OWNER)"
```

## 3. FILE OWNERSHIP — ABSOLUTE (this is how three instances don't collide)
- **You write ONLY `bricks/cell_scgpt.py`** (+ your own scratch/scoring script).
- **DO NOT edit** `spine/`, `run_demo.py`, `bricks/cell_transfer.py`, or any other
  brick — those are `vps`-owned. `run_demo.py` is `vps`'s; `vps` wires your model
  in AFTER you've scored it. Leave wiring instructions in
  `claims/B2-SCGPT/NOTES` (import path, class name, constructor) and `vps` does
  the integration. You never touch the wire.

## 4. THE INTERFACE YOU MUST IMPLEMENT (freeze — do not renegotiate)
Your model is a harness `PerturbationModel`:
```python
class ScGPTCellModel:
    name = "cell:scgpt-..."
    def predict(self, control_mean: np.ndarray, cell_type: str) -> np.ndarray:
        """Return predicted PERTURBED (IFN-β) mean expression, shape [n_genes]."""
```
That's it. `predict` takes a cell type's control mean-expression vector and
returns its predicted IFN-β-stimulated mean. The harness computes the delta and
scores it. **LOCTO is mandatory:** when predicting cell type X, your model must
NOT have trained on X's perturbed data. If scGPT is zero-shot, that's automatic;
if you fit anything, hold X out.

## 5. HOW TO SCORE YOURSELF (no dependency on vps)
```python
from data.kang import to_benchmark
from bricks.baselines import GlobalMeanShiftNull
bench = to_benchmark()                       # 8 cell types, real Kang IFN-β
model = ScGPTCellModel(...)                  # yours
r_model = bench.evaluate(model)["delta_pearson"].mean()
r_bar   = bench.evaluate(GlobalMeanShiftNull(bench.global_mean_delta()))["delta_pearson"].mean()
print(f"model={r_model:.4f}  bar={r_bar:.4f}  BEATS={r_model > r_bar}")
```
`bench.profiles[ct].control_mean` and `.perturbed_mean` give you the raw vectors;
`bench.gene_names` the gene order. The data is log1p-normalized already.

## 6. ENVIRONMENT — READ THIS, IT SAVES YOU AN HOUR
- **CPU only, no GPU.** Use scGPT's SMALL checkpoint; expect slow inference. Plan
  batch sizes accordingly. The blood + brain checkpoints are the MS-relevant ones.
- **DO NOT install torch/scGPT into the main `.venv`.** Its deps (torch, scanpy
  pins, possibly flash-attn) will fight the pinned analysis stack (numpy 2.2.6)
  and could break the harness for everyone. **Make a SEPARATE venv:**
  `python3 -m venv .venv-scgpt` (gitignored already via `.venv` pattern — add
  `.venv-scgpt/` to `.gitignore` if needed).
- **Clean hand-off pattern to dodge ABI hell:** in `.venv-scgpt`, run scGPT,
  compute predicted perturbed means per cell type, and **save them to
  `data/cache/scgpt_pred.npz`** (gitignored). Then `bricks/cell_scgpt.py`, running
  in the MAIN `.venv`, just loads that npz and serves it through `predict()`.
  Decouples the heavy model env from the pinned scoring env entirely.
- scGPT weights: `bowang-lab/scGPT` (GitHub) / HF checkpoints. If download or
  install fights the box for more than ~45 min, say so in NOTES and fall back to
  scGPT **gene embeddings as features** for a transfer model — still a legitimate
  "scGPT-based" attempt.

## 7. GUARDRAILS (identical for every instance)
- Author **always** `cedric@condaydigital.com`. **NEVER** the gmail. Set it:
  `git config user.email "cedric@condaydigital.com"; git config user.name "Cedric Conday"`.
- **Push held until Cedric says go.** Local commits only. No remote.
- No weights/datasets committed (they're gitignored; keep it that way).
- **Built ≠ validated.** If your model beats the bar, that is a real backtest
  result on ONE dataset (IFN-β) — it is NOT "we can predict MS drugs." Label it.
- Commit when you have a scored result (win or honest loss).

## 8. DEFINITION OF DONE
`bricks/cell_scgpt.py` implements `PerturbationModel`; a scoring run prints your
LOCTO `delta_pearson` vs the 0.85 bar with `BEATS=True/False`; `claims/B2-SCGPT/
NOTES` tells `vps` how to wire it; and you state the honest verdict — the number,
and if it lost, why. Then `echo "done: $(date -u +%FT%TZ)" >> claims/B2-SCGPT/OWNER`
and commit.

## 9. CONTEXT (optional, for grounding)
Full research + brick map: `~/repos/ms-twin/docs/RESEARCH_FINDINGS.md` (scGPT is
B1 there: MIT, brain+blood checkpoints, the commercial-clean cell brick). Why the
bar is 0.85 and hard: the interferon signature is largely shared across immune
cell types, so beating "everyone responds the same" requires real cell-type-
specific signal. That's your target.
