"""B2 — Cell brick (scGPT): similarity-weighted cross-cell-type transfer.

THE QUESTION
------------
Given how *other* immune cell types respond to IFN-beta, predict a held-out cell
type's response — well enough to beat the mean-shift nulls. Beating them is the
only way to earn the words "cell-type specific": IFN-beta's interferon signature
is largely shared, so reproducing the average response is easy and means nothing.

Two bars, both reported (see bricks/baselines.py for why they differ):
  leave-one-out mean shift  0.8166   THE CANONICAL BAR — the fair one
  global mean shift         0.8498   leaky in the null's favour; harder
This model scores 0.8732 and clears both.

THE MODEL
---------
The mean-shift null averages every training cell type's delta with equal weight.
That is exactly its weakness: it insists a megakaryocyte responds like a
monocyte. This model keeps the averaging but *weights* it — training cell types
whose resting (control) transcriptome resembles the held-out cell type's get
more say:

    w_j  = softmax_j( sim(X, j) / T )          over training cell types j
    dhat = sum_j w_j * unit(delta_j)           unit-normalized, so one
                                               large-amplitude cell type
                                               cannot dominate the direction
    pred = max(dhat * scale, -control_X)       predicted expression stays >= 0

`sim` is a cell-state similarity computed from **control data only** — either
scGPT <cls> cell embeddings (`similarity="scgpt"`) or plain Pearson correlation
of the control mean profiles (`similarity="control"`). Both are legitimate under
LOCTO: neither ever touches the held-out cell type's perturbed data.

The final `max(..., -control_X)` is a floor, not a fudge: the data is
log1p-normalized and therefore non-negative, so a gene sitting at zero in the
held-out cell type cannot go down. The mean-shift null happily predicts negative
expression there and pays for it.

HONESTY RULES OBSERVED HERE
---------------------------
- **LOCTO.** Predicting cell type X uses training deltas from the other cell
  types only. Fit per held-out cell type, never once on everything.
- **Nested selection.** The temperature T is chosen by an *inner* leave-one-out
  over the 7 training cell types, separately inside each outer fold. Picking T
  by looking at the outer LOCTO score would be fitting the benchmark; the number
  this brick reports is not selected on the fold it is reported for.
- **A negative control is part of the scorecard.** Scrambling which similarity
  goes with which training delta collapses the model to ~0.78, below even the
  plain leave-one-out mean shift. The gain comes from the similarity, not from
  a leak.
- **We clear the fair bar AND the leaky one.** The global mean shift averages
  the true delta over ALL cell types including the held-out one, so it has seen
  data this model never gets — leaky in the null's favour. We did not quietly
  swap to the friendlier comparison: both are printed, canonical first.
- **One fold cannot be scored at all.** Megakaryocytes (63/69 cells) have a
  measured-delta reliability of 0.045, so that fold scores noise for every model
  here including the nulls. The scorecard prints the aggregate restricted to the
  measurable cell types next to the headline, because the unrestricted mean
  rewards whoever got luckier on noise.

WHAT scGPT DID AND DID NOT CONTRIBUTE — full account in claims/B2-SCGPT/NOTES.

    similarity source              aggregate    7 measurable cell types
    loo mean shift (CANONICAL BAR)    0.8166              0.8689
    global mean shift (leaky)         0.8498              0.9032
    scGPT-blood cell embeddings       0.8696              0.9247
    scGPT-human cell embeddings       0.8708              0.9266
    control-profile correlation       0.8732              0.9293
    noise ceiling (perfect model)     0.8925              0.9898

scGPT is real here and it works: the encoder validates at 94% 1-NN cell-type
accuracy on Kang control cells, and its cell embeddings beat the bar as the
transfer metric. What it does NOT do is beat plain Pearson correlation of the
control profiles — a metric that costs nothing (paired over 8 folds, mean
difference -0.0036, Wilcoxon p=0.55: no detectable difference). scGPT's
pretrained *gene* embeddings contributed nothing at all: two smoother
formulations (rank-512 kernel and sparse kNN) both degraded the score, and the
nested selection zeroed them out in all 8 folds.

Scope of that negative result, so it does not widen in the retelling: scGPT was
tested as a cell-state SIMILARITY METRIC and as a GENE-EMBEDDING SMOOTHER. It
was NOT tested in the GEARS-style fine-tuned perturbation mode it is actually
promoted for, because that needs the MLM decoder head, which did not validate
(see scripts/scgpt_model.py). Supported: "scGPT buys nothing over np.corrcoef as
a cell-state metric on this benchmark." NOT supported: "scGPT doesn't work for
perturbation prediction." 

So the honest label on this brick is **cell-state transfer beats the bar**, not
"scGPT beats the bar". The default similarity is therefore the cheap one; the
scGPT path is kept, reproducible, and reported, because "we tried the foundation
model and it did not add over a correlation" is the finding, and it should stay
falsifiable rather than be deleted.

Interfaces:
  - PerturbationModel:  predict(control_mean, cell_type) -> perturbed_mean
  - Stage:              __call__(state) -> state  writes state["cell_delta"]
"""

from __future__ import annotations

import os

import numpy as np

from backtest.harness import PerturbationBenchmark

_CACHE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "data", "cache")

# Temperatures offered to the inner leave-one-out. Wide enough to include the
# near-uniform end (T=1 is close to the mean-shift null) and the near-nearest-
# neighbour end, so the selection can honestly land on "don't weight at all".
_TEMPS = (1.0, 0.5, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.015, 0.01, 0.007)


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float).ravel()
    b = np.asarray(b, float).ravel()
    if a.std() < 1e-12 or b.std() < 1e-12:
        return 0.0
    r = float(np.corrcoef(a, b)[0, 1])
    return r if np.isfinite(r) else 0.0


def _unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v / n if n > 1e-12 else v


def scgpt_embeddings_path(checkpoint: str = "blood") -> str:
    return os.path.join(_CACHE, f"scgpt_embeddings_{checkpoint}.npz")


def load_scgpt_cell_similarity(cell_types: list[str], checkpoint: str = "blood"):
    """Cosine similarity between scGPT control-cell embeddings, in `cell_types` order.

    Produced by `scripts/scgpt_embed.py` running in `.venv-scgpt` — the heavy
    torch environment never has to be importable here. Returns None (so the
    caller can fall back) if the cache is absent or does not cover these cell
    types, rather than silently substituting a different metric.
    """
    path = scgpt_embeddings_path(checkpoint)
    if not os.path.exists(path):
        return None
    # No allow_pickle: the cache holds only plain numeric and unicode-string
    # arrays, so it loads without enabling arbitrary-object deserialization.
    d = np.load(path)
    names = [str(x) for x in d["cell_types"]]
    if not set(cell_types) <= set(names):
        return None
    emb = d["cell_emb_ctrl"].astype(float)
    pos = {n: i for i, n in enumerate(names)}
    m = np.vstack([emb[pos[ct]] for ct in cell_types])
    m /= np.maximum(np.linalg.norm(m, axis=1, keepdims=True), 1e-12)
    return m @ m.T


class ScGPTCellModel:
    """Similarity-weighted LOCTO perturbation transfer over a PerturbationBenchmark.

    Args:
      benchmark:  the Kang backtest.
      similarity: "control" — Pearson correlation of control mean profiles
                              (DEFAULT: it scores as well as scGPT here and
                              costs nothing — no checkpoint, no torch, no
                              second venv. See the module docstring);
                  "scgpt"   — scGPT <cls> embeddings of control cells (needs the
                              npz from scripts/scgpt_embed.py);
                  "auto"    — scgpt if cached, else control.
      checkpoint: which scGPT checkpoint's embeddings to use ("blood"/"human").
      temperature: fixed T, or None to select it per fold by inner leave-one-out
                   over the training cell types (the honest default).
      floor:      clamp predicted expression at >= 0.
    """

    def __init__(self, benchmark: PerturbationBenchmark, similarity: str = "control",
                 checkpoint: str = "blood", temperature: float | None = None,
                 floor: bool = True) -> None:
        self.benchmark = benchmark
        self.checkpoint = checkpoint
        self.floor = floor
        cts = benchmark.cell_types
        self.cell_types = cts
        self._idx = {ct: i for i, ct in enumerate(cts)}

        self._control = np.vstack([benchmark.profiles[ct].control_mean for ct in cts])
        self._delta = np.vstack([benchmark.profiles[ct].true_delta for ct in cts])

        sim = None
        if similarity in ("scgpt", "auto"):
            sim = load_scgpt_cell_similarity(cts, checkpoint)
            if sim is None and similarity == "scgpt":
                raise FileNotFoundError(
                    f"no scGPT embedding cache at {scgpt_embeddings_path(checkpoint)}; "
                    f"run:  . .venv-scgpt/bin/activate && "
                    f"python scripts/scgpt_embed.py --ckpt data/cache/scgpt/{checkpoint}")
        if sim is None:
            sim = np.corrcoef(self._control)
            self.similarity = "control"
        else:
            self.similarity = "scgpt"
        self._sim = sim

        # Per held-out cell type: pick T on the training cell types only, then
        # build that fold's prediction. Nothing here reads the held-out delta.
        self._pred_delta: dict[str, np.ndarray] = {}
        self._chosen_T: dict[str, float] = {}
        for ct in cts:
            held = self._idx[ct]
            train = [j for j in range(len(cts)) if j != held]
            t = temperature if temperature is not None else self._select_T(train)
            self._chosen_T[ct] = t
            self._pred_delta[ct] = self._combine(train, held, t)

    name_prefix = "cell:scgpt"

    @property
    def name(self) -> str:
        return f"cell:{self.similarity}-sim-transfer(LOCTO,{self.checkpoint})"

    # --- core ---------------------------------------------------------------
    def _combine(self, train: list[int], held: int, temperature: float) -> np.ndarray:
        s = self._sim[held, train]
        w = np.exp((s - s.max()) / temperature)
        w /= w.sum()
        d = np.tensordot(w, np.vstack([_unit(self._delta[j]) for j in train]), axes=1)
        # Restore a plausible amplitude before the floor: the unit-normalized
        # combination has no scale of its own, and the floor is scale-sensitive.
        # Pearson itself is scale-invariant, so this only matters for the clamp.
        ref = np.linalg.norm(self._delta[train].mean(axis=0))
        d = d * (ref / max(np.linalg.norm(d), 1e-12))
        if self.floor:
            d = np.maximum(d, -self._control[held])
        return d

    def _select_T(self, train: list[int]) -> float:
        """Inner leave-one-out *within* the training cell types — never the outer fold."""
        best_t, best_score = _TEMPS[0], -np.inf
        for t in _TEMPS:
            scores = []
            for h in train:
                sub = [j for j in train if j != h]
                scores.append(_pearson(self._delta[h], self._combine(sub, h, t)))
            m = float(np.mean(scores))
            if m > best_score:
                best_score, best_t = m, t
        return best_t

    # --- PerturbationModel interface ---------------------------------------
    def predict(self, control_mean: np.ndarray, cell_type: str) -> np.ndarray:
        """Return predicted IFN-beta-perturbed mean expression, shape [n_genes]."""
        control_mean = np.asarray(control_mean, float).ravel()
        return control_mean + self._pred_delta[cell_type]

    # --- Stage interface (spine) -------------------------------------------
    def __call__(self, state: dict) -> dict:
        state["cell_delta"] = {ct: self._pred_delta[ct].copy() for ct in self.cell_types}
        return state


# ---------------------------------------------------------------------------
def _scrambled_similarity_control(bench: PerturbationBenchmark, n_rep: int = 20) -> float:
    """Negative control: keep everything, break only the similarity->delta pairing.

    If the model's gain survives this, the gain was never coming from the
    similarity and something is leaking. It should fall to roughly the plain
    leave-one-out mean shift or below.
    """
    rng = np.random.default_rng(0)
    scores = []
    for _ in range(n_rep):
        m = ScGPTCellModel(bench, similarity="control")
        k = len(m.cell_types)
        for ct in m.cell_types:
            held = m._idx[ct]
            train = [j for j in range(k) if j != held]
            shuffled = list(rng.permutation(train))
            s = m._sim.copy()
            s[held, train] = m._sim[held, shuffled]
            saved, m._sim = m._sim, s
            m._pred_delta[ct] = m._combine(train, held, m._chosen_T[ct])
            m._sim = saved
        scores.append(bench.evaluate(m)["delta_pearson"].mean())
    return float(np.mean(scores))


def _scorecard() -> None:
    from data.kang import to_benchmark
    from bricks.baselines import (GlobalMeanShiftNull, IdentityNull,
                                  LeaveOneOutMeanShiftNull)

    print("B2-SCGPT — LOCTO scorecard on real Kang IFN-beta data.\n")
    bench = to_benchmark()
    r_bar = bench.evaluate(GlobalMeanShiftNull(bench.global_mean_delta()))["delta_pearson"]
    r_ident = bench.evaluate(IdentityNull())["delta_pearson"]

    models = [ScGPTCellModel(bench, similarity="control")]
    for ckpt in ("blood", "human"):
        try:
            models.append(ScGPTCellModel(bench, similarity="scgpt", checkpoint=ckpt))
        except FileNotFoundError:
            print(f"[scGPT-{ckpt} variant not cached — run scripts/scgpt_embed.py "
                  f"in .venv-scgpt to reproduce it]")
    scores = {m.name: bench.evaluate(m)["delta_pearson"] for m in models}
    labels = [m.similarity if m.similarity == "control" else f"scgpt-{m.checkpoint}"
              for m in models]

    print(f"\n{'cell type':<20}{'bar':>8}" + "".join(f"{lab:>16}" for lab in labels))
    for ct in bench.cell_types:
        row = f"  {ct:<18}{r_bar[ct]:8.3f}"
        for m in models:
            v = scores[m.name][ct]
            row += f"{v:13.3f}{' +' if v > r_bar[ct] else ' -'}"
        print(row)

    r_canon = bench.evaluate(LeaveOneOutMeanShiftNull(bench))["delta_pearson"]
    r_canon_floor = bench.evaluate(
        LeaveOneOutMeanShiftNull(bench, floor=True))["delta_pearson"]

    print("\n-- nulls, worst to best --")
    print(f"  identity (no change)             {r_ident.mean():.4f}")
    print(f"  scrambled-similarity control     {_scrambled_similarity_control(bench):.4f}"
          f"   (our model with the metric broken)")
    print(f"  loo mean shift  [CANONICAL BAR]  {r_canon.mean():.4f}"
          f"   ({r_canon_floor.mean():.4f} with the >=0 floor)")
    print(f"  global mean shift (leaky)        {r_bar.mean():.4f}"
          f"   <-- harder secondary; it averages in the held-out delta")

    print("\n-- models --")
    best = None
    for m, lab in zip(models, labels):
        agg = scores[m.name].mean()
        wins = int((scores[m.name] > r_bar).sum())
        print(f"  {lab:<32} {agg:.4f}   beats bar in {wins}/{len(r_bar)} cell types")
        if best is None or agg > best[1]:
            best = (lab, agg)

    if bench.reliability:
        print("\n-- how much of this is even reachable --")
        ceilings = [bench.ceiling(ct) for ct in bench.cell_types]
        print(f"  ceiling for a PERFECT model      {sum(ceilings)/len(ceilings):.4f}")
        noisy = bench.unreliable_cell_types
        if noisy:
            keep = [ct for ct in bench.cell_types if ct not in noisy]
            for ct in noisy:
                print(f"  {ct} has delta reliability "
                      f"{bench.reliability[ct]:.3f} — that fold scores noise, for")
                print("  every model here including the nulls, and drags every aggregate.")
            print(f"  Restricted to the {len(keep)} measurable cell types:")
            print(f"    canonical bar {r_canon.loc[keep].mean():.4f}"
                  f"   leaky bar {r_bar.loc[keep].mean():.4f}"
                  f"   ceiling {sum(bench.ceiling(c) for c in keep)/len(keep):.4f}")
            for m, lab in zip(models, labels):
                print(f"    {lab:<30} {scores[m.name].loc[keep].mean():.4f}")

    primary = scores[models[0].name].mean()
    print(f"\nselected temperature per fold (inner LOO): "
          f"{sorted(set(models[0]._chosen_T.values()))}")
    print(f"\nmodel={primary:.4f}  bar={r_bar.mean():.4f}  BEATS={primary > r_bar.mean()}")
    print(f"(best variant: {best[0]} at {best[1]:.4f})")


if __name__ == "__main__":
    _scorecard()
