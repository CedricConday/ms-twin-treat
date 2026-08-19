"""Backtest harness for perturbation prediction — ms-twin-treat.

The edge of the whole project: trust no prediction until it replays known history.

Given a *control* cell population and a *known* perturbed outcome (e.g. PBMCs
before and after IFN-beta), we score whether a model reproduces the real
per-cell-type perturbation — and, crucially, whether it beats the null that
predicts nothing changed. A model that cannot beat "predict no change" has
learned nothing about the biology, however alive the demo looks.

Design notes
------------
- v0 operates on per-cell-type *mean expression vectors*. That is the right
  altitude for a first honest loop: cheap, laptop-scale, and it makes the
  scored quantity — the perturbation delta — explicit and interpretable.
- The scored quantity is the *delta*:  perturbed_mean - control_mean, a vector
  over genes, computed per cell type. A model predicts the perturbed mean (or,
  equivalently, the delta). We compare predicted delta to true delta.
- Two nulls are provided (see baselines.py). The one that matters
  scientifically is GLOBAL MEAN SHIFT: "the perturbation does the same thing in
  every cell type." A model only earns the claim of cell-type specificity by
  beating that, not merely by beating "no change."
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
import pandas as pd
from scipy import stats


# A model is anything that, given a control mean-expression vector for one cell
# type, predicts that cell type's perturbed mean-expression vector.
class PerturbationModel(Protocol):
    name: str

    def predict(self, control_mean: np.ndarray, cell_type: str) -> np.ndarray:
        """Return predicted perturbed mean expression (shape: [n_genes])."""
        ...


def _safe_pearson(a: np.ndarray, b: np.ndarray) -> float:
    """Pearson r that returns 0.0 (not NaN) when a vector is constant.

    A constant prediction — e.g. the identity null's all-zero delta — has no
    linear relationship to anything; scoring it 0 is the honest reading, and it
    keeps the null from silently poisoning aggregate scores with NaN.
    """
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return 0.0
    r, _ = stats.pearsonr(a, b)
    return float(r) if np.isfinite(r) else 0.0


def score_delta(true_delta: np.ndarray, pred_delta: np.ndarray, top_k: int = 50) -> dict:
    """Score one predicted perturbation delta against the truth.

    Returns a dict of interpretable metrics:
      delta_pearson : linear agreement of predicted vs true delta (headline)
      delta_mse     : mean squared error on the delta
      delta_cosine  : direction agreement, magnitude-invariant
      deg_precision_at_k : of the K genes that truly move most, the fraction the
                           model also ranks in its top-K movers (did it find the
                           interferon-stimulated genes, not just the average shift?)
    """
    true_delta = np.asarray(true_delta, dtype=float).ravel()
    pred_delta = np.asarray(pred_delta, dtype=float).ravel()
    if true_delta.shape != pred_delta.shape:
        raise ValueError(f"shape mismatch: true {true_delta.shape} vs pred {pred_delta.shape}")

    mse = float(np.mean((true_delta - pred_delta) ** 2))

    denom = np.linalg.norm(true_delta) * np.linalg.norm(pred_delta)
    cosine = float(np.dot(true_delta, pred_delta) / denom) if denom > 1e-12 else 0.0

    k = min(top_k, true_delta.size)
    true_top = set(np.argsort(np.abs(true_delta))[-k:])
    pred_top = set(np.argsort(np.abs(pred_delta))[-k:])
    deg_prec = len(true_top & pred_top) / k if k > 0 else 0.0

    return {
        "delta_pearson": _safe_pearson(true_delta, pred_delta),
        "delta_mse": mse,
        "delta_cosine": cosine,
        "deg_precision_at_k": deg_prec,
    }


@dataclass
class CellTypeProfile:
    """Per-cell-type control/perturbed mean-expression, and the true delta."""

    cell_type: str
    control_mean: np.ndarray
    perturbed_mean: np.ndarray

    @property
    def true_delta(self) -> np.ndarray:
        return self.perturbed_mean - self.control_mean


class PerturbationBenchmark:
    """A backtest: known control -> known perturbed, scored per cell type.

    Build it from paired mean-expression matrices (cell_type x gene) for the
    control and perturbed conditions. `evaluate(model)` returns a per-cell-type
    score table; `report(model)` adds the null comparison that decides whether
    the model has actually earned trust.
    """

    def __init__(
        self,
        gene_names: list[str],
        control_means: dict[str, np.ndarray],
        perturbed_means: dict[str, np.ndarray],
        reliability: dict[str, float] | None = None,
    ) -> None:
        shared = sorted(set(control_means) & set(perturbed_means))
        if not shared:
            raise ValueError("no cell types shared between control and perturbed")
        self.gene_names = list(gene_names)
        self.profiles: dict[str, CellTypeProfile] = {
            ct: CellTypeProfile(ct, np.asarray(control_means[ct], float),
                                np.asarray(perturbed_means[ct], float))
            for ct in shared
        }
        # Per-cell-type reliability of the MEASURED delta (see `ceiling`).
        self.reliability: dict[str, float] = dict(reliability or {})

    # --- how much of the score is even reachable -----------------------------
    # A score means nothing without the ceiling it is scored against. The delta
    # we measure carries sampling noise, badly so where a cell type is rare, and
    # no model can correlate with noise. `reliability` is the split-half
    # reliability of the full-sample delta (see data/kang.py); its square root is
    # the highest delta_pearson any model could achieve against that measurement.
    LOW_RELIABILITY = 0.5   # below this, a fold is scoring noise, not biology

    def ceiling(self, cell_type: str) -> float | None:
        r = self.reliability.get(cell_type)
        return float(np.sqrt(max(r, 0.0))) if r is not None else None

    @property
    def unreliable_cell_types(self) -> list[str]:
        """Folds whose measured delta is too noisy to score anything against."""
        return [ct for ct in self.profiles
                if self.reliability.get(ct, 1.0) < self.LOW_RELIABILITY]

    @property
    def cell_types(self) -> list[str]:
        return list(self.profiles)

    def global_mean_delta(self) -> np.ndarray:
        """Average true delta across cell types — the shared perturbation signal.

        This is what the GlobalMeanShift null predicts for every cell type. A
        model must beat it to claim it captured cell-type-*specific* response.
        """
        return np.mean([p.true_delta for p in self.profiles.values()], axis=0)

    def evaluate(self, model: PerturbationModel) -> pd.DataFrame:
        rows = []
        for ct, prof in self.profiles.items():
            pred_perturbed = model.predict(prof.control_mean, ct)
            pred_delta = np.asarray(pred_perturbed, float).ravel() - prof.control_mean
            scores = score_delta(prof.true_delta, pred_delta)
            row = {"cell_type": ct, "model": model.name, **scores}
            ceil = self.ceiling(ct)
            if ceil is not None:
                # Score in context: 0.48 against a ceiling of 0.33 is a noise
                # fold, not a failure; 0.91 against a ceiling of 0.99 is a miss.
                row["ceiling"] = ceil
                row["frac_of_ceiling"] = (scores["delta_pearson"] / ceil
                                          if ceil > 1e-9 else float("nan"))
            rows.append(row)
        return pd.DataFrame(rows).set_index("cell_type")

    def report(self, model: PerturbationModel, nulls: list[PerturbationModel]) -> dict:
        """Score the model and each null; decide if the model beats the nulls.

        The headline verdict is on `delta_pearson`, aggregated (mean) across
        cell types. `beats_all_nulls` is the honest gate: pass it, and the model
        has demonstrated something a trivial predictor cannot.
        """
        table = self.evaluate(model)
        model_score = table["delta_pearson"].mean()
        null_scores = {n.name: self.evaluate(n)["delta_pearson"].mean() for n in nulls}
        out = {
            "model": model.name,
            "model_delta_pearson": float(model_score),
            "null_delta_pearson": {k: float(v) for k, v in null_scores.items()},
            "beats_all_nulls": bool(all(model_score > v + 1e-9 for v in null_scores.values())),
        }
        if self.reliability:
            noisy = self.unreliable_cell_types
            out["ceiling"] = float(np.mean([self.ceiling(ct) for ct in self.profiles]))
            out["unreliable_cell_types"] = noisy
            if noisy:
                # The aggregate silently averages in folds nothing can score.
                # Report the restricted number too, so a comparison is not
                # dominated by which model got luckier on noise.
                keep = [ct for ct in self.profiles if ct not in noisy]
                out["model_delta_pearson_reliable_only"] = float(
                    table.loc[keep, "delta_pearson"].mean())
                out["ceiling_reliable_only"] = float(
                    np.mean([self.ceiling(ct) for ct in keep]))
        return out
