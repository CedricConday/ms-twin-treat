"""Kang 2018 (GSE96583) loader — the v0 backtest anchor.

~14k human PBMCs, 8 donors, +/- IFN-beta stimulation, cell-type labeled.
IFN-beta is a first-line MS disease-modifying therapy, so this is the closest
open "real MS drug on immune cells" dataset we have. The known biology — a
strong, well-characterized interferon-stimulated-gene response — is what the
harness backtests a model against.

Source: figshare file 34464122 (the processed AnnData used by pertpy.data.kang_2018).
Cached under data/cache/ (gitignored — no datasets committed to the repo).

This module turns raw counts into the object the harness consumes: per-cell-type
mean expression for control and for IFN-beta, on log-normalized data.
"""

from __future__ import annotations

import json
import os

import numpy as np
import scanpy as sc
from anndata import AnnData

from backtest.harness import PerturbationBenchmark

_FIGSHARE_URL = "https://ndownloader.figshare.com/files/34464122"
_CACHE = os.path.join(os.path.dirname(__file__), "cache")
_H5AD = os.path.join(_CACHE, "kang_2018.h5ad")
_RELIABILITY = os.path.join(_CACHE, "kang_reliability.json")


def _load_reliability_cache() -> dict[str, float]:
    try:
        with open(_RELIABILITY) as fh:
            return {str(k): float(v) for k, v in json.load(fh).items()}
    except (OSError, ValueError):
        return {}


def _save_reliability_cache(reliability: dict[str, float]) -> None:
    try:
        os.makedirs(_CACHE, exist_ok=True)
        with open(_RELIABILITY, "w") as fh:
            json.dump(reliability, fh, indent=2, sort_keys=True)
    except OSError:
        pass    # a missing cache costs time on the next run, nothing more


def fetch_kang() -> AnnData:
    """Download (once, cached) and return the Kang AnnData."""
    os.makedirs(_CACHE, exist_ok=True)
    adata = sc.read(_H5AD, backup_url=_FIGSHARE_URL)
    return adata


def _detect_col(adata: AnnData, candidates: list[str], values_hint: set[str]) -> str:
    """Find the obs column holding a known field, by name or by its values."""
    for c in candidates:
        if c in adata.obs.columns:
            return c
    # Fall back to any obs column whose value set overlaps the hint.
    for c in adata.obs.columns:
        vals = {str(v).lower() for v in adata.obs[c].unique()[:20]}
        if vals & values_hint:
            return c
    raise KeyError(f"could not locate a column among {candidates} or matching {values_hint}; "
                   f"available: {list(adata.obs.columns)}")


def _normalize(adata: AnnData) -> AnnData:
    """Log-normalize if the matrix still looks like raw counts."""
    x = adata.X
    xmax = float(x.max())
    looks_like_counts = xmax > 30 and np.allclose(x.data if hasattr(x, "data") else x,
                                                  np.round(x.data if hasattr(x, "data") else x))
    if looks_like_counts:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    return adata


def _split_half_reliability(x_ctrl: np.ndarray, x_stim: np.ndarray,
                            n_rep: int, rng: np.random.Generator) -> float:
    """Reliability of this cell type's MEASURED delta, in [0, 1].

    Split the control cells in half and the stimulated cells in half, compute a
    delta from each half, correlate the two. That is the reliability of a
    half-sample delta; Spearman-Brown lifts it to the full sample. The square
    root of the result is the best delta_pearson any model could score here —
    see PerturbationBenchmark.ceiling.

    This exists because an aggregate score silently averages folds that cannot
    be scored at all. On Kang, Megakaryocytes (63/69 cells) land near zero.
    """
    rs = []
    for _ in range(n_rep):
        ca, cb = np.array_split(rng.permutation(len(x_ctrl)), 2)
        sa, sb = np.array_split(rng.permutation(len(x_stim)), 2)
        da = x_stim[sa].mean(axis=0) - x_ctrl[ca].mean(axis=0)
        db = x_stim[sb].mean(axis=0) - x_ctrl[cb].mean(axis=0)
        if np.std(da) < 1e-12 or np.std(db) < 1e-12:
            continue
        rs.append(float(np.corrcoef(da, db)[0, 1]))
    if not rs:
        return 0.0
    r_half = float(np.mean(rs))
    return float(np.clip(2 * r_half / (1 + r_half), 0.0, 1.0))   # Spearman-Brown


def to_benchmark(adata: AnnData | None = None, min_cells: int = 50,
                 with_reliability: bool = True, n_rep: int = 20,
                 seed: int = 0) -> PerturbationBenchmark:
    """Build a PerturbationBenchmark: per-cell-type control vs IFN-beta means.

    Cell types with fewer than `min_cells` in either condition are dropped — a
    mean over a handful of cells is noise, and the honest move is to not score
    what we cannot estimate.

    `with_reliability` also estimates, per cell type, how much of the measured
    delta is real signal rather than sampling noise, so every score can be read
    against the ceiling it is scored against. It is cached to
    data/cache/kang_reliability.json, so the cost is paid once. Turn it off only
    if you genuinely want a score with no ceiling attached.
    """
    if adata is None:
        adata = fetch_kang()
    adata = _normalize(adata.copy())

    cond_col = _detect_col(adata, ["label", "condition", "stim"],
                           {"ctrl", "stim", "control", "stimulated"})
    ct_col = _detect_col(adata, ["cell_type", "cell", "celltype", "cell_abbr"], set())

    vals = {str(v).lower() for v in adata.obs[cond_col].unique()}
    ctrl_key = next(v for v in adata.obs[cond_col].unique() if str(v).lower() in {"ctrl", "control"})
    stim_key = next(v for v in adata.obs[cond_col].unique()
                    if str(v).lower() in {"stim", "stimulated"})

    gene_names = list(adata.var_names)
    cached = _load_reliability_cache() if with_reliability else {}
    rng = np.random.default_rng(seed)

    control_means, perturbed_means, reliability = {}, {}, {}
    for ct in adata.obs[ct_col].unique():
        ct_mask = adata.obs[ct_col] == ct
        ctrl = adata[ct_mask & (adata.obs[cond_col] == ctrl_key)]
        stim = adata[ct_mask & (adata.obs[cond_col] == stim_key)]
        if ctrl.n_obs < min_cells or stim.n_obs < min_cells:
            continue
        name = str(ct)
        control_means[name] = np.asarray(ctrl.X.mean(axis=0)).ravel()
        perturbed_means[name] = np.asarray(stim.X.mean(axis=0)).ravel()
        if with_reliability:
            if name in cached:
                reliability[name] = cached[name]
            else:
                dense = lambda a: np.asarray(  # noqa: E731
                    a.todense() if hasattr(a, "todense") else a, dtype=np.float64)
                reliability[name] = _split_half_reliability(
                    dense(ctrl.X), dense(stim.X), n_rep, rng)

    if not control_means:
        raise RuntimeError("no cell type had enough cells in both conditions")
    if with_reliability and reliability != cached:
        _save_reliability_cache(reliability)
    return PerturbationBenchmark(gene_names, control_means, perturbed_means,
                                 reliability=reliability or None)
