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

import os

import numpy as np
import scanpy as sc
from anndata import AnnData

from backtest.harness import PerturbationBenchmark

_FIGSHARE_URL = "https://ndownloader.figshare.com/files/34464122"
_CACHE = os.path.join(os.path.dirname(__file__), "cache")
_H5AD = os.path.join(_CACHE, "kang_2018.h5ad")


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


def to_benchmark(adata: AnnData | None = None, min_cells: int = 50) -> PerturbationBenchmark:
    """Build a PerturbationBenchmark: per-cell-type control vs IFN-beta means.

    Cell types with fewer than `min_cells` in either condition are dropped — a
    mean over a handful of cells is noise, and the honest move is to not score
    what we cannot estimate.
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
    control_means, perturbed_means = {}, {}
    for ct in adata.obs[ct_col].unique():
        ct_mask = adata.obs[ct_col] == ct
        ctrl = adata[ct_mask & (adata.obs[cond_col] == ctrl_key)]
        stim = adata[ct_mask & (adata.obs[cond_col] == stim_key)]
        if ctrl.n_obs < min_cells or stim.n_obs < min_cells:
            continue
        control_means[str(ct)] = np.asarray(ctrl.X.mean(axis=0)).ravel()
        perturbed_means[str(ct)] = np.asarray(stim.X.mean(axis=0)).ravel()

    if not control_means:
        raise RuntimeError("no cell type had enough cells in both conditions")
    return PerturbationBenchmark(gene_names, control_means, perturbed_means)
