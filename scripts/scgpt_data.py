"""Kang -> scGPT input tensors: log-normalize, per-cell quantile binning (51
bins), non-zero-gene tokenization with a prepended <cls>, exactly as scGPT's
Preprocessor + tokenize_and_pad_batch do."""
from __future__ import annotations
import numpy as np, scanpy as sc, anndata as ad

N_BINS, MAX_LEN, MASK_VALUE, PAD_VALUE, CLS_VALUE = 51, 1200, -1, -2, 0.0


def load_kang(path="data/cache/kang_2018.h5ad"):
    a = ad.read_h5ad(path)
    x = a.X
    xmax = float(x.max())
    data = x.data if hasattr(x, "data") else x
    if xmax > 30 and np.allclose(data, np.round(data)):
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    return a


def _digitize(x, bins, rng):
    left = np.digitize(x, bins)
    right = np.digitize(x, bins, right=True)
    d = rng.random(len(x)) * (right - left) + left
    return np.ceil(d).astype(np.int64)


def bin_rows(X, rng):
    """Per-cell quantile binning of non-zero values into 1..N_BINS-1."""
    X = np.asarray(X.todense() if hasattr(X, "todense") else X, dtype=np.float32)
    out = np.zeros_like(X, dtype=np.float32)
    for i in range(X.shape[0]):
        nz = np.nonzero(X[i])[0]
        if nz.size == 0:
            continue
        vals = X[i, nz]
        bins = np.quantile(vals, np.linspace(0, 1, N_BINS - 1))
        out[i, nz] = _digitize(vals, bins, rng).astype(np.float32)
    return out


def tokenize(binned, gene_ids, cls_id, pad_id, rng, max_len=MAX_LEN):
    """Non-zero genes only, <cls> first, random-truncate/pad to max_len."""
    n = binned.shape[0]
    gi = np.full((n, max_len), pad_id, dtype=np.int64)
    gv = np.full((n, max_len), PAD_VALUE, dtype=np.float32)
    for i in range(n):
        idx = np.nonzero(binned[i])[0]
        idx = idx[gene_ids[idx] >= 0]                     # keep in-vocab genes
        if idx.size > max_len - 1:
            idx = rng.choice(idx, max_len - 1, replace=False)
        L = idx.size + 1
        gi[i, 0] = cls_id
        gv[i, 0] = CLS_VALUE
        gi[i, 1:L] = gene_ids[idx]
        gv[i, 1:L] = binned[i, idx]
    return gi, gv


def map_genes(var_names, vocab):
    """Kang gene symbol -> scGPT vocab id, -1 if absent."""
    return np.array([vocab.get(str(g), -1) for g in var_names], dtype=np.int64)
