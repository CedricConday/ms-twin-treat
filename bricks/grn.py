"""B3 — Gene-regulatory network brick (v0).

Infers a co-expression network from the real Kang expression matrix. Two paths:

  DEFAULT (Stage, demo-safe): Spearman co-expression among the top-variance
  genes, thresholded to the strongest edges. Real inference from real data,
  runs in seconds, deterministic, no cluster to hang on. Edges are *undirected*
  co-expression -- a true GRN needs directionality (TF->target) from motif
  enrichment, so these are labelled method="spearman-coexpression", validated=False.

  UPGRADE (opt-in, off the demo path): grnboost2_edges() uses arboreto/pyscenic
  (tree-based, directed from a TF list). Heavier (dask); call it explicitly when
  you want the real SCENIC-style network, not in the auto demo.

Writes state["grn_edges"] = list[(gene_a, gene_b, weight)] + meta.
"""

from __future__ import annotations

import numpy as np
import scanpy as sc
from scipy.stats import rankdata

from data.kang import fetch_kang


def _load_expr(n_genes: int, n_cells: int, seed: int = 0):
    """Return (gene_names, X[n_cells, n_genes]) of top-variance genes, log-normalized."""
    adata = fetch_kang()
    if float(adata.X.max()) > 30:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    var = X.var(axis=0)
    top = np.argsort(var)[-n_genes:]
    genes = [str(g) for g in np.asarray(adata.var_names)[top]]
    Xg = X[:, top]
    if Xg.shape[0] > n_cells:
        rng = np.random.default_rng(seed)
        Xg = Xg[rng.choice(Xg.shape[0], n_cells, replace=False)]
    return genes, Xg


def infer_grn(n_genes: int = 100, n_cells: int = 2000, top_edges: int = 200,
              seed: int = 0) -> list[tuple[str, str, float]]:
    """Spearman co-expression edges among top-variance genes (undirected)."""
    genes, Xg = _load_expr(n_genes, n_cells, seed)
    ranks = np.apply_along_axis(rankdata, 0, Xg)          # rank per gene -> Spearman via Pearson
    corr = np.corrcoef(ranks, rowvar=False)
    corr = np.nan_to_num(corr)
    iu = np.triu_indices_from(corr, k=1)                  # upper triangle, no self-edges
    weights = corr[iu]
    order = np.argsort(np.abs(weights))[::-1][:top_edges]
    return [(genes[iu[0][k]], genes[iu[1][k]], float(weights[k])) for k in order]


def grnboost2_edges(n_genes: int = 100, n_cells: int = 2000, seed: int = 0):
    """UPGRADE path: directed GRN via arboreto GRNBoost2. Heavier; opt-in only."""
    import pandas as pd
    from arboreto.algo import grnboost2
    genes, Xg = _load_expr(n_genes, n_cells, seed)
    df = pd.DataFrame(Xg, columns=genes)
    net = grnboost2(expression_data=df, seed=seed)
    return [(r.TF, r.target, float(r.importance)) for r in net.itertuples()]


class GRNBrick:
    """Stage: infer a co-expression network from Kang, write grn_edges."""

    name = "grn:spearman-coexpression"

    def __init__(self, n_genes: int = 100, top_edges: int = 200) -> None:
        self.n_genes = n_genes
        self.top_edges = top_edges

    def __call__(self, state: dict) -> dict:
        edges = infer_grn(n_genes=self.n_genes, top_edges=self.top_edges)
        state["grn_edges"] = edges
        state["grn_meta"] = {"validated": False, "method": self.name,
                             "n_edges": len(edges), "directed": False}
        return state


if __name__ == "__main__":
    print("B3 GRN — inferring co-expression network from Kang IFN-beta data...")
    edges = infer_grn(n_genes=100, top_edges=200)
    print(f"  inferred {len(edges)} edges among top-variance genes")
    print("  strongest 8:")
    for a, b, w in edges[:8]:
        print(f"    {a:<12} -- {b:<12}  r={w:+.3f}")
    print("  (undirected co-expression; validated=False; GRNBoost2 upgrade available)")
