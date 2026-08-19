"""Extract scGPT embeddings for the Kang backtest — runs in `.venv-scgpt` ONLY.

This is the heavy half of the B2-SCGPT lane. It runs the pretrained scGPT
encoder over Kang PBMCs and writes a small npz that the scoring brick
(`bricks/cell_scgpt.py`, main pinned `.venv`) reads with plain numpy. That split
is deliberate: torch/scanpy never touch the pinned analysis stack.

What it writes to data/cache/scgpt_embeddings_<ckpt>.npz:

  cell_types      (K,)      cell-type names, sorted
  cell_emb_ctrl   (K, 512)  scGPT <cls> embedding, averaged over CONTROL cells
                            only, L2-normalized. Control-only is what makes the
                            downstream model leave-one-cell-type-out honest: the
                            embedding of a held-out cell type never sees its
                            IFN-beta-perturbed data.
  n_cells         (K,)      how many control cells went into each centroid
  gene_names      (G,)      Kang gene order (matches the harness)
  gene_emb        (G, 512)  scGPT pretrained gene embedding per Kang gene,
                            L2-normalized; all-zero row for the ~12% of Kang
                            genes absent from scGPT's vocabulary
  gene_in_vocab   (G,) bool which genes actually got an embedding

Usage:
    . .venv-scgpt/bin/activate
    python scripts/scgpt_embed.py --ckpt data/cache/scgpt/blood --n-per-type 128
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scgpt_data as D  # noqa: E402
import scgpt_model as M  # noqa: E402


def cell_type_centroids(model, adata, gene_ids, vocab, n_per_type, seed, batch=8):
    """Mean L2-normalized <cls> embedding per cell type, over control cells."""
    rng = np.random.default_rng(seed)
    ctrl = adata[adata.obs["label"] == "ctrl"]
    cell_types = sorted(str(c) for c in ctrl.obs["cell_type"].unique())
    cents, counts = [], []
    t0 = time.time()
    for ct in cell_types:
        idx = np.where(ctrl.obs["cell_type"].values == ct)[0]
        pick = rng.choice(idx, min(n_per_type, len(idx)), replace=False)
        embs = []
        for s in range(0, len(pick), batch):
            X = ctrl[pick[s:s + batch]].X
            binned = D.bin_rows(X, rng)
            gi, gv = D.tokenize(binned, gene_ids, vocab["<cls>"], vocab["<pad>"], rng)
            with torch.no_grad():
                embs.append(model.cell_embeddings(torch.from_numpy(gi),
                                                  torch.from_numpy(gv)).numpy())
        e = np.vstack(embs)
        c = e.mean(axis=0)
        cents.append(c / max(np.linalg.norm(c), 1e-8))
        counts.append(len(e))
        print(f"  {ct:<20} n={len(e):>4}  {time.time() - t0:6.0f}s", flush=True)
    return cell_types, np.vstack(cents), np.array(counts)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", default="data/cache/scgpt/blood")
    ap.add_argument("--n-per-type", type=int, default=128)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-dir", default="data/cache")
    args = ap.parse_args()

    tag = os.path.basename(args.ckpt.rstrip("/"))
    torch.set_num_threads(max(1, (os.cpu_count() or 2)))

    model, vocab, ckpt_args, missing, unexpected = M.load(args.ckpt)
    assert not missing, f"checkpoint is missing encoder weights: {missing}"
    print(f"[{tag}] loaded scGPT: {ckpt_args['nlayers']}L d={ckpt_args['embsize']} "
          f"ntoken={ckpt_args['ntoken']}; unused heads={sorted(unexpected)}")

    adata = D.load_kang()
    gene_ids = D.map_genes(adata.var_names, vocab)
    in_vocab = gene_ids >= 0
    print(f"[{tag}] gene vocab coverage: {int(in_vocab.sum())}/{len(gene_ids)} Kang genes")

    # Gene embeddings: a straight lookup into the pretrained embedding table —
    # no forward pass, so nothing here can be distorted by the re-implementation.
    sd = torch.load(f"{args.ckpt}/best_model.pt", map_location="cpu", weights_only=True)
    table = sd["encoder.embedding.weight"].numpy()
    gene_emb = np.zeros((len(gene_ids), table.shape[1]), dtype=np.float32)
    gene_emb[in_vocab] = table[gene_ids[in_vocab]]
    norms = np.linalg.norm(gene_emb, axis=1, keepdims=True)
    gene_emb = np.divide(gene_emb, np.maximum(norms, 1e-8))

    print(f"[{tag}] encoding control cells (<={args.n_per_type} per cell type)...")
    cell_types, cell_emb, n_cells = cell_type_centroids(
        model, adata, gene_ids, vocab, args.n_per_type, args.seed)

    os.makedirs(args.out_dir, exist_ok=True)
    out = os.path.join(args.out_dir, f"scgpt_embeddings_{tag}.npz")
    np.savez_compressed(
        out,
        cell_types=np.array(cell_types),
        cell_emb_ctrl=cell_emb.astype(np.float32),
        n_cells=n_cells,
        gene_names=np.array([str(g) for g in adata.var_names]),
        gene_emb=gene_emb.astype(np.float32),
        gene_in_vocab=in_vocab,
        checkpoint=np.array(tag),
        n_per_type=np.array(args.n_per_type),
    )
    print(f"[{tag}] wrote {out}")


if __name__ == "__main__":
    main()
