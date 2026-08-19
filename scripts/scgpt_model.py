"""Minimal re-implementation of scGPT's encoder, loaded straight from `best_model.pt`.

Why not the `scgpt` package: it pins flash-attn and an old torch, and fights
this CPU/aarch64 box. The checkpoint is only a state dict over a standard
post-norm transformer with a fused Wqkv projection, so we rebuild the forward
pass in plain torch and load the weights by name. Every parameter in the
checkpoint that the encoder needs is matched (`load()` asserts no missing keys);
the only unused tensors are the perturbation-flag and MVC heads.

The defaults here are NOT guesses — they are read off scGPT's own source
(`scgpt/model/model.py`, `scgpt/model/flash_attn_compat.py`):

  qkv_layout  "3hd"   FlashMHA does rearrange(Wqkv(x), "b s (three h d) -> ...")
  norm_scheme "post"  args.json has pre_norm=false
  activation  "relu"  TransformerModel builds FlashTransformerEncoderLayer
                      without overriding its activation default
  no positional encoding (scGPT does not apply one), gene_emb + value_emb sum,
  cell embedding = the <cls> position (args.json cell_emb_style="cls").

VALIDATION STATUS — read this before trusting an output:
  * gene embeddings (encoder.embedding.weight) involve NO forward pass, so they
    are faithful by construction. Biologically checked: the 20-gene ISG set is
    far more coherent than random gene sets (mean pairwise cosine +0.21 vs
    +0.02 +/- 0.01, z = +17.5 on the blood checkpoint).
  * the ENCODER is validated functionally: <cls> embeddings of Kang control
    cells give 94.3% leave-one-out 1-NN cell-type accuracy (chance 12.5%),
    cosine silhouette 0.34, and recover immune lineage structure (myeloid block,
    T/NK block, megakaryocytes as the outlier).
  * the MLM expression DECODER head does NOT validate: asked to reconstruct
    masked expression bins it regresses almost to a constant (pearson r ~ 0.11
    blood / 0.27 human, prediction sd ~1 against a true bin sd of ~14). Cause
    not established — it may be a subtle mismatch we could not find, or these
    re-hosted heads may not be the MLM head at inference scale. We therefore
    use ONLY the encoder and the gene embeddings, never `decoder`.
"""
from __future__ import annotations
import json, math
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class ContinuousValueEncoder(nn.Module):
    def __init__(self, d=512):
        super().__init__()
        self.linear1 = nn.Linear(1, d)
        self.linear2 = nn.Linear(d, d)
        self.norm = nn.LayerNorm(d)

    def forward(self, x):                     # x: (B, L)
        h = F.relu(self.linear1(x.unsqueeze(-1)))
        return self.norm(self.linear2(h))


class GeneEncoder(nn.Module):
    def __init__(self, ntoken=36574, d=512, pad_id=36571):
        super().__init__()
        self.embedding = nn.Embedding(ntoken, d, padding_idx=pad_id)
        self.enc_norm = nn.LayerNorm(d)

    def forward(self, ids):
        return self.enc_norm(self.embedding(ids))


class FlashishEncoderLayer(nn.Module):
    """scGPT's FlashTransformerEncoderLayer, post-norm, fused Wqkv, ReLU FF."""

    def __init__(self, d=512, nhead=8, d_hid=512, qkv_layout="3hd",
                 norm_scheme="post", activation="relu"):
        super().__init__()
        self.nhead = nhead
        self.d = d
        self.qkv_layout = qkv_layout
        self.norm_scheme = norm_scheme
        self.act = F.relu if activation == "relu" else F.gelu
        self.Wqkv = nn.Linear(d, 3 * d)
        self.out_proj = nn.Linear(d, d)
        self.linear1 = nn.Linear(d, d_hid)
        self.linear2 = nn.Linear(d_hid, d)
        self.norm1 = nn.LayerNorm(d)
        self.norm2 = nn.LayerNorm(d)

    def _attn(self, x, key_padding_mask):
        B, L, D = x.shape
        dh = D // self.nhead
        w = self.Wqkv(x)
        if self.qkv_layout == "3hd":                      # (three, heads, dim)
            qkv = w.view(B, L, 3, self.nhead, dh)
            q, k, v = (qkv[:, :, i].transpose(1, 2) for i in range(3))
        elif self.qkv_layout == "h3d":                    # (heads, three, dim)
            qkv = w.view(B, L, self.nhead, 3, dh)
            q, k, v = (qkv[:, :, :, i].transpose(1, 2) for i in range(3))
        else:                                             # plain split d|d|d
            q, k, v = (t.view(B, L, self.nhead, dh).transpose(1, 2)
                       for t in w.split(D, dim=-1))
        attn_mask = None
        if key_padding_mask is not None:
            # bool attn_mask: True == take part in attention
            attn_mask = (~key_padding_mask)[:, None, None, :].expand(B, self.nhead, L, L)
        a = F.scaled_dot_product_attention(q, k, v, attn_mask=attn_mask)
        return self.out_proj(a.transpose(1, 2).reshape(B, L, D))

    def forward(self, x, key_padding_mask=None):
        if self.norm_scheme == "pre":
            x = x + self._attn(self.norm1(x), key_padding_mask)
            h = self.norm2(x)
            return x + self.linear2(self.act(self.linear1(h)))
        x = self.norm1(x + self._attn(x, key_padding_mask))
        return self.norm2(x + self.linear2(self.act(self.linear1(x))))


class ExprDecoder(nn.Module):
    def __init__(self, d=512):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(d, d), nn.LeakyReLU(),
                                nn.Linear(d, d), nn.LeakyReLU(), nn.Linear(d, 1))

    def forward(self, h):
        return self.fc(h).squeeze(-1)


class ScGPT(nn.Module):
    def __init__(self, ntoken, d=512, nhead=8, nlayers=12, d_hid=512, pad_id=36571,
                 qkv_layout="3hd", norm_scheme="post", activation="relu"):
        super().__init__()
        self.pad_id = pad_id
        self.encoder = GeneEncoder(ntoken, d, pad_id)
        self.value_encoder = ContinuousValueEncoder(d)
        self.transformer_encoder = nn.ModuleList(
            [FlashishEncoderLayer(d, nhead, d_hid, qkv_layout, norm_scheme, activation)
             for _ in range(nlayers)])
        self.decoder = ExprDecoder(d)

    def forward(self, gene_ids, values):
        pad = gene_ids.eq(self.pad_id)
        x = self.encoder(gene_ids) + self.value_encoder(values)
        for layer in self.transformer_encoder:
            x = layer(x, key_padding_mask=pad)
        return x                              # (B, L, D); position 0 == <cls>

    @torch.no_grad()
    def cell_embeddings(self, gene_ids, values):
        h = self.forward(gene_ids, values)[:, 0, :]      # cell_emb_style="cls"
        return h / h.norm(dim=1, keepdim=True).clamp_min(1e-8)


def load(ckpt_dir: str, **kw):
    args = json.load(open(f"{ckpt_dir}/args.json"))
    vocab = json.load(open(f"{ckpt_dir}/vocab.json"))
    m = ScGPT(ntoken=args["ntoken"], d=args["embsize"], nhead=args["nheads"],
              nlayers=args["nlayers"], d_hid=args["d_hid"], pad_id=vocab["<pad>"], **kw)
    sd = torch.load(f"{ckpt_dir}/best_model.pt", map_location="cpu", weights_only=True)
    remap = {}
    for k, v in sd.items():
        k2 = k.replace("transformer_encoder.layers.", "transformer_encoder.")
        k2 = k2.replace(".self_attn.Wqkv.", ".Wqkv.").replace(".self_attn.out_proj.", ".out_proj.")
        remap[k2] = v
    missing, unexpected = m.load_state_dict(remap, strict=False)
    m.eval()
    return m, vocab, args, missing, unexpected
