"""Deck figures — honest, validated palette. Static PNGs for a slide deck.

Two figures, both light-surface, colorblind-safe (dataviz reference palette):
  fig1_backtest.png  — per-cell-type: the mean-shift bar vs our cell model.
                       We fall just short everywhere; Megakaryocytes are the tell.
  fig2_arms.png      — mean lesion_proxy per real trial arm. Worked arms separate
                       from untreated (status green); the harmed arm (APL, status
                       red) reads like untreated because the model can't flag it yet.

All numbers are verified in RESULTS.md. Nothing here is evidence about MS.
Run:  python -m results.figures   ->  writes results/figures/*.png
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# --- dataviz reference palette (validated) ---------------------------------
SURFACE = "#fcfcfb"; INK = "#0b0b0b"; SEC = "#52514e"; MUTED = "#898781"
GRID = "#e1e0d9"; BASE = "#c3c2b7"
BLUE = "#2a78d6"; ORANGE = "#eb6834"; GOOD = "#0ca30c"; CRIT = "#d03b3b"

OUT = os.path.join(os.path.dirname(__file__), "figures")
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "font.family": "sans-serif", "font.size": 11,
    "text.color": INK, "axes.labelcolor": SEC, "axes.edgecolor": BASE,
    "xtick.color": MUTED, "ytick.color": SEC,
})


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(BASE)
    ax.spines["bottom"].set_color(BASE)
    ax.tick_params(length=0)


# --- Figure 1: the honest backtest -----------------------------------------
def fig_backtest():
    # verified per-cell-type delta_pearson, sorted by the bar (mean-shift) desc
    rows = [
        ("Dendritic cells", 0.932, 0.899), ("FCGR3A+ Mono", 0.930, 0.900),
        ("CD14+ Mono", 0.912, 0.870), ("B cells", 0.905, 0.881),
        ("CD4 T", 0.888, 0.868), ("NK", 0.887, 0.864),
        ("CD8 T", 0.870, 0.850), ("Megakaryocytes", 0.476, 0.432),
    ]
    labels = [r[0] for r in rows]
    bar = [r[1] for r in rows]
    model = [r[2] for r in rows]
    y = range(len(rows))
    h = 0.38

    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    ax.barh([i + h / 2 + 0.02 for i in y], bar, height=h, color=ORANGE, label="mean-shift null (the bar)")
    ax.barh([i - h / 2 - 0.02 for i in y], model, height=h, color=BLUE, label="our cell model (LOCTO)")
    ax.set_yticks(list(y)); ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("delta_pearson  (higher = better; 1.0 = perfect)")
    ax.axvline(0.85, color=MUTED, lw=1, ls=(0, (4, 3)))
    ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="lower right", fontsize=9.5)
    _style(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.86])
    fig.text(0.5, 0.965, "Backtest: our cell model vs the null it must beat",
             ha="center", va="top", fontsize=13, color=INK, weight="bold")
    fig.text(0.5, 0.905, "Kang 2018 IFN-β PBMCs · aggregate 0.82 model vs 0.85 bar — does NOT beat it yet (honest)",
             ha="center", va="top", fontsize=9.5, color=SEC)
    p = os.path.join(OUT, "fig1_backtest.png")
    fig.savefig(p, dpi=200, facecolor=SURFACE); plt.close(fig)
    return p


# --- Figure 2: arm separation, and the gap we don't hide -------------------
def fig_arms():
    # verified mean lesion_proxy per arm (results/experiment.py)
    arms = [
        ("untreated", 32.7, MUTED, ""),
        ("IFN-β", 23.7, GOOD, "(−9.0)"),
        ("glatiramer acetate", 24.8, GOOD, "(−7.9)"),
        ("APL CGP77116", 32.7, CRIT, ""),
    ]
    labels = [a[0] for a in arms]
    vals = [a[1] for a in arms]
    colors = [a[2] for a in arms]
    notes = [a[3] for a in arms]
    y = range(len(arms))

    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    ax.barh(list(y), vals, height=0.6, color=colors)
    ax.invert_yaxis()
    ax.set_yticks(list(y)); ax.set_yticklabels(labels)
    ax.set_xlim(0, 40)
    ax.set_xlabel("mean lesion proxy across the virtual cohort  (lower = better; proxy, not clinical)")
    ax.axvline(32.7, color=BASE, lw=1, ls=(0, (4, 3)))
    for i, (v, n) in enumerate(zip(vals, notes)):
        ax.text(v + 0.6, i, f"{v:.1f}  {n}".rstrip(), va="center", color=SEC, fontsize=10)
    ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
    legend = [Patch(color=GOOD, label="worked (approved DMT)"),
              Patch(color=CRIT, label="harmed (Phase II halted)"),
              Patch(color=MUTED, label="untreated control")]
    _style(ax)
    fig.tight_layout(rect=[0, 0.10, 1, 0.84])
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.5, 0.015))
    fig.text(0.5, 0.965, "Virtual cohort across real trial arms",
             ha="center", va="top", fontsize=13, color=INK, weight="bold")
    fig.text(0.5, 0.90, "The pipeline separates the therapies that worked — and openly can't yet catch the one that harmed patients",
             ha="center", va="top", fontsize=9.5, color=SEC)
    p = os.path.join(OUT, "fig2_arms.png")
    fig.savefig(p, dpi=200, facecolor=SURFACE); plt.close(fig)
    return p


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    print("wrote:", fig_backtest())
    print("wrote:", fig_arms())
    print("Static PNGs, light surface, validated colorblind-safe palette. Proxies, not clinical endpoints.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
