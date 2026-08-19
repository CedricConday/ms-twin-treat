"""Deck figures — honest, validated palette. Static PNGs for a slide deck.

  fig1_backtest.png — the cell-model backtest on Kang IFN-β (aggregate, LOCTO):
                      we now beat the null (control-sim 0.87 > bar 0.85), and
                      scGPT (0.87) does NOT earn it (≈ np.corrcoef, p=0.55).
  fig2_arms.png     — mean lesion_proxy per real trial arm. Worked arms separate
                      from untreated; the harmed arm (APL) reads like untreated
                      because the model can't flag it yet.

All numbers verified in RESULTS.md. Nothing here is evidence about MS.
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


# --- Figure 1: the cell-model backtest, aggregate --------------------------
def fig_backtest():
    # verified aggregate delta_pearson (LOCTO) on Kang IFN-β; sg's scorecard.
    # (label, score, color, note)
    rows = [
        ("noise ceiling", 0.9075, BASE, "attainable limit"),
        ("control-similarity transfer", 0.8732, BLUE, "beats the bar — winner"),
        ("scGPT-blood embeddings", 0.8696, MUTED, "≈ free correlation (p=0.55)"),
        ("global-mean-shift null (the bar)", 0.8498, ORANGE, "the bar to beat"),
        ("affine transfer (v0)", 0.8204, BASE, "below the bar"),
        ("identity null", 0.0000, BASE, "captures nothing"),
    ]
    labels = [r[0] for r in rows]
    vals = [r[1] for r in rows]
    colors = [r[2] for r in rows]
    notes = [r[3] for r in rows]
    y = range(len(rows))

    fig, ax = plt.subplots(figsize=(10.6, 4.8))
    ax.barh(list(y), vals, height=0.62, color=colors)
    ax.invert_yaxis()
    ax.set_yticks(list(y)); ax.set_yticklabels(labels)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("aggregate delta_pearson  (higher = better; 1.0 = perfect)")
    ax.axvline(0.8498, color=MUTED, lw=1, ls=(0, (4, 3)))
    for i, (v, n) in enumerate(zip(vals, notes)):
        ax.text(v + 0.008, i, f"{v:.3f}   {n}", va="center", color=SEC, fontsize=8.5)
    ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
    _style(ax)
    fig.tight_layout(rect=[0, 0, 1, 0.85])
    fig.text(0.5, 0.965, "The cell model now beats the null — and scGPT didn't earn it",
             ha="center", va="top", fontsize=13, color=INK, weight="bold")
    fig.text(0.5, 0.905, "Kang IFN-β · leave-one-cell-type-out · scGPT clears the bar but ties free correlation (p=0.55)",
             ha="center", va="top", fontsize=9, color=SEC)
    p = os.path.join(OUT, "fig1_backtest.png")
    fig.savefig(p, dpi=200, facecolor=SURFACE); plt.close(fig)
    return p


# --- Figure 2: arm separation, and the gap we don't hide -------------------
def fig_arms():
    arms = [
        ("untreated", 32.7, MUTED, ""),
        ("IFN-β", 23.7, GOOD, "(−9.0)"),
        ("glatiramer acetate", 24.8, GOOD, "(−7.9)"),
        ("APL CGP77116", 36.1, CRIT, "(+3.4)"),
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
    fig.text(0.5, 0.90, "Separates the therapies that worked — and now flags the harmful arm (mechanism-encoded, not discovered)",
             ha="center", va="top", fontsize=9, color=SEC)
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
