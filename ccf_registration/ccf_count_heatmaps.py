"""Cell-count heatmaps: area x condition and area x mouse, with the raw N printed in cell.

Each column gets its own colour normalisation. A shared scale would be meaningless because
the columns hold different things (total cells 8,682 vs excited 696) and the conditions
sampled different numbers of cells (43k face poke vs 80k melittin), so colour encodes rank
within a column and the printed number is the value you actually read.

Outputs -> figures_v2/
  fig11_count_heatmap_by_condition.png
  fig12_count_heatmap_by_mouse.png
"""
import collections
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures_v2")
CNT = os.path.join(HERE, "counts")
CONDS = ["Face poke", "Melittin", "Melittin+morphine"]


def read(p):
    with open(os.path.join(CNT, p), newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in ("n_cells", "n_excited", "n_inhibited", "n_nonsig"):
            r[k] = int(r[k])
    return rows


leaf = read("counts_by_area_leaf.csv")
pm = read("counts_by_area_and_mouse.csv")

# areas present in all three conditions with a usable sample
by = {(r["condition"], r["area"]): r for r in leaf}
areas = sorted({r["area"] for r in leaf if r["n_cells"] >= 300},
               key=lambda a: -by.get(("Face poke", a), {"n_cells": 0})["n_cells"])
areas = [a for a in areas if all((c, a) in by for c in CONDS)]


def grid(ax, M, rowlab, collab, cmap, title, fmt="{:,.0f}", fs=7.5):
    """One heatmap, each COLUMN normalised independently."""
    Z = np.full(M.shape, np.nan)
    for c in range(M.shape[1]):
        col = M[:, c].astype(float)
        lo, hi = np.nanmin(col), np.nanmax(col)
        Z[:, c] = (col - lo) / max(hi - lo, 1e-9)
    ax.imshow(Z, cmap=cmap, aspect="auto", vmin=0, vmax=1)
    for r in range(M.shape[0]):
        for c in range(M.shape[1]):
            if np.isfinite(M[r, c]):
                ax.text(c, r, fmt.format(M[r, c]), ha="center", va="center", fontsize=fs,
                        color="white" if Z[r, c] > .62 else "black")
    ax.set_xticks(range(len(collab)))
    ax.set_xticklabels(collab, fontsize=8.5, rotation=30, ha="right")
    ax.set_yticks(range(len(rowlab)))
    ax.set_yticklabels(rowlab, fontsize=8.5)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)


# ============================================================ FIG 11
fig, axes = plt.subplots(1, 3, figsize=(17, 9), sharey=True)
for ax, cond in zip(axes, CONDS):
    M = np.array([[by[(cond, a)]["n_cells"], by[(cond, a)]["n_excited"],
                   by[(cond, a)]["n_inhibited"], by[(cond, a)]["n_nonsig"]] for a in areas],
                 float)
    n = M[:, 0].sum(); e = M[:, 1].sum(); i = M[:, 2].sum()
    up = "excited" if cond == "Face poke" else "up"
    dn = "inhibited" if cond == "Face poke" else "down"
    grid(ax, M, areas, ["cells", up, dn, "no change"], "viridis",
         f"{cond}\n{int(n):,} cells   {int(e):,} {up}   {int(i):,} {dn}")
axes[0].set_ylabel("Allen area  (ordered by face-poke cell count)", fontsize=10)
fig.suptitle("Fig 11  |  Cell counts, area x condition\n"
             "colour is normalised within each column; read the printed number, not the shade",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .90])
fig.savefig(os.path.join(OUT, "fig11_count_heatmap_by_condition.png"), dpi=160)
plt.close(fig)
print("fig11 done")

# ============================================================ FIG 12 : per mouse
groups = ["SSp", "VIS", "MO", "RSP", "SSs", "AUD"]
pmd = {(r["condition"], r["area_group"], r["mouse"]): r for r in pm}

fig, axes = plt.subplots(2, 3, figsize=(16, 8), gridspec_kw=dict(hspace=.45, wspace=.30))
for col, cond in enumerate(CONDS):
    mice = sorted({r["mouse"] for r in pm if r["condition"] == cond})
    for row, (key, cmap, name) in enumerate(
            [("n_excited", "Reds", "excited / up"), ("n_inhibited", "Blues", "inhibited / down")]):
        ax = axes[row, col]
        M = np.array([[pmd[(cond, g, m)][key] if (cond, g, m) in pmd else np.nan
                       for m in mice] for g in groups], float)
        grid(ax, M, groups, mice, cmap, f"{cond}\n{name}", fs=8)
fig.suptitle("Fig 12  |  Cell counts per mouse, area group x animal\n"
             "colour normalised within each animal's column; face poke m201 is large because "
             "it contributed 6 of the 18 recordings",
             fontsize=13.5, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .89])
fig.savefig(os.path.join(OUT, "fig12_count_heatmap_by_mouse.png"), dpi=160)
plt.close(fig)
print("fig12 done")
print("->", OUT)
