"""Cell counts painted onto the brain: each Allen area filled by its count, N printed on it.

A choropleth, not a smoothed density. Every pixel of an area gets the same colour, so the
number you read is the exact count for that area - unlike fig1/fig5 where the smoothing
kernel mixes neighbouring cells across borders.

Outputs -> figures_v2/
  fig13_brain_counts_facepoke.png
  fig14_brain_counts_melittin.png
"""
import csv
import os

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures_v2")
CNT = os.path.join(HERE, "counts")

A = sio.loadmat(os.path.join(HERE, "allen_dorsal_map.mat"))
LABEL = A["label_img"]
BOUND = A["boundary"].astype(bool)
BREGMA = A["bregma_rc"].ravel() + 1
ACR = [str(x[0]) for x in A["acronyms"].ravel()]
IDS = A["ids"].ravel()
ACR2ID = {a: int(i) for i, a in zip(IDS, ACR)}
CTX = {"MOs", "MOp", "SSp-m", "SSp-n", "SSp-bfd", "SSp-ul", "SSp-ll", "SSp-tr", "SSp-un",
       "SSs", "VISp", "VISa", "VISam", "VISal", "VISpm", "VISrl", "VISl", "VISpor",
       "RSPd", "RSPv", "RSPagl", "AUDp", "AUDd", "AUDv", "ACAd", "TEa", "VISli", "AUDpo"}
CTX_MASK = np.isin(LABEL, [ACR2ID[a] for a in CTX if a in ACR2ID])
XLIM, YLIM = (60, 420), (400, 60)


def read(p):
    with open(os.path.join(CNT, p), newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        for k in ("n_cells", "n_excited", "n_inhibited", "n_nonsig"):
            r[k] = int(r[k])
    return rows


leaf = read("counts_by_area_leaf.csv")


MIN_LABEL_PX = 900     # below this an area is too small to hold readable text


def label_pos(area):
    """Medoid of the area in the LEFT hemisphere only.

    Counts are pooled over both hemispheres, so the same number belongs on both sides.
    Printing it twice doubled the clutter and made the lateral areas unreadable, so only
    one side is labelled.
    """
    if area not in ACR2ID:
        return []
    mid = int(round(BREGMA[1]))
    m = np.zeros(LABEL.shape, bool)
    m[:, :mid] = LABEL[:, :mid] == ACR2ID[area]
    if m.sum() < MIN_LABEL_PX:
        return []
    rr, cc = np.nonzero(m)
    d = (rr - rr.mean()) ** 2 + (cc - cc.mean()) ** 2
    k = int(np.argmin(d))
    return [(cc[k], rr[k])]


def paint(ax, vals, title, cmap, cbl, diverging=False, fmt="{:,.0f}", show_acr=True):
    """vals: {area: number}. Fill each area with its value."""
    img = np.full(LABEL.shape, np.nan)
    for a, v in vals.items():
        if a in ACR2ID:
            img[LABEL == ACR2ID[a]] = v
    ax.imshow(np.where(CTX_MASK, .93, np.nan), cmap="gray", vmin=0, vmax=1)
    if diverging:
        lim = np.nanmax(np.abs(list(vals.values())))
        im = ax.imshow(img, cmap=cmap, norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim))
    else:
        im = ax.imshow(img, cmap=cmap, vmin=0, vmax=np.nanmax(list(vals.values())))
    ax.imshow(np.ma.masked_where(~BOUND, BOUND), cmap="gray_r", alpha=.45,
              interpolation="nearest")
    hi = np.nanmax(np.abs(list(vals.values())))
    for a, v in vals.items():
        for (x, y) in label_pos(a):
            frac = abs(v) / max(hi, 1e-9)
            col = "white" if (not diverging and frac > .55) else "black"
            txt = (a + "\n" if show_acr else "") + fmt.format(v)
            ax.text(x, y, txt, ha="center", va="center", fontsize=7.6, color=col,
                    linespacing=1.05, fontweight="bold")
    ax.plot(BREGMA[1], BREGMA[0], "k+", ms=9, mew=1.5)
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=11, fontweight="bold")
    cb = plt.colorbar(im, ax=ax, fraction=.045, pad=.02)
    cb.set_label(cbl, fontsize=8)
    cb.ax.tick_params(labelsize=7)


def build(cond):
    rows = [r for r in leaf if r["condition"] == cond and r["area"] in CTX]
    return ({r["area"]: r["n_cells"] for r in rows},
            {r["area"]: r["n_excited"] for r in rows},
            {r["area"]: r["n_inhibited"] for r in rows},
            {r["area"]: r["n_excited"] - r["n_inhibited"] for r in rows})


# ============================================================ FIG 13 face poke
n_, e_, i_, d_ = build("Face poke")
tot = sum(n_.values())
fig, axes = plt.subplots(2, 2, figsize=(15, 14))
paint(axes[0, 0], n_, "Cells sampled", "Greys", "cells")
paint(axes[0, 1], e_, "Excited cells", "Reds", "cells")
paint(axes[1, 0], i_, "Inhibited cells", "Blues", "cells")
paint(axes[1, 1], d_, "Excited minus inhibited", "bwr", "cells", diverging=True,
      fmt="{:+,.0f}")
fig.suptitle(f"Fig 13  |  FACE POKE cell counts on the brain    {tot:,} cells, "
             f"{sum(e_.values()):,} excited, {sum(i_.values()):,} inhibited, 4 mice\n"
             "each Allen area filled by its own count; anterior up, bregma marked",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .93])
fig.savefig(os.path.join(OUT, "fig13_brain_counts_facepoke.png"), dpi=160)
plt.close(fig)
print("fig13 done")

# ============================================================ FIG 14 melittin
fig, axes = plt.subplots(2, 4, figsize=(24, 13))
for row, cond in enumerate(["Melittin", "Melittin+morphine"]):
    n_, e_, i_, d_ = build(cond)
    paint(axes[row, 0], n_, f"{cond}\nCells sampled", "Greys", "cells")
    paint(axes[row, 1], e_, f"{cond}\nUp vs baseline", "Reds", "cells")
    paint(axes[row, 2], i_, f"{cond}\nDown vs baseline", "Blues", "cells")
    paint(axes[row, 3], d_, f"{cond}\nUp minus down", "bwr", "cells",
          diverging=True, fmt="{:+,.0f}")
fig.suptitle("Fig 14  |  MELITTIN cell counts on the brain    "
             "melittin 4 mice, melittin + morphine 2 mice\n"
             "each Allen area filled by its own count; the two rows sampled different "
             "numbers of cells, so compare the up-minus-down column, not the raw totals",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .92])
fig.savefig(os.path.join(OUT, "fig14_brain_counts_melittin.png"), dpi=150)
plt.close(fig)
print("fig14 done")
print("->", OUT)
