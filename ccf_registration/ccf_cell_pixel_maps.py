"""Cell-resolution maps in CCF space - the aligned version of the raw take-home map.

Two kinds of panel, neither of them smoothed:

  cell scatter   every cell drawn at its own CCF position, coloured by response amplitude.
                 This is the direct counterpart of the raw "Mean map after alignment"
                 figure, but in Allen coordinates with area borders on top.
  pixel counts   cells binned into 250 um squares, so each square shows how many cells sit
                 there. Answers "how many cells per pixel" without any kernel.

Outputs -> figures_v2/
  fig15_cellmap_facepoke.png
  fig16_cellmap_melittin.png
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

A = sio.loadmat(os.path.join(HERE, "allen_dorsal_map.mat"))
LABEL = A["label_img"]
BOUND = A["boundary"].astype(bool)
RES = float(A["res_um"][0][0])          # 25 um per atlas pixel
BREGMA = A["bregma_rc"].ravel() + 1
ACR = [str(x[0]) for x in A["acronyms"].ravel()]
IDS = A["ids"].ravel()
CTX = {"MOs", "MOp", "SSp-m", "SSp-n", "SSp-bfd", "SSp-ul", "SSp-ll", "SSp-tr", "SSp-un",
       "SSs", "VISp", "VISa", "VISam", "VISal", "VISpm", "VISrl", "VISl", "VISpor",
       "RSPd", "RSPv", "RSPagl", "AUDp", "AUDd", "AUDv", "ACAd", "TEa", "VISli", "AUDpo"}
CTX_MASK = np.isin(LABEL, [int(i) for i, a in zip(IDS, ACR) if a in CTX])
XLIM, YLIM = (60, 420), (400, 60)
BIN_UM = 250
BIN = int(round(BIN_UM / RES))          # 10 atlas px


def load(path, ampcol="amp"):
    out = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            out.append((float(r["atlas_col"]), float(r["atlas_row"]),
                        float(r[ampcol]), r["is_exc"] == "1", r["is_inh"] == "1",
                        r.get("cohort", "")))
    return np.array([(x, y, a) for x, y, a, _, _, _ in out]), \
           np.array([e for _, _, _, e, _, _ in out]), \
           np.array([i for _, _, _, _, i, _ in out]), \
           np.array([c for _, _, _, _, _, c in out], dtype=object)


def backdrop(ax):
    ax.imshow(np.where(CTX_MASK, .96, np.nan), cmap="gray", vmin=0, vmax=1)
    ax.imshow(np.ma.masked_where(~BOUND, BOUND), cmap="gray_r", alpha=.55,
              interpolation="nearest")
    ax.plot(BREGMA[1], BREGMA[0], "k+", ms=10, mew=1.6)
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(.5, .975, "A", transform=ax.transAxes, ha="center", color="#C0392B",
            fontsize=15, fontweight="bold")
    ax.text(.5, .008, "P", transform=ax.transAxes, ha="center", color="#2C6FBB",
            fontsize=15, fontweight="bold")


def cell_scatter(ax, XYA, vlim, title, cbl):
    backdrop(ax)
    # draw the weakest first so strong responders are not buried
    o = np.argsort(np.abs(XYA[:, 2]))
    sc = ax.scatter(XYA[o, 0], XYA[o, 1], c=XYA[o, 2], cmap="bwr",
                    norm=TwoSlopeNorm(vcenter=0, vmin=-vlim, vmax=vlim),
                    s=1.6, linewidths=0, alpha=.85)
    ax.set_title(title, fontsize=11, fontweight="bold")
    cb = plt.colorbar(sc, ax=ax, fraction=.045, pad=.02)
    cb.set_label(cbl, fontsize=8)
    cb.ax.tick_params(labelsize=7)


def pixel_counts(ax, XYA, mask, title, cmap, annotate=False):
    """Bin cells into BIN_UM squares and show the count in each."""
    h, w = LABEL.shape
    nr, nc = h // BIN + 1, w // BIN + 1
    br = (XYA[mask, 1] // BIN).astype(int)
    bc = (XYA[mask, 0] // BIN).astype(int)
    ok = (br >= 0) & (br < nr) & (bc >= 0) & (bc < nc)
    grid = np.zeros((nr, nc))
    np.add.at(grid, (br[ok], bc[ok]), 1)
    # only show bins that overlap cortex
    ctx_small = np.zeros((nr, nc), bool)
    rr, cc = np.nonzero(CTX_MASK)
    ctx_small[rr // BIN, cc // BIN] = True
    grid[~ctx_small] = np.nan
    grid[grid == 0] = np.nan

    big = np.kron(grid, np.ones((BIN, BIN)))[:h, :w]
    ax.imshow(np.where(CTX_MASK, .96, np.nan), cmap="gray", vmin=0, vmax=1)
    im = ax.imshow(big, cmap=cmap, vmin=0, vmax=np.nanpercentile(grid, 99))
    ax.imshow(np.ma.masked_where(~BOUND, BOUND), cmap="gray_r", alpha=.35,
              interpolation="nearest")
    ax.plot(BREGMA[1], BREGMA[0], "k+", ms=10, mew=1.6)
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.set_xticks([]); ax.set_yticks([])
    if annotate:
        hi = np.nanpercentile(grid, 99)
        for r in range(nr):
            for c in range(nc):
                v = grid[r, c]
                if np.isfinite(v) and v >= 1:
                    x, y = c * BIN + BIN / 2, r * BIN + BIN / 2
                    if XLIM[0] < x < XLIM[1] and YLIM[1] < y < YLIM[0]:
                        ax.text(x, y, f"{int(v)}", ha="center", va="center", fontsize=4.2,
                                color="white" if v > .6 * hi else "black")
    tot = int(np.nansum(grid))
    ax.set_title(f"{title}\n{tot:,} cells in {int(np.isfinite(grid).sum())} bins",
                 fontsize=11, fontweight="bold")
    cb = plt.colorbar(im, ax=ax, fraction=.045, pad=.02)
    cb.set_label(f"cells per {BIN_UM} um bin", fontsize=8)
    cb.ax.tick_params(labelsize=7)


# ============================================================ FIG 15 face poke
XYA, exc, inh, _ = load(os.path.join(HERE, "n18", "cell_ccf_all.csv"))
fig, axes = plt.subplots(2, 2, figsize=(15, 14))
cell_scatter(axes[0, 0], XYA, 0.3,
             f"Every cell, coloured by poke amplitude\n{len(XYA):,} cells, 18 recordings, 4 mice",
             "poke amplitude (z)")
pixel_counts(axes[0, 1], XYA, np.ones(len(XYA), bool), "All cells per bin", "viridis",
             annotate=True)
pixel_counts(axes[1, 0], XYA, exc, "Excited cells per bin", "Reds", annotate=True)
pixel_counts(axes[1, 1], XYA, inh, "Inhibited cells per bin", "Blues", annotate=True)
fig.suptitle("Fig 15  |  FACE POKE at cell resolution, in Allen coordinates\n"
             f"no smoothing anywhere: left = one dot per cell, others = counts in "
             f"{BIN_UM} um bins",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .93])
fig.savefig(os.path.join(OUT, "fig15_cellmap_facepoke.png"), dpi=170)
plt.close(fig)
print("fig15 done")

# ============================================================ FIG 16 melittin
XYA, exc, inh, coh = load(os.path.join(HERE, "melittin_ccf", "cell_ccf_all.csv"), "amp_pct")
only = coh != "melittin_morphine"
morp = ~only
fig, axes = plt.subplots(2, 3, figsize=(19, 13))
for row, (m, name) in enumerate([(only, "Melittin only  (8 rec, 4 mice)"),
                                 (morp, "Melittin + morphine  (4 rec, 2 mice)")]):
    cell_scatter(axes[row, 0], XYA[m], 100,
                 f"{name}\nevery cell, {int(m.sum()):,} total",
                 "change vs baseline (%)")
    pixel_counts(axes[row, 1], XYA, m & exc, f"{name}\nup cells per bin", "Reds")
    pixel_counts(axes[row, 2], XYA, m & inh, f"{name}\ndown cells per bin", "Blues")
fig.suptitle("Fig 16  |  MELITTIN at cell resolution, in Allen coordinates\n"
             f"one dot per cell, then counts in {BIN_UM} um bins; the two rows sampled "
             "different numbers of cells",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .92])
fig.savefig(os.path.join(OUT, "fig16_cellmap_melittin.png"), dpi=160)
plt.close(fig)
print("fig16 done")
print("->", OUT)
