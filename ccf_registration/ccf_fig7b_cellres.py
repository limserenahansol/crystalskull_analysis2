"""Fig 7b rebuilt at cell resolution - no smoothing in the brain panels.

The previous version drew a block-binned, interpolated median surface. This one shows the
cells themselves.

A median cannot be drawn per cell: one cell has one amplitude, not a median. The
cell-resolution counterpart of the median is the SIGN of each cell, because the median is
negative exactly when more than half the cells are negative. So:

  panel A   every cell coloured by its amplitude          (the raw picture)
  panel B   every cell coloured by sign only              (the median, made visible)
  panel C   anterior vs posterior third, median per recording
  panel D   median amplitude by Allen area, anterior to posterior

Output -> figures_v2/fig7b_facepoke_median_cellres.png
"""
import collections
import csv
import os

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures_v2")

A = sio.loadmat(os.path.join(HERE, "allen_dorsal_map.mat"))
LABEL = A["label_img"]
BOUND = A["boundary"].astype(bool)
BREGMA = A["bregma_rc"].ravel() + 1
ACR = [str(x[0]) for x in A["acronyms"].ravel()]
IDS = A["ids"].ravel()
CTX = {"MOs", "MOp", "SSp-m", "SSp-n", "SSp-bfd", "SSp-ul", "SSp-ll", "SSp-tr", "SSp-un",
       "SSs", "VISp", "VISa", "VISam", "VISal", "VISpm", "VISrl", "VISl", "VISpor",
       "RSPd", "RSPv", "RSPagl", "AUDp", "AUDd", "AUDv", "ACAd", "TEa"}
CTX_MASK = np.isin(LABEL, [int(i) for i, a in zip(IDS, ACR) if a in CTX])
XLIM, YLIM = (60, 420), (400, 60)
RED, BLUE = "#C0392B", "#2C6FBB"

rows = list(csv.DictReader(open(os.path.join(HERE, "n18", "cell_ccf_all.csv"))))
x = np.array([float(r["atlas_col"]) for r in rows])
y = np.array([float(r["atlas_row"]) for r in rows])
amp = np.array([float(r["amp"]) for r in rows])
ap = np.array([float(r["AP_mm"]) for r in rows])
lab = np.array([r["label"] for r in rows], dtype=object)
area = np.array([r["area_acronym"] for r in rows], dtype=object)


def backdrop(ax):
    ax.imshow(np.where(CTX_MASK, .96, np.nan), cmap="gray", vmin=0, vmax=1)
    ax.imshow(np.ma.masked_where(~BOUND, BOUND), cmap="gray_r", alpha=.5,
              interpolation="nearest")
    ax.plot(BREGMA[1], BREGMA[0], "k+", ms=10, mew=1.6)
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.set_xticks([]); ax.set_yticks([])
    ax.text(.5, .975, "A", transform=ax.transAxes, ha="center", color=RED,
            fontsize=15, fontweight="bold")
    ax.text(.5, .008, "P", transform=ax.transAxes, ha="center", color=BLUE,
            fontsize=15, fontweight="bold")


fig = plt.figure(figsize=(14.5, 13.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=.24, wspace=.30)

# ---- A: every cell, amplitude
ax = fig.add_subplot(gs[0, 0])
backdrop(ax)
o = np.argsort(np.abs(amp))
sc = ax.scatter(x[o], y[o], c=amp[o], cmap="bwr",
                norm=TwoSlopeNorm(vcenter=0, vmin=-.3, vmax=.3),
                s=1.5, linewidths=0, alpha=.85)
ax.set_title(f"Every cell, poke amplitude\n{len(amp):,} cells, 18 recordings, 4 mice",
             fontsize=11, fontweight="bold")
cb = plt.colorbar(sc, ax=ax, fraction=.045, pad=.02)
cb.set_label("poke amplitude (z)", fontsize=8.5)
cb.ax.tick_params(labelsize=7)

# ---- B: every cell, sign only
ax = fig.add_subplot(gs[0, 1])
backdrop(ax)
pos = amp > 0
# Draw in random order. Plotting all blue then all red buries the blue and makes
# 44% red look like the majority - the opposite of what the panel is meant to show.
rng = np.random.default_rng(0)
sh = rng.permutation(len(amp))
ax.scatter(x[sh], y[sh], c=np.where(pos[sh], RED, BLUE), s=1.5, linewidths=0,
           alpha=.6)
ax.set_title("Every cell, sign only\n"
             f"{100*pos.mean():.0f}% above zero, {100*(~pos).mean():.0f}% below",
             fontsize=11, fontweight="bold")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=RED, label="amp > 0", ms=6),
                   Line2D([], [], marker="o", ls="", color=BLUE, label="amp < 0", ms=6)],
          loc="lower left", frameon=False, fontsize=8.5)

# ---- C: A vs P, median per recording
ax = fig.add_subplot(gs[1, 0])
post, ant = [], []
for L in sorted(set(lab)):
    m = lab == L
    lo, hi = np.percentile(ap[m], [33.3, 66.7])
    post.append(np.median(amp[m & (ap <= lo)]))
    ant.append(np.median(amp[m & (ap >= hi)]))
post, ant = np.array(post), np.array(ant)
_, pv = wilcoxon(ant, post)
for p_, a_ in zip(post, ant):
    ax.plot([0, 1], [p_, a_], "-", color="0.65", lw=.9, zorder=1)
ax.scatter(np.zeros_like(post), post, s=42, c=BLUE, zorder=3, edgecolor="w", lw=.6)
ax.scatter(np.ones_like(ant), ant, s=42, c=RED, zorder=3, edgecolor="w", lw=.6)
ax.plot([0, 1], [np.mean(post), np.mean(ant)], "k-", lw=3, zorder=4)
ax.axhline(0, color="k", ls=":", lw=1)
ax.set_xlim(-.35, 1.35); ax.set_xticks([0, 1])
ax.set_xticklabels(["Posterior\nthird", "Anterior\nthird"])
ax.set_ylabel("median poke amplitude (z)")
ax.set_title(f"{int((ant > post).sum())}/{len(ant)} recordings A > P\nsignrank p = {pv:.1e}",
             fontsize=11, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)

# ---- D: median by area
ax = fig.add_subplot(gs[1, 1])
byarea = collections.defaultdict(list)
for a_, ap_, am_ in zip(area, ap, amp):
    if a_ in CTX:
        byarea[a_].append((ap_, am_))
ar = [(np.mean([t[0] for t in v]), k, np.median([t[1] for t in v]), len(v))
      for k, v in byarea.items() if len(v) >= 300]
ar.sort(reverse=True)
yy = np.arange(len(ar))
vals = [t[2] for t in ar]
ax.barh(yy, vals, color=[RED if v > 0 else BLUE for v in vals], edgecolor="white")
ax.set_yticks(yy)
ax.set_yticklabels([f"{t[1]}  ({t[0]:+.1f})" for t in ar], fontsize=8.5)
ax.invert_yaxis()
ax.axvline(0, color="k", lw=1)
ax.set_xlabel("median poke amplitude (z)")
ax.set_title("By Allen area, anterior -> posterior\n(brackets = AP mm from bregma)",
             fontsize=11, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)

fig.suptitle("Fig 7b  |  Face poke, MEDIAN amplitude at cell resolution\n"
             "brain panels are the cells themselves, no smoothing",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .93])
fig.savefig(os.path.join(OUT, "fig7b_facepoke_median_cellres.png"), dpi=170)
plt.close(fig)

print(f"cells above zero: {100*pos.mean():.1f}%   median amp {np.median(amp):+.4f}")
print(f"A vs P (median): {int((ant > post).sum())}/{len(ant)}, "
      f"post {post.mean():+.4f}  ant {ant.mean():+.4f}, p={pv:.2e}")
print("->", os.path.join(OUT, "fig7b_facepoke_median_cellres.png"))
