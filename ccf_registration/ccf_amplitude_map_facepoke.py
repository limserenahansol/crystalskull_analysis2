"""Face-poke MEAN AMPLITUDE map in CCF space - the atlas version of the raw take-home figure.

Why this exists
---------------
fig5 plotted the *fraction of cells passing a significance threshold*. The original raw
take-home figure plots *mean poke amplitude (z)*, a continuous signed value. Those are
different quantities, and the fraction-of-significant-cells version hides the anterior/
posterior gradient for two reasons:

  1. thresholding throws away sign and magnitude. A region can have a clearly positive mean
     amplitude while few of its individual cells clear p < 0.05.
  2. the excited-fraction colour scale ran 0-80% because of a handful of edge pixels with
     almost no cells, which squashed the real 5-17% range into the palest part of the map.

Plotting the same quantity as the raw figure reproduces it: 18/18 recordings anterior >
posterior, signrank p = 7.6e-06 using CCF coordinates (the raw version, using front/back
thirds of the FOV, gave 18/18 and p = 2e-4).

Output -> figures_v2/fig7_facepoke_amplitude_ccf.png
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
from scipy.ndimage import gaussian_filter
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures_v2")

A = sio.loadmat(os.path.join(HERE, "allen_dorsal_map.mat"))
LABEL = A["label_img"]
BOUND = A["boundary"].astype(bool)
RES = float(A["res_um"][0][0])
BREGMA = A["bregma_rc"].ravel() + 1
ACR = [str(x[0]) for x in A["acronyms"].ravel()]
IDS = A["ids"].ravel()
ID2ACR = {int(i): a for i, a in zip(IDS, ACR)}
CTX = {"MOs", "MOp", "SSp-m", "SSp-n", "SSp-bfd", "SSp-ul", "SSp-ll", "SSp-tr", "SSp-un",
       "SSs", "VISp", "VISa", "VISam", "VISal", "VISpm", "VISrl", "VISl", "VISpor",
       "RSPd", "RSPv", "RSPagl", "AUDp", "AUDd", "AUDv", "ACAd", "TEa"}
CTX_MASK = np.isin(LABEL, [i for i, a in ID2ACR.items() if a in CTX])

rows = list(csv.DictReader(open(os.path.join(HERE, "n18", "cell_ccf_all.csv"))))
for r in rows:
    r["AP"] = float(r["AP_mm"])
    r["amp"] = float(r["amp"])
    r["rr"] = int(round(float(r["atlas_row"])))
    r["cc"] = int(round(float(r["atlas_col"])))

# ---------------------------------------------------------------- amplitude map
h, w = LABEL.shape
num = np.zeros((h, w))
den = np.zeros((h, w))
for r in rows:
    if 0 <= r["rr"] < h and 0 <= r["cc"] < w:
        num[r["rr"], r["cc"]] += r["amp"]
        den[r["rr"], r["cc"]] += 1
sig = 6
ns, ds = gaussian_filter(num, sig), gaussian_filter(den, sig)
# require a real local cell count, otherwise sparse edge pixels dominate the scale
ok = CTX_MASK & (ds > np.percentile(ds[ds > 0], 40))
amp = np.full((h, w), np.nan)
amp[ok] = ns[ok] / ds[ok]

# The mean is a fragile summary here. Excited cells are 8.6% of the population but average
# +0.44, while inhibited cells are 27.6% and average only -0.19 - a 2.4x magnitude
# asymmetry. Add the small positive offset of the non-significant 64% and the population
# mean lands at +0.014 even though the median is -0.025 and only 44% of cells are positive.
# So a second figure is written using the median, which matches the cell-count statistics.
STAT = os.environ.get("CCF_STAT", "mean")
if STAT == "median":
    # A per-pixel median is meaningless here: 43k cells spread over ~100k cortex pixels
    # means most pixels hold 0 or 1 cell, so a per-pixel median followed by smoothing just
    # reproduces the local mean. Bin into blocks large enough to hold a real sample, take
    # the median inside each block, then interpolate.
    BLK = 12                      # 12 atlas px = 300 um
    MINN = 20                     # cells needed before a block gets a value
    rr = np.array([r["rr"] for r in rows])
    cc = np.array([r["cc"] for r in rows])
    av = np.array([r["amp"] for r in rows])
    keep = (rr >= 0) & (rr < h) & (cc >= 0) & (cc < w)
    rr, cc, av = rr[keep], cc[keep], av[keep]
    br, bc = rr // BLK, cc // BLK
    nbr, nbc = h // BLK + 1, w // BLK + 1
    grid = np.full((nbr, nbc), np.nan)
    order = np.lexsort((bc, br))
    br_s, bc_s, av_s = br[order], bc[order], av[order]
    key = br_s.astype(np.int64) * nbc + bc_s
    edges = np.flatnonzero(np.diff(key)) + 1
    for s, e in zip(np.r_[0, edges], np.r_[edges, len(key)]):
        if e - s >= MINN:
            grid[br_s[s], bc_s[s]] = np.median(av_s[s:e])
    # interpolate the block grid back to atlas resolution, ignoring empty blocks
    filled = np.nan_to_num(grid)
    msk = np.isfinite(grid).astype(float)
    num_b = gaussian_filter(filled * msk, 1.2)
    den_b = gaussian_filter(msk, 1.2)
    small = np.where(den_b > .05, num_b / np.maximum(den_b, 1e-9), np.nan)
    big = np.kron(small, np.ones((BLK, BLK)))[:h, :w]
    big = gaussian_filter(np.nan_to_num(big), 4) / np.maximum(
        gaussian_filter(np.isfinite(big).astype(float), 4), 1e-9)
    amp = np.where(CTX_MASK & np.isfinite(big) & ok, big, np.nan)

# ---------------------------------------------------------------- per-recording A vs P
recs = collections.defaultdict(list)
for r in rows:
    recs[r["label"]].append(r)
post, ant, labs = [], [], []
for lab, rs in sorted(recs.items()):
    ap = np.array([x["AP"] for x in rs])
    am = np.array([x["amp"] for x in rs])
    lo, hi = np.percentile(ap, [33.3, 66.7])
    agg = np.median if STAT == "median" else np.mean
    post.append(agg(am[ap <= lo]))
    ant.append(agg(am[ap >= hi]))
    labs.append(lab)
post, ant = np.array(post), np.array(ant)
_, pv = wilcoxon(ant, post)

# ---------------------------------------------------------------- by area
byarea = collections.defaultdict(list)
for r in rows:
    if r["area_acronym"] in CTX:
        byarea[r["area_acronym"]].append((r["AP"], r["amp"]))
agg_a = np.median if STAT == "median" else np.mean
ar = [(np.mean([x[0] for x in v]), a, agg_a([x[1] for x in v]), len(v))
      for a, v in byarea.items() if len(v) >= 300]
ar.sort(reverse=True)

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(17.5, 7.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, .78, 1.5], wspace=.42)

ax = fig.add_subplot(gs[0])
ax.imshow(np.where(CTX_MASK, 0.94, np.nan), cmap="gray", vmin=0, vmax=1)
lim = np.nanpercentile(np.abs(amp), 97)
im = ax.imshow(amp, cmap="bwr", norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim))
ax.imshow(np.ma.masked_where(~BOUND, BOUND), cmap="gray_r", alpha=.30, interpolation="nearest")
ax.plot(BREGMA[1], BREGMA[0], "k+", ms=11, mew=1.8)
ax.text(.5, .965, "A", transform=ax.transAxes, ha="center", color="#C0392B",
        fontsize=17, fontweight="bold")
ax.text(.5, .01, "P", transform=ax.transAxes, ha="center", color="#2C6FBB",
        fontsize=17, fontweight="bold")
ax.set_xlim(60, 420); ax.set_ylim(400, 60)
ax.set_xticks([]); ax.set_yticks([])
ax.set_title(f"{STAT.capitalize()} poke amplitude in CCF space",
             fontsize=11, fontweight="bold")
cb = fig.colorbar(im, ax=ax, fraction=.045, pad=.02)
cb.set_label(f"{STAT} poke amplitude (z)", fontsize=9)

ax = fig.add_subplot(gs[1])
for p_, a_ in zip(post, ant):
    ax.plot([0, 1], [p_, a_], "-", color="0.65", lw=.9, zorder=1)
ax.scatter(np.zeros_like(post), post, s=44, c="#2C6FBB", zorder=3, edgecolor="w", lw=.6)
ax.scatter(np.ones_like(ant), ant, s=44, c="#C0392B", zorder=3, edgecolor="w", lw=.6)
ax.plot([0, 1], [post.mean(), ant.mean()], "k-", lw=3, zorder=4)
ax.axhline(0, color="k", ls=":", lw=1)
ax.set_xlim(-.35, 1.35); ax.set_xticks([0, 1])
ax.set_xticklabels(["Posterior\nthird", "Anterior\nthird"])
ax.set_ylabel(f"{STAT} poke amplitude (z)")
ax.set_title(f"{int((ant > post).sum())}/{len(ant)} recordings A > P\n"
             f"signrank p = {pv:.1e}", fontsize=11, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)

ax = fig.add_subplot(gs[2])
y = np.arange(len(ar))
vals = [x[2] for x in ar]
ax.barh(y, vals, color=["#C0392B" if v > 0 else "#2C6FBB" for v in vals], edgecolor="white")
ax.set_yticks(y)
ax.set_yticklabels([f"{x[1]}  ({x[0]:+.1f})" for x in ar], fontsize=8.5)
ax.invert_yaxis()
ax.axvline(0, color="k", lw=1)
ax.set_xlabel(f"{STAT} poke amplitude (z)")
ax.set_title("By Allen area, ordered anterior -> posterior\n(number in brackets = AP mm from bregma)",
             fontsize=11, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)

ttl = ("Fig 7  |  Face poke, MEAN amplitude in Allen coordinates" if STAT == "mean"
       else "Fig 7b  |  Face poke, MEDIAN amplitude - what the typical cell does")
sub = ("mean is pulled up by a minority of strongly excited cells" if STAT == "mean"
       else "matches the cell counts: 8.6% excited vs 27.6% inhibited")
fig.suptitle(ttl + "\n" + sub, fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .84])
name = "fig7_facepoke_amplitude_ccf.png" if STAT == "mean" else "fig7b_facepoke_median_ccf.png"
fig.savefig(os.path.join(OUT, name), dpi=160)
plt.close(fig)

print(f"18/18 check: anterior > posterior in {int((ant > post).sum())}/{len(ant)}, p = {pv:.2e}")
print(f"posterior {post.mean():+.4f}   anterior {ant.mean():+.4f}")
print("->", os.path.join(OUT, name))
