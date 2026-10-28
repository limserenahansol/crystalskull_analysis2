"""Standalone CCF heatmaps, one figure per experiment.

Face poke and melittin are different measurements and must not share a colour scale:
  face poke  peri-poke vs pre-poke, stimulus-locked, same cells within a session
  melittin   injection session vs baseline session, independent EXTRACT IDs
A sustained session-level shift is far easier to detect than a phasic response, so the
melittin percentages run ~3x higher for reasons that have nothing to do with biology.
Each figure therefore gets its own scale, set from its own data.

Outputs -> figures_v2/
  fig5_heatmap_facepoke.png
  fig6_heatmap_melittin.png
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
XLIM, YLIM = (60, 420), (400, 60)
MIN_CELLS = 300


def load(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


poke = load(os.path.join(HERE, "n18", "cell_ccf_all.csv"))
for r in poke:
    mv = r["movie"]
    r["mouse"] = "m1" if mv.startswith("0004749770_m1_") else mv.split("_")[0]
mel = load(os.path.join(HERE, "melittin_ccf", "cell_ccf_all.csv"))


def maps(rows, sigma=6):
    """Smoothed % excited, % inhibited and their difference, in atlas space."""
    h, w = LABEL.shape
    tot = np.zeros((h, w)); e = np.zeros((h, w)); i = np.zeros((h, w))
    for r in rows:
        y, x = int(round(float(r["atlas_row"]))), int(round(float(r["atlas_col"])))
        if not (0 <= y < h and 0 <= x < w):
            continue
        tot[y, x] += 1
        if r["is_exc"] == "1":
            e[y, x] += 1
        if r["is_inh"] == "1":
            i[y, x] += 1
    ts, es, is_ = (gaussian_filter(z, sigma) for z in (tot, e, i))
    ok = CTX_MASK & (ts > np.percentile(ts[ts > 0], 25))
    pe = np.full((h, w), np.nan); pi = np.full((h, w), np.nan)
    pe[ok] = 100 * es[ok] / ts[ok]
    pi[ok] = 100 * is_[ok] / ts[ok]
    return pe, pi, pe - pi


def draw(ax, m, cmap, title, norm=None, vmin=None, vmax=None, cbl="% of local cells"):
    ax.imshow(np.where(CTX_MASK, 0.94, np.nan), cmap="gray", vmin=0, vmax=1)
    im = ax.imshow(m, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax)
    ax.imshow(np.ma.masked_where(~BOUND, BOUND), cmap="gray_r", alpha=.30,
              interpolation="nearest")
    ax.plot(BREGMA[1], BREGMA[0], "k+", ms=10, mew=1.6)
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=10.5, fontweight="bold")
    cb = plt.colorbar(im, ax=ax, fraction=.045, pad=.02)
    cb.set_label(cbl, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return im


def region_table(rows):
    tot = collections.Counter(); e = collections.Counter(); i = collections.Counter()
    for r in rows:
        a = r["area_acronym"]
        if a not in CTX:
            continue
        tot[a] += 1
        if r["is_exc"] == "1":
            e[a] += 1
        if r["is_inh"] == "1":
            i[a] += 1
    return tot, e, i


def region_strip(axes, cols, areas, title, fmt="{:.1f}"):
    """Areas along x, one row per metric, each row on its OWN colour scale.

    A single shared scale across rows is misleading here: % excited runs 5-17, % inhibited
    18-36 and the net -30 to -5, so one scale renders the excited row uniformly pale.
    """
    for ax, (name, d, cmap, diverging) in zip(axes, cols):
        v = np.array([d.get(a, np.nan) for a in areas], float)
        if diverging:
            lim = np.nanmax(np.abs(v))
            norm = TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim)
            im = ax.imshow(v[None, :], cmap=cmap, aspect="auto", norm=norm)
        else:
            im = ax.imshow(v[None, :], cmap=cmap, aspect="auto")
        for c, val in enumerate(v):
            if np.isfinite(val):
                lo, hi = np.nanmin(v), np.nanmax(v)
                rel = (val - lo) / max(hi - lo, 1e-9)
                ax.text(c, 0, fmt.format(val), ha="center", va="center", fontsize=7,
                        color="white" if rel > 0.65 else "black", rotation=90)
        ax.set_yticks([0]); ax.set_yticklabels([name], fontsize=9, fontweight="bold")
        ax.set_xticks(range(len(areas)))
        ax.set_xticklabels(areas if ax is axes[-1] else [""] * len(areas),
                           fontsize=8, rotation=70, ha="right")
        plt.colorbar(im, ax=ax, fraction=.012, pad=.006).ax.tick_params(labelsize=6)
    axes[0].set_title(title, fontsize=11, fontweight="bold", pad=8)


# =============================================================== FIG 5 : face poke
pe, pi, pn = maps(poke)
tot, e, i = region_table(poke)
areas = sorted([a for a, n in tot.items() if n >= MIN_CELLS],
               key=lambda a: -(e[a] / tot[a]))
n_mice = len(set(r["mouse"] for r in poke))
n_rec = len(set(r["label"] for r in poke))

fig = plt.figure(figsize=(15, 10.5))
gs = fig.add_gridspec(2, 3, height_ratios=[2.5, 1.35], hspace=.30, wspace=.18)
draw(fig.add_subplot(gs[0, 0]), pe, "Reds", "Excited\n(peri-poke > pre-poke)")
draw(fig.add_subplot(gs[0, 1]), pi, "Blues", "Inhibited\n(peri-poke < pre-poke)")
lim = np.nanpercentile(np.abs(pn), 98)
draw(fig.add_subplot(gs[0, 2]), pn, "RdBu_r", "Net (excited - inhibited)",
     norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim), cbl="percentage points")
sub = gs[1, :].subgridspec(3, 1, hspace=.12)
region_strip([fig.add_subplot(sub[k]) for k in range(3)],
             [("% excited", {a: 100 * e[a] / tot[a] for a in areas}, "Reds", False),
              ("% inhibited", {a: 100 * i[a] / tot[a] for a in areas}, "Blues", False),
              ("net (pp)", {a: 100 * (e[a] - i[a]) / tot[a] for a in areas}, "RdBu_r", True)],
             areas, "By Allen area  (sorted by % excited)")
fig.suptitle(f"FACE POKE  -  {len(poke):,} cells, {n_rec} recordings, {n_mice} mice "
             f"(13 + 18 Aug 2025)\n"
             "stimulus-locked contrast on the same cells within a session; "
             "anterior up, bregma marked",
             fontsize=13.5, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .90])
fig.savefig(os.path.join(OUT, "fig5_heatmap_facepoke.png"), dpi=160)
plt.close(fig)
print("fig5 face poke done")

# =============================================================== FIG 6 : melittin
only = [r for r in mel if r["cohort"] != "melittin_morphine"]
morp = [r for r in mel if r["cohort"] == "melittin_morphine"]
pe1, pi1, pn1 = maps(only)
pe2, pi2, pn2 = maps(morp)
# shared scale WITHIN this experiment, since the two arms are the same measurement
emax = np.nanpercentile(np.concatenate([pe1[np.isfinite(pe1)], pe2[np.isfinite(pe2)]]), 99)
imax = np.nanpercentile(np.concatenate([pi1[np.isfinite(pi1)], pi2[np.isfinite(pi2)]]), 99)

t1, e1, i1 = region_table(only)
t2, e2, i2 = region_table(morp)
areas_m = sorted([a for a in set(t1) | set(t2)
                  if t1.get(a, 0) >= MIN_CELLS and t2.get(a, 0) >= MIN_CELLS],
                 key=lambda a: -(e1[a] / t1[a]))

fig = plt.figure(figsize=(15, 13))
gs = fig.add_gridspec(3, 3, height_ratios=[2.2, 2.2, 1.35], hspace=.28, wspace=.18)
lim = max(np.nanpercentile(np.abs(pn1), 98), np.nanpercentile(np.abs(pn2), 98))
for row, (lab, pe_, pi_, pn_) in enumerate(
        [("Melittin only  (8 rec, 4 mice)", pe1, pi1, pn1),
         ("Melittin + morphine  (4 rec, 2 mice)", pe2, pi2, pn2)]):
    draw(fig.add_subplot(gs[row, 0]), pe_, "Reds", f"{lab}\nUp vs baseline",
         vmin=0, vmax=emax)
    draw(fig.add_subplot(gs[row, 1]), pi_, "Blues", f"{lab}\nDown vs baseline",
         vmin=0, vmax=imax)
    draw(fig.add_subplot(gs[row, 2]), pn_, "RdBu_r", f"{lab}\nNet (up - down)",
         norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim), cbl="percentage points")

sub = gs[2, :].subgridspec(3, 1, hspace=.12)
region_strip([fig.add_subplot(sub[k]) for k in range(3)],
             [("mel % up", {a: 100 * e1[a] / t1[a] for a in areas_m}, "Reds", False),
              ("mel+mor % up", {a: 100 * e2[a] / t2[a] for a in areas_m}, "Reds", False),
              ("delta (pp)", {a: 100 * (e2[a] / t2[a] - e1[a] / t1[a]) for a in areas_m},
               "RdBu_r", True)],
             areas_m, "By Allen area  (% up, sorted by melittin-only)")
fig.suptitle(f"MELITTIN  -  {len(mel):,} cells, 12 recordings, 6 mice\n"
             "injection session vs that animal's own baseline session; cells are NOT tracked "
             "across the injection (independent EXTRACT IDs)",
             fontsize=13.5, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .91])
fig.savefig(os.path.join(OUT, "fig6_heatmap_melittin.png"), dpi=160)
plt.close(fig)
print("fig6 melittin done")

print(f"\nface poke  : excited {100*sum(e.values())/sum(tot.values()):.1f}%  "
      f"inhibited {100*sum(i.values())/sum(tot.values()):.1f}%")
print(f"melittin   : up {100*sum(e1.values())/sum(t1.values()):.1f}%  "
      f"down {100*sum(i1.values())/sum(t1.values()):.1f}%")
print(f"mel+morph  : up {100*sum(e2.values())/sum(t2.values()):.1f}%  "
      f"down {100*sum(i2.values())/sum(t2.values()):.1f}%")
print("->", OUT)
