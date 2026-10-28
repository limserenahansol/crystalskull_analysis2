"""Melittin MEAN AMPLITUDE map in CCF space - the melittin counterpart of fig7.

Quantity plotted is amp_pct = 100 * (injection mean dF/F - baseline session mean) / baseline
mean, i.e. percent change from that animal's own baseline. The raw amp cannot be pooled:
expression differs enough between recordings that raw per-recording means span -7.3 to +3.5
in arbitrary dF/F units, while the normalised version spans -26.7% to +4.8%.

Like fig7 this plots a continuous signed value rather than the fraction of cells clearing a
significance threshold, so it is directly comparable in *kind* to the face-poke amplitude
map - though not in magnitude, since melittin is a session-level contrast and face poke is
a stimulus-locked one.

Output -> figures_v2/fig8_melittin_amplitude_ccf.png
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
from scipy.stats import wilcoxon, mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures_v2")

A = sio.loadmat(os.path.join(HERE, "allen_dorsal_map.mat"))
LABEL = A["label_img"]
BOUND = A["boundary"].astype(bool)
BREGMA = A["bregma_rc"].ravel() + 1
ACR = [str(x[0]) for x in A["acronyms"].ravel()]
IDS = A["ids"].ravel()
ID2ACR = {int(i): a for i, a in zip(IDS, ACR)}
CTX = {"MOs", "MOp", "SSp-m", "SSp-n", "SSp-bfd", "SSp-ul", "SSp-ll", "SSp-tr", "SSp-un",
       "SSs", "VISp", "VISa", "VISam", "VISal", "VISpm", "VISrl", "VISl", "VISpor",
       "RSPd", "RSPv", "RSPagl", "AUDp", "AUDd", "AUDv", "ACAd", "TEa"}
CTX_MASK = np.isin(LABEL, [i for i, a in ID2ACR.items() if a in CTX])
RED, BLUE = "#C0392B", "#2C6FBB"

rows = list(csv.DictReader(open(os.path.join(HERE, "melittin_ccf", "cell_ccf_all.csv"))))
for r in rows:
    r["AP"] = float(r["AP_mm"])
    r["amp"] = float(r["amp_pct"])
    r["rr"] = int(round(float(r["atlas_row"])))
    r["cc"] = int(round(float(r["atlas_col"])))

only = [r for r in rows if r["cohort"] != "melittin_morphine"]
morp = [r for r in rows if r["cohort"] == "melittin_morphine"]
h, w = LABEL.shape


def amp_map(rs, sig=6):
    num = np.zeros((h, w)); den = np.zeros((h, w))
    for r in rs:
        if 0 <= r["rr"] < h and 0 <= r["cc"] < w:
            num[r["rr"], r["cc"]] += r["amp"]
            den[r["rr"], r["cc"]] += 1
    ns, ds = gaussian_filter(num, sig), gaussian_filter(den, sig)
    ok = CTX_MASK & (ds > np.percentile(ds[ds > 0], 40))
    m = np.full((h, w), np.nan)
    m[ok] = ns[ok] / ds[ok]
    return m


def ap_split(rs):
    recs = collections.defaultdict(list)
    for r in rs:
        recs[r["label"]].append(r)
    p, a = [], []
    for lab, v in sorted(recs.items()):
        ap = np.array([x["AP"] for x in v]); am = np.array([x["amp"] for x in v])
        lo, hi = np.percentile(ap, [33.3, 66.7])
        p.append(am[ap <= lo].mean()); a.append(am[ap >= hi].mean())
    return np.array(p), np.array(a)


def by_area(rs):
    d = collections.defaultdict(list)
    for r in rs:
        if r["area_acronym"] in CTX:
            d[r["area_acronym"]].append((r["AP"], r["amp"]))
    return {a: (np.mean([x[0] for x in v]), np.mean([x[1] for x in v]), len(v))
            for a, v in d.items() if len(v) >= 300}


m1, m2 = amp_map(only), amp_map(morp)
lim = np.nanpercentile(np.abs(np.concatenate([m1[np.isfinite(m1)], m2[np.isfinite(m2)]])), 97)
a1, a2 = by_area(only), by_area(morp)
areas = sorted(set(a1) & set(a2), key=lambda x: -a1[x][0])   # anterior first


def draw_map(ax, m, title, vlim, cmap="bwr", cbl="mean change (% of baseline)"):
    ax.imshow(np.where(CTX_MASK, 0.94, np.nan), cmap="gray", vmin=0, vmax=1)
    im = ax.imshow(m, cmap=cmap, norm=TwoSlopeNorm(vcenter=0, vmin=-vlim, vmax=vlim))
    ax.imshow(np.ma.masked_where(~BOUND, BOUND), cmap="gray_r", alpha=.30,
              interpolation="nearest")
    ax.plot(BREGMA[1], BREGMA[0], "k+", ms=10, mew=1.6)
    ax.text(.5, .965, "A", transform=ax.transAxes, ha="center", color=RED,
            fontsize=15, fontweight="bold")
    ax.text(.5, .01, "P", transform=ax.transAxes, ha="center", color=BLUE,
            fontsize=15, fontweight="bold")
    ax.set_xlim(60, 420); ax.set_ylim(400, 60)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=10.5, fontweight="bold")
    cb = plt.colorbar(im, ax=ax, fraction=.045, pad=.02)
    cb.set_label(cbl, fontsize=8)
    cb.ax.tick_params(labelsize=7)


fig = plt.figure(figsize=(17.5, 14))
gs = fig.add_gridspec(3, 3, width_ratios=[1.15, .78, 1.5], wspace=.42, hspace=.34)

for row, (lab, rs, m) in enumerate([("Melittin only  (8 rec, 4 mice)", only, m1),
                                    ("Melittin + morphine  (4 rec, 2 mice)", morp, m2)]):
    draw_map(fig.add_subplot(gs[row, 0]), m, f"{lab}\nmean change vs own baseline", lim)

    p, a = ap_split(rs)
    ax = fig.add_subplot(gs[row, 1])
    for p_, a_ in zip(p, a):
        ax.plot([0, 1], [p_, a_], "-", color="0.65", lw=.9, zorder=1)
    ax.scatter(np.zeros_like(p), p, s=44, c=BLUE, zorder=3, edgecolor="w", lw=.6)
    ax.scatter(np.ones_like(a), a, s=44, c=RED, zorder=3, edgecolor="w", lw=.6)
    ax.plot([0, 1], [p.mean(), a.mean()], "k-", lw=3, zorder=4)
    ax.axhline(0, color="k", ls=":", lw=1)
    ax.set_xlim(-.35, 1.35); ax.set_xticks([0, 1])
    ax.set_xticklabels(["Posterior\nthird", "Anterior\nthird"])
    ax.set_ylabel("mean change (% of baseline)")
    try:
        _, pv = wilcoxon(a, p)
        ptxt = f"signrank p = {pv:.3f}"
    except ValueError:
        ptxt = "n too small"
    ax.set_title(f"{int((a > p).sum())}/{len(a)} recordings A > P\n{ptxt}",
                 fontsize=10.5, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    d = a1 if row == 0 else a2
    ax = fig.add_subplot(gs[row, 2])
    y = np.arange(len(areas))
    vals = [d[x][1] for x in areas]
    ax.barh(y, vals, color=[RED if v > 0 else BLUE for v in vals], edgecolor="white")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{x}  ({a1[x][0]:+.1f})" for x in areas], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, color="k", lw=1)
    ax.set_xlabel("mean change (% of baseline)")
    ax.set_title("By Allen area, anterior -> posterior", fontsize=10.5, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

# ---- difference row
diff = m2 - m1
dlim = np.nanpercentile(np.abs(diff[np.isfinite(diff)]), 97)
draw_map(fig.add_subplot(gs[2, 0]), diff,
         "Difference\n(melittin + morphine)  minus  (melittin only)", dlim, cmap="PuOr_r",
         cbl="percentage points")

ax = fig.add_subplot(gs[2, 1])
ax.axis("off")
p1, aa1 = ap_split(only)
p2, aa2 = ap_split(morp)
txt = "\n".join([
    "Two things are happening at once.",
    "",
    "1. A cortex-wide DROP, which is what",
    "   separates melittin from face poke:",
    f"     melittin     {np.mean(np.concatenate([p1, aa1])):+6.1f} %",
    f"     + morphine   {np.mean(np.concatenate([p2, aa2])):+6.1f} %",
    "   Face poke sits near zero overall.",
    "",
    "2. An anterior > posterior gradient ON",
    "   TOP of that drop, same direction as",
    "   face poke:",
    f"     melittin    {int((aa1 > p1).sum())}/{len(aa1)} recordings",
    f"     + morphine  {int((aa2 > p2).sum())}/{len(aa2)} recordings",
    f"     post {p1.mean():+.1f} %   ant {aa1.mean():+.1f} %",
    "",
    "The gradient is NOT specific to face",
    "poke. The offset is specific to melittin.",
    "",
    "n = 4 vs 2 mice, no vehicle arm:",
    "the difference map is descriptive.",
])
ax.text(0, .98, txt, va="top", ha="left", fontsize=8.6, color="#555555",
        family="monospace", linespacing=1.45, transform=ax.transAxes)

ax = fig.add_subplot(gs[2, 2])
y = np.arange(len(areas))
dv = [a2[x][1] - a1[x][1] for x in areas]
ax.barh(y, dv, color=["#B5651D" if v > 0 else "#5F3DC4" for v in dv], edgecolor="white")
ax.set_yticks(y)
ax.set_yticklabels([f"{x}  ({a1[x][0]:+.1f})" for x in areas], fontsize=8)
ax.invert_yaxis()
ax.axvline(0, color="k", lw=1)
ax.set_xlabel("morphine minus melittin (percentage points)")
ax.set_title("Difference by area", fontsize=10.5, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)

fig.suptitle("Fig 8  |  Melittin: mean amplitude in Allen coordinates\n"
             "injection session vs that animal's own baseline; same quantity type as the "
             "face-poke amplitude map (fig7), different measurement",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .94])
fig.savefig(os.path.join(OUT, "fig8_melittin_amplitude_ccf.png"), dpi=150)
plt.close(fig)

print(f"melittin only     A vs P: {int((aa1 > p1).sum())}/{len(aa1)}, "
      f"post {p1.mean():+.1f}%  ant {aa1.mean():+.1f}%")
print(f"melittin+morphine A vs P: {int((aa2 > p2).sum())}/{len(aa2)}, "
      f"post {p2.mean():+.1f}%  ant {aa2.mean():+.1f}%")
print(f"overall mean change: melittin {np.mean(np.concatenate([p1, aa1])):+.1f}%, "
      f"+morphine {np.mean(np.concatenate([p2, aa2])):+.1f}%")
print("->", os.path.join(OUT, "fig8_melittin_amplitude_ccf.png"))
