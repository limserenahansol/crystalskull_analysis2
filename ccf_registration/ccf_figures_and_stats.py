"""Superimposed CCF figures and per-region statistics for all three conditions.

Inputs (already produced):
  n18/cell_ccf_all.csv           face poke, 18 recordings, 4 mice   (peri-poke vs pre-poke)
  melittin_ccf/cell_ccf_all.csv  melittin +/- morphine, 12 recordings, 6 mice
                                 (injection session vs baseline session)

Outputs -> figures_v2/

A note on what is comparable. Face poke is a stimulus-locked, within-session contrast on the
same cells. Melittin is a between-session contrast with independent EXTRACT IDs. The
percentages are therefore NOT interchangeable across paradigms - a sustained session-level
shift is far easier to detect than a phasic response, which is most of why the melittin
fractions are ~3x larger. Compare regional *patterns* within a paradigm, and compare
melittin against melittin+morphine, which do share a measurement.
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
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures_v2")
os.makedirs(OUT, exist_ok=True)

RED, BLUE, GREY = "#C0392B", "#2C6FBB", "#555555"
COND_COL = {"Face poke": "#7B4FA8", "Melittin": "#0D5C63", "Melittin+morphine": "#D08C1E"}
MIN_CELLS = 300          # per region, pooled, to be reported at all


# ---------------------------------------------------------------- atlas
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


def load(path, cond_from):
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            r["_cond"] = cond_from(r)
            rows.append(r)
    return rows


poke = load(os.path.join(HERE, "n18", "cell_ccf_all.csv"), lambda r: "Face poke")
mel = load(os.path.join(HERE, "melittin_ccf", "cell_ccf_all.csv"),
           lambda r: "Melittin+morphine" if r["cohort"] == "melittin_morphine" else "Melittin")

# The face-poke recording index puts the FOV label (m1f, m2c, m3_m1f) in its "mouse"
# column, not the animal. The animal is the leading token of the movie name, and m201 was
# imaged on both days. Without this the figures claim 3 mice when there are 4.
for r in poke:
    mv = r["movie"]
    animal = "m1" if mv.startswith("0004749770_m1_") else mv.split("_")[0]
    day = "0813" if "d250813" in mv else "0818"
    r["mouse"] = animal
    r["mouse_day"] = f"{animal}-{day}"
for r in mel:
    r["mouse_day"] = r["mouse"]

allrows = poke + mel
CONDS = ["Face poke", "Melittin", "Melittin+morphine"]

for c in CONDS:
    s = [r for r in allrows if r["_cond"] == c]
    print(f"{c:20s} {len(s):7d} cells, {len(set(r['mouse'] for r in s))} mice, "
          f"{len(set(r['label'] for r in s))} recordings")


# ============================================================ FIG 1: CCF response maps
def response_map(rows, key="is_exc", sigma_px=6):
    """Density of a cell class in atlas space, normalised by total cell density."""
    h, w = LABEL.shape
    tot = np.zeros((h, w))
    sel = np.zeros((h, w))
    for r in rows:
        y, x = int(round(float(r["atlas_row"]))), int(round(float(r["atlas_col"])))
        if not (0 <= y < h and 0 <= x < w):
            continue
        tot[y, x] += 1
        if r[key] in ("1", "true", "True"):
            sel[y, x] += 1
    from scipy.ndimage import gaussian_filter
    tot_s = gaussian_filter(tot, sigma_px)
    sel_s = gaussian_filter(sel, sigma_px)
    out = np.full((h, w), np.nan)
    ok = tot_s > np.percentile(tot_s[tot_s > 0], 25) if (tot_s > 0).any() else tot_s > 0
    out[ok] = 100 * sel_s[ok] / tot_s[ok]
    return out, tot_s


fig, axes = plt.subplots(2, 3, figsize=(15, 11))
for j, cond in enumerate(CONDS):
    sub = [r for r in allrows if r["_cond"] == cond]
    for i, (key, cmap, ttl) in enumerate([("is_exc", "Reds", "Excited / up"),
                                          ("is_inh", "Blues", "Inhibited / down")]):
        ax = axes[i, j]
        m, _ = response_map(sub, key)
        ax.imshow(np.where(CTX_MASK, 0.93, np.nan), cmap="gray", vmin=0, vmax=1)
        im = ax.imshow(np.where(CTX_MASK, m, np.nan), cmap=cmap)
        ax.imshow(np.ma.masked_where(~BOUND, BOUND), cmap="gray_r", alpha=.35,
                  interpolation="nearest")
        ax.plot(BREGMA[1], BREGMA[0], "k+", ms=9, mew=1.5)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlim(60, 420); ax.set_ylim(400, 60)
        n_mice = len(set(r["mouse"] for r in sub))
        ax.set_title(f"{cond}\n{ttl}  ({len(sub):,} cells, {n_mice} mice)",
                     fontsize=10, color=COND_COL[cond], fontweight="bold")
        cb = fig.colorbar(im, ax=ax, fraction=.04, pad=.02)
        cb.set_label("% of local cells", fontsize=8)
        cb.ax.tick_params(labelsize=7)
fig.suptitle("Fig 1  |  Where responsive cells sit in the Allen CCF, by condition\n"
             "anterior up, bregma marked; grey = dorsal cortex, thin lines = area borders",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .93])
fig.savefig(os.path.join(OUT, "fig1_ccf_maps_by_condition.png"), dpi=160)
plt.close(fig)
print("fig1 done")


# ============================================================ region statistics
def per_region(rows):
    """Return {area: (n_all, n_exc, n_inh)} and per-mouse breakdown."""
    tot = collections.Counter()
    exc = collections.Counter()
    inh = collections.Counter()
    bymouse = collections.defaultdict(lambda: collections.Counter())
    for r in rows:
        a = r["area_acronym"]
        if a not in CTX:
            continue
        tot[a] += 1
        bymouse[(r["mouse"], a)]["n"] += 1
        if r["is_exc"] == "1":
            exc[a] += 1
            bymouse[(r["mouse"], a)]["e"] += 1
        if r["is_inh"] == "1":
            inh[a] += 1
            bymouse[(r["mouse"], a)]["i"] += 1
    return tot, exc, inh, bymouse


stats = {c: per_region([r for r in allrows if r["_cond"] == c]) for c in CONDS}

areas = sorted({a for c in CONDS for a, n in stats[c][0].items() if n >= MIN_CELLS})
print(f"\n{len(areas)} regions with >= {MIN_CELLS} cells in at least one condition")

# enrichment = (% excited in area) / (% excited overall), per condition
rows_out = []
for a in areas:
    rec = {"area": a}
    for c in CONDS:
        tot, exc, inh, _ = stats[c]
        n = tot.get(a, 0)
        if n < MIN_CELLS:
            rec[f"{c}_n"] = n
            rec[f"{c}_exc"] = np.nan
            rec[f"{c}_inh"] = np.nan
            rec[f"{c}_enr"] = np.nan
            continue
        base_e = sum(exc.values()) / max(sum(tot.values()), 1)
        rec[f"{c}_n"] = n
        rec[f"{c}_exc"] = 100 * exc.get(a, 0) / n
        rec[f"{c}_inh"] = 100 * inh.get(a, 0) / n
        rec[f"{c}_enr"] = (exc.get(a, 0) / n) / base_e if base_e > 0 else np.nan
    rows_out.append(rec)

with open(os.path.join(OUT, "region_statistics.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
    w.writeheader()
    for r in rows_out:
        w.writerow(r)


# ============================================================ FIG 2: enrichment by region
# One common row order across all three panels, so a row means the same region in each and
# you can read across. Sorting each panel independently under sharey=True silently
# mismatched the labels and the bars.
order = [r["area"] for r in sorted(
    rows_out,
    key=lambda r: (r.get("Face poke_enr") if np.isfinite(r.get("Face poke_enr", np.nan))
                   else -np.inf))]
xmin = min(v for r in rows_out for v in
           [r.get(f"{c}_enr", np.nan) for c in CONDS] if np.isfinite(v))
xmax = max(v for r in rows_out for v in
           [r.get(f"{c}_enr", np.nan) for c in CONDS] if np.isfinite(v))

fig, axes = plt.subplots(1, 3, figsize=(16, 7), sharey=True, sharex=True)
by_area = {r["area"]: r for r in rows_out}
for j, cond in enumerate(CONDS):
    ax = axes[j]
    vals = [by_area[a].get(f"{cond}_enr", np.nan) for a in order]
    y = np.arange(len(order))
    cols = [RED if (np.isfinite(v) and v > 1) else BLUE for v in vals]
    plot_v = [(v - 1) if np.isfinite(v) else 0 for v in vals]
    ax.barh(y, plot_v, left=1, color=cols, edgecolor="white")
    for i, v in enumerate(vals):
        if not np.isfinite(v):
            ax.text(1, i, "  n/a", va="center", fontsize=7, color=GREY)
    ax.axvline(1, color="k", lw=1)
    ax.set_xlim(max(0, xmin - .1), xmax + .1)
    ax.set_xlabel("excitation enrichment\n(area rate / overall rate)")
    ax.set_title(cond, color=COND_COL[cond], fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_yticks(np.arange(len(order)))
axes[0].set_yticklabels(order, fontsize=9)
fig.suptitle("Fig 2  |  Which regions are over-represented among excited cells\n"
             "1.0 = the region responds at the overall rate; area size is divided out",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .90])
fig.savefig(os.path.join(OUT, "fig2_region_enrichment.png"), dpi=160)
plt.close(fig)
print("fig2 done")


# ============================================================ FIG 3: melittin vs +morphine
# proper unit of analysis is the mouse
def mouse_rates(cond, area, key="e"):
    _, _, _, bym = stats[cond]
    out = []
    for (m, a), c in bym.items():
        if a == area and c["n"] >= 50:
            out.append(100 * c[key] / c["n"])
    return out


comp = []
for a in areas:
    x = mouse_rates("Melittin", a)
    y = mouse_rates("Melittin+morphine", a)
    if len(x) >= 3 and len(y) >= 2:
        p = mannwhitneyu(x, y, alternative="two-sided", method="exact").pvalue
        comp.append((a, np.mean(x), np.mean(y), np.mean(y) - np.mean(x), p, len(x), len(y)))

comp.sort(key=lambda t: t[3], reverse=True)
with open(os.path.join(OUT, "melittin_vs_morphine_by_region.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["area", "melittin_exc_pct", "morphine_exc_pct", "delta_pct", "p_exact",
                "n_mice_mel", "n_mice_mor"])
    w.writerows(comp)

fig, ax = plt.subplots(figsize=(9, max(5, .38 * len(comp) + 2)))
y = np.arange(len(comp))
delta = [c[3] for c in comp]
cols = [RED if d > 0 else BLUE for d in delta]
ax.barh(y, delta, color=cols, edgecolor="white")
for i, c in enumerate(comp):
    if c[4] < 0.05:
        ax.text(c[3] + (0.4 if c[3] > 0 else -0.4), i, "*", va="center",
                ha="left" if c[3] > 0 else "right", fontsize=13, fontweight="bold")
ax.set_yticks(y)
ax.set_yticklabels([c[0] for c in comp], fontsize=9)
ax.axvline(0, color="k", lw=1)
ax.set_xlabel("change in % excited cells,  melittin + morphine  minus  melittin alone\n"
              "(per-mouse means; * = exact Mann-Whitney p < 0.05, n = 4 vs 2 mice)")
ax.set_title("Fig 3  |  Does morphine change the melittin response, and where?",
             fontsize=13, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig3_melittin_vs_morphine_by_region.png"), dpi=160)
plt.close(fig)
print("fig3 done")

print(f"\nregions where morphine most increased excitation:")
for c in comp[:6]:
    print(f"  {c[0]:9s} {c[1]:5.1f}% -> {c[2]:5.1f}%   delta {c[3]:+5.1f}   p={c[4]:.3f}")
print(f"\nsignificant at p<0.05: {sum(1 for c in comp if c[4] < 0.05)} of {len(comp)} regions "
      f"(no multiple-comparison correction; with n=4 vs 2 mice treat all of this as descriptive)")
print("\n->", OUT)
