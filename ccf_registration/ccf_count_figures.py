"""Cell-count figures: absolute numbers, not percentages.

Bar length encodes the number of cells, so a long bar means the FOV sampled that region
heavily, not that the region responded more. The percentage is printed alongside for the
comparison that is actually fair between regions.

Outputs -> figures_v2/
  fig9_counts_facepoke.png
  fig10_counts_melittin.png
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

RED, BLUE, GREY = "#C0392B", "#2C6FBB", "#D8D8D8"


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


leaf = read(os.path.join(CNT, "counts_by_area_leaf.csv"))
par = read(os.path.join(CNT, "counts_by_area_parent.csv"))
pm = read(os.path.join(CNT, "counts_by_area_and_mouse.csv"))
for t in (leaf, par, pm):
    for r in t:
        for k in ("n_cells", "n_excited", "n_inhibited", "n_nonsig"):
            r[k] = int(r[k])


def stacked(ax, recs, title, min_n=200, label_key="area"):
    recs = [r for r in recs if r["n_cells"] >= min_n]
    recs.sort(key=lambda r: r["n_cells"])
    y = np.arange(len(recs))
    e = np.array([r["n_excited"] for r in recs])
    i = np.array([r["n_inhibited"] for r in recs])
    z = np.array([r["n_nonsig"] for r in recs])
    ax.barh(y, e, color=RED, label="excited", edgecolor="white", lw=.5)
    ax.barh(y, i, left=e, color=BLUE, label="inhibited", edgecolor="white", lw=.5)
    ax.barh(y, z, left=e + i, color=GREY, label="no change", edgecolor="white", lw=.5)
    tot = e + i + z
    for k, r in enumerate(recs):
        # counts inside their own segment when it is wide enough, else skip
        if e[k] > .045 * tot.max():
            ax.text(e[k] / 2, k, f"{e[k]:,}", ha="center", va="center", fontsize=7,
                    color="white", fontweight="bold")
        if i[k] > .045 * tot.max():
            ax.text(e[k] + i[k] / 2, k, f"{i[k]:,}", ha="center", va="center", fontsize=7,
                    color="white", fontweight="bold")
        ax.text(tot[k] + .012 * tot.max(), k,
                f"n={tot[k]:,}   {100*i[k]/tot[k]:.0f}% inh", va="center", fontsize=7.5,
                color="#333333")
    ax.set_yticks(y)
    ax.set_yticklabels([r[label_key] for r in recs], fontsize=9)
    ax.set_xlim(0, tot.max() * 1.30)
    ax.set_xlabel("number of cells")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)


def per_mouse(ax, recs, cond, groups, title):
    recs = [r for r in recs if r["condition"] == cond and r["area_group"] in groups]
    mice = sorted({r["mouse"] for r in recs})
    d = {(r["area_group"], r["mouse"]): r for r in recs}
    nG, nM = len(groups), len(mice)
    width = .8 / nM
    x = np.arange(nG)
    for j, m in enumerate(mice):
        off = (j - (nM - 1) / 2) * width
        exc = [d[(g, m)]["n_excited"] if (g, m) in d else 0 for g in groups]
        inh = [d[(g, m)]["n_inhibited"] if (g, m) in d else 0 for g in groups]
        ax.bar(x + off, exc, width * .95, color=RED, edgecolor="white", lw=.4)
        ax.bar(x + off, inh, width * .95, bottom=exc, color=BLUE, edgecolor="white", lw=.4)
        for k in range(nG):
            if exc[k] + inh[k] > 0:
                ax.text(x[k] + off, exc[k] + inh[k], m, rotation=90, ha="center",
                        va="bottom", fontsize=6, color="#555555")
    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=10)
    ax.set_ylabel("number of responsive cells")
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.margins(y=.16)


# ============================================================ FIG 9 : face poke
fp_leaf = [r for r in leaf if r["condition"] == "Face poke"]
fp_par = [r for r in par if r["condition"] == "Face poke"]
tot_n = sum(r["n_cells"] for r in fp_par)
tot_e = sum(r["n_excited"] for r in fp_par)
tot_i = sum(r["n_inhibited"] for r in fp_par)

fig = plt.figure(figsize=(16, 9.5))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15], height_ratios=[1.05, 1], hspace=.34,
                      wspace=.26)
stacked(fig.add_subplot(gs[:, 0]), fp_leaf, "By Allen area")
stacked(fig.add_subplot(gs[0, 1]), fp_par,
        "Rolled up  (SSp = all SSp-*, VIS = all visual, MO = MOp+MOs)",
        label_key="area_group")
per_mouse(fig.add_subplot(gs[1, 1]), pm, "Face poke", ["SSp", "VIS", "MO", "RSP"],
          "Responsive cells per mouse  (bar = excited + inhibited)")
fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for c in (RED, BLUE, GREY)],
           labels=["excited", "inhibited", "no change"], loc="lower center", ncol=3,
           frameon=False, fontsize=10, bbox_to_anchor=(.5, .005))
fig.suptitle(f"Fig 9  |  FACE POKE cell counts    "
             f"{tot_n:,} cells   {tot_e:,} excited ({100*tot_e/tot_n:.1f}%)   "
             f"{tot_i:,} inhibited ({100*tot_i/tot_n:.1f}%)   4 mice, 18 recordings",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, .035, 1, .945])
fig.savefig(os.path.join(OUT, "fig9_counts_facepoke.png"), dpi=160)
plt.close(fig)
print("fig9 done")

# ============================================================ FIG 10 : melittin
fig = plt.figure(figsize=(16, 10.5))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.15], hspace=.34, wspace=.26)
for row, cond in enumerate(["Melittin", "Melittin+morphine"]):
    lf = [r for r in leaf if r["condition"] == cond]
    pr = [r for r in par if r["condition"] == cond]
    n = sum(r["n_cells"] for r in pr)
    e = sum(r["n_excited"] for r in pr)
    i = sum(r["n_inhibited"] for r in pr)
    stacked(fig.add_subplot(gs[row, 0]), lf,
            f"{cond}  -  by Allen area\n{n:,} cells, {e:,} up, {i:,} down", min_n=400)
    stacked(fig.add_subplot(gs[row, 1]), pr, f"{cond}  -  rolled up",
            min_n=100, label_key="area_group")
fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for c in (RED, BLUE, GREY)],
           labels=["up vs baseline", "down vs baseline", "no change"], loc="lower center",
           ncol=3, frameon=False, fontsize=10, bbox_to_anchor=(.5, .005))
fig.suptitle("Fig 10  |  MELITTIN cell counts    "
             "melittin 80,281 cells / 4 mice     melittin + morphine 61,083 cells / 2 mice",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, .03, 1, .95])
fig.savefig(os.path.join(OUT, "fig10_counts_melittin.png"), dpi=160)
plt.close(fig)
print("fig10 done")
print("->", OUT)
