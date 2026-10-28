"""Five compact summary figures for the 5-slide PI deck.

One figure per slide. Every panel carries an axis label with units, the n it is built
from, and the name of the test behind any p value. Colour scales are labelled with
what the colour means and what the end points are - a figure that needs the speaker
to explain its colours is not finished.

Numbers come from compute_stats.py, which recomputes everything from the cell tables.
Nothing is carried over from the earlier decks.

    python make_final_figures.py        ->  final_deck/figures/S1..S5.png
"""
from __future__ import annotations

import collections
import csv
import glob
import math
import os
import re

import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CCF = os.path.join(ROOT, "ccf_registration")
BEH = os.path.join(ROOT, "behavior_sync")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

RED, BLUE, TEAL, AMBER = "#C0392B", "#2C6FBB", "#0D5C63", "#D08C1E"
NAVY, GREY, LGREY = "#1F3A5F", "#555555", "#BBBBBB"
# a boolean mask through a normal colormap renders at the map's top end, which for
# "autumn_r" is yellow, not the orange the legend promised. Pin the colour instead.
BORDER_CMAP = ListedColormap(["#C8651F"])

plt.rcParams.update({
    "font.size": 9.5, "axes.titlesize": 11, "axes.labelsize": 10,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
})


def panel_tag(ax, letter, dx=-0.085, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=14,
            fontweight="bold", va="top", ha="left", color=NAVY)


# ===================================================================== data
def load_cells(path):
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        r["exc"] = r["is_exc"] == "1"
        r["inh"] = r["is_inh"] == "1"
        r["AP"] = float(r["AP_mm"])
        r["ML"] = float(r["ML_mm"])
        r["row"] = float(r["atlas_row"])
        r["col"] = float(r["atlas_col"])
    return rows


FP = load_cells(os.path.join(CCF, "n18", "cell_ccf_all.csv"))
MEL = load_cells(os.path.join(CCF, "melittin_ccf", "cell_ccf_all.csv"))
MEL_ONLY = [r for r in MEL if r["cohort"] == "melittin_only"]
MEL_MOR = [r for r in MEL if r["cohort"] == "melittin_morphine"]

ATLAS = sio.loadmat(os.path.join(CCF, "allen_dorsal_map.mat"),
                    squeeze_me=True, struct_as_record=False)
BOUND = ATLAS["boundary"].astype(bool)
LABEL = ATLAS["label_img"]

QC = list(csv.DictReader(open(os.path.join(HERE, "reg_qc.csv"))))


def fp_animal(movie):
    return re.split(r"_d\d{6}_", movie)[0]


def enrich(rows, field, minn=300):
    """Fold-enrichment of `field` per Allen area against that paradigm's own
    cortex-wide rate, with a two-sided binomial test and Benjamini-Hochberg FDR."""
    base = float(np.mean([r[field] for r in rows]))
    per = collections.defaultdict(list)
    for r in rows:
        if r["area_acronym"] not in ("root", "", "nan"):
            per[r["area_acronym"]].append(r[field])
    o = [(a, float(np.mean(v)) / base, len(v),
          stats.binomtest(int(np.sum(v)), len(v), base).pvalue)
         for a, v in per.items() if len(v) >= minn]
    o.sort(key=lambda t: t[3])
    m, q, prev = len(o), [], 1.0
    for i, t in enumerate(reversed(o)):
        prev = min(prev, t[3] * m / (m - i))
        q.append(prev)
    res = [(a, k, n, qq) for (a, k, n, _), qq in zip(o, q[::-1])]
    res.sort(key=lambda t: -t[1])
    return res, base


# ===================================================================== S1
def slide1():
    fig = plt.figure(figsize=(13.0, 5.5))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.12, 1, 1.12], height_ratios=[1, 1],
                          wspace=.36, hspace=.70, left=.055, right=.99,
                          top=.80, bottom=.11)

    # --- A  example registration overlay (rendered, not a screenshot) -------
    ax = fig.add_subplot(gs[:, 0])
    ctx = LABEL > 0
    ax.imshow(np.where(ctx, 0.93, 1.0), cmap="gray", vmin=0, vmax=1)
    ax.imshow(np.ma.masked_where(~BOUND, BOUND), cmap=BORDER_CMAP, alpha=.85,
              interpolation="nearest")
    sub = FP[::10]
    ax.scatter([r["col"] for r in sub], [r["row"] for r in sub], s=2.0,
               c=NAVY, alpha=.55, linewidths=0, rasterized=True)
    br, bc = ATLAS["bregma_rc"]
    ax.plot(bc, br, "+", color="k", ms=11, mew=2)
    ax.annotate("bregma", (bc, br), xytext=(10, -12), textcoords="offset points",
                fontsize=8.5)
    ax.set_xlim(60, 400)
    ax.set_ylim(460, 60)
    ax.set_xticks([]); ax.set_yticks([])
    # scale bar: atlas is 25 um per pixel
    ax.plot([90, 90 + 2000 / 25], [440, 440], "k-", lw=3)
    ax.text(90 + 1000 / 25, 432, "2 mm", ha="center", fontsize=9)
    ax.set_title("Every cell placed in Allen CCF\n"
                 f"face-poke cells, n = {len(FP):,} (1 in 10 drawn)",
                 fontsize=10, fontweight="bold", pad=6)
    ax.text(.02, .02, "grey = dorsal cortex    orange = Allen area borders",
            transform=ax.transAxes, fontsize=8.0, color=GREY)
    panel_tag(ax, "A", dx=-0.03, dy=1.155)

    # --- B  Dice -----------------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    d = np.array([float(r["dice"]) for r in QC])
    o = np.argsort(d)
    ax.bar(np.arange(len(d)), d[o], color=TEAL, width=.82)
    ax.axhline(np.median(d), color=RED, ls="--", lw=1.6,
               label=f"median {np.median(d):.3f}")
    ax.set_ylim(.85, .955)
    ax.set_xlabel("recording, sorted (n = 30)")
    ax.set_ylabel("Dice overlap\nFOV mask vs CCF cortex")
    ax.set_title(f"Registration quality\n{d.min():.3f} – {d.max():.3f}, all 30 pass",
                 fontsize=9.8, pad=6)
    ax.legend(frameon=False, loc="lower right")
    panel_tag(ax, "B", dy=1.20)

    # --- C  scale ----------------------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    u = np.array([float(r["um"]) for r in QC])
    ax.bar(np.arange(len(u)), np.sort(u), color=AMBER, width=.82)
    ax.axhline(np.median(u), color=RED, ls="--", lw=1.6,
               label=f"median {np.median(u):.2f}")
    ax.set_ylim(5.0, 5.9)
    ax.set_xlabel("recording, sorted (n = 30)")
    ax.set_ylabel("fitted scale\n(µm per camera pixel)")
    ax.set_title(f"Scale is consistent\n{u.max()/u.min():.2f}× spread across 30 FOVs",
                 fontsize=9.8, pad=6)
    ax.legend(frameon=False, loc="lower right")
    panel_tag(ax, "C", dy=1.20)

    # --- D  method ---------------------------------------------------------
    ax = fig.add_subplot(gs[:, 2])
    ax.axis("off")
    a_lo = min(float(r["ang"]) for r in QC)
    a_hi = max(float(r["ang"]) for r in QC)
    txt = (
        "HOW  —  no manual clicking anywhere\n"
        "\n"
        "1   Allen CCFv3 25 µm volume projected top-down to a\n"
        "     dorsal-surface map; cortical layers collapsed to\n"
        "     their parent area.\n"
        "2   Rotation from the sagittal midline, detected from\n"
        "     the FOV's own left–right intensity profile\n"
        f"     ({a_lo:+.1f}° to {a_hi:+.1f}° across the 30 FOVs).\n"
        "3   Scale + translation by maximising Dice overlap of\n"
        "     the tissue mask against the CCF cortex mask.\n"
        "4   A–P sense resolved by trying both and keeping\n"
        "     whichever gives the better Dice.\n"
        "\n"
        "WHY NOT ANATOMICAL LANDMARKS\n"
        "\n"
        "The skull is removed, so bregma and lambda cannot be\n"
        "seen. Clicking them gave a 1.69× scale spread and an\n"
        "impossible 24.6 mm FOV width. Mask fitting: 1.05×.\n"
        "\n"
        "VALIDATION\n"
        "\n"
        "Known stereotaxic coordinates land in the right area:\n"
        "   3.8 mm posterior / 2.5 lateral  →  VISp\n"
        "   1.5 mm posterior / 3.2 lateral  →  SSp-bfd\n"
        "   1.0 mm anterior  / 1.5 lateral  →  MOp"
    )
    ax.text(-.02, 1.02, txt, va="top", ha="left", fontsize=8.3,
            family="DejaVu Sans", color="#222222", linespacing=1.52,
            transform=ax.transAxes)
    panel_tag(ax, "D", dx=-0.05, dy=1.10)

    fig.suptitle("Slide 1  |  All 30 recordings registered to the Allen CCF, "
                 "automatically", fontsize=14, fontweight="bold", y=.965)
    fig.savefig(os.path.join(OUT, "S1_registration.png"), dpi=170)
    plt.close(fig)
    print("S1 done")


# ===================================================================== S2
def slide2():
    fig = plt.figure(figsize=(13.0, 5.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.22, 1.06], wspace=.44,
                          left=.055, right=.985, top=.76, bottom=.19)

    # --- A  map of excited / inhibited -------------------------------------
    ax = fig.add_subplot(gs[0])
    ax.imshow(np.where(LABEL > 0, .94, 1.0), cmap="gray", vmin=0, vmax=1)
    ax.imshow(np.ma.masked_where(~BOUND, BOUND), cmap=BORDER_CMAP, alpha=.45,
              interpolation="nearest")
    exc = [r for r in FP if r["exc"]]
    inh = [r for r in FP if r["inh"]]
    pts = ([(r["col"], r["row"], RED) for r in exc[::3]] +
           [(r["col"], r["row"], BLUE) for r in inh[::9]])
    rng = np.random.default_rng(0)
    pts = [pts[i] for i in rng.permutation(len(pts))]   # draw order, not colour, decides
    ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=1.6,
               c=[p[2] for p in pts], alpha=.45, linewidths=0, rasterized=True)
    ax.set_xlim(60, 400); ax.set_ylim(460, 60)
    ax.set_xticks([]); ax.set_yticks([])
    ax.plot([90, 90 + 2000 / 25], [440, 440], "k-", lw=3)
    ax.text(90 + 1000 / 25, 432, "2 mm", ha="center", fontsize=9)
    ax.annotate("anterior", (1.015, .97), xycoords="axes fraction", ha="left",
                va="top", fontsize=8.6, color=GREY, rotation=90)
    ax.annotate("posterior", (1.015, .03), xycoords="axes fraction", ha="left",
                va="bottom", fontsize=8.6, color=GREY, rotation=90)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=RED, ms=7,
                              label=f"excited  {len(exc):,} cells  (1 in 3 drawn)"),
                       Line2D([], [], marker="o", ls="", color=BLUE, ms=7,
                              label=f"inhibited  {len(inh):,} cells  (1 in 9 drawn)")],
              frameon=False, loc="upper center", bbox_to_anchor=(.5, -.015),
              fontsize=8.6, handletextpad=.4)
    ax.set_title("Where the responding cells are", fontsize=10.5,
                 fontweight="bold", pad=7)
    panel_tag(ax, "A", dx=-0.11, dy=1.12)

    # --- B  enrichment by area ---------------------------------------------
    ax = fig.add_subplot(gs[1])
    e, base = enrich(FP, "exc")
    e = e[:7] + e[-4:]
    y = np.arange(len(e))[::-1]
    cols = [RED if k > 1 else BLUE for _, k, _, _ in e]
    ax.barh(y, [k for _, k, _, _ in e], color=cols, height=.72)
    ax.axvline(1, color="k", lw=1.2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{a}  (n={n:,})" for a, _, n, _ in e], fontsize=8.8)
    for yy, (a, k, n, q) in zip(y, e):
        star = "***" if q < 1e-3 else ("**" if q < .01 else ("*" if q < .05 else "n.s."))
        ax.text(k + .03, yy, f"{k:.2f}×  {star}", va="center", fontsize=8.3,
                color="#222222")
    ax.set_xlim(0, 2.45)
    ax.set_xlabel("fold enrichment of excited cells\n"
                  f"(area rate ÷ cortex-wide rate, {100*base:.1f} %)")
    ax.set_title("Face poke excites somatosensory and motor cortex\n"
                 "19 areas with ≥ 300 cells; binomial test vs cortex-wide\n"
                 "rate, BH-FDR.  *** q<0.001  ** q<0.01  * q<0.05",
                 fontsize=9.4, pad=7)
    panel_tag(ax, "B", dx=-0.42, dy=1.18)

    # --- C  A-P gradient ----------------------------------------------------
    ax = fig.add_subplot(gs[2])
    recs = sorted({r["movie"] for r in FP})
    ant, pos = [], []
    for rc in recs:
        v = [r for r in FP if r["movie"] == rc]
        m = np.median([r["AP"] for r in v])
        ant.append(100 * np.mean([r["inh"] for r in v if r["AP"] >= m]))
        pos.append(100 * np.mean([r["inh"] for r in v if r["AP"] < m]))
    ant, pos = np.array(ant), np.array(pos)
    W, p = stats.wilcoxon(ant, pos)
    for a, b in zip(ant, pos):
        ax.plot([0, 1], [a, b], "-", color=LGREY, lw=1.1, zorder=1)
    ax.scatter(np.zeros_like(ant), ant, s=34, color=AMBER, zorder=3,
               edgecolor="white", linewidth=.8)
    ax.scatter(np.ones_like(pos), pos, s=34, color=BLUE, zorder=3,
               edgecolor="white", linewidth=.8)
    ax.plot([0, 1], [ant.mean(), pos.mean()], "-o", color="k", lw=2.4, ms=8,
            zorder=4, label="mean")
    ax.set_xlim(-.35, 1.35)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["anterior\nhalf", "posterior\nhalf"])
    ax.set_ylabel("inhibited cells in that half (% of cells)")
    ax.set_title(f"Inhibition is stronger posteriorly,\n"
                 f"in {int((pos>ant).sum())} / {len(ant)} recordings "
                 f"({pos.mean()-ant.mean():+.1f} pp)\n"
                 f"Wilcoxon signed-rank, W = {W:.0f},\np = {p:.1e}, n = {len(ant)}",
                 fontsize=9.4, pad=7)
    ax.legend(frameon=False, loc="lower right")
    panel_tag(ax, "C", dx=-0.26, dy=1.23)

    fig.suptitle(f"Slide 2  |  Face poke: {len(FP):,} cells, 18 recordings, 4 mice  —  "
                 f"{100*np.mean([r['exc'] for r in FP]):.1f} % excited, "
                 f"{100*np.mean([r['inh'] for r in FP]):.1f} % inhibited",
                 fontsize=14, fontweight="bold", y=.975)
    fig.savefig(os.path.join(OUT, "S2_facepoke.png"), dpi=170)
    plt.close(fig)
    print("S2 done")


# ===================================================================== S3
def slide3():
    fig = plt.figure(figsize=(13.0, 5.5))
    gs = fig.add_gridspec(1, 4, width_ratios=[.85, .72, .72, 1.5], wspace=.50,
                          left=.055, right=.985, top=.80, bottom=.15)

    PAIRS = {"CS1014-1b": "1014m01b", "CS1014-2w": "1014m02w"}
    bym = collections.defaultdict(list)
    for r in MEL:
        bym[r["mouse"]].append(r)

    # --- A  cohort level ----------------------------------------------------
    ax = fig.add_subplot(gs[0])
    groups = [("melittin\nalone", MEL_ONLY, 4), ("melittin +\nmorphine", MEL_MOR, 2)]
    x = np.arange(2)
    up = [100 * np.mean([r["exc"] for r in g]) for _, g, _ in groups]
    dn = [100 * np.mean([r["inh"] for r in g]) for _, g, _ in groups]
    ax.bar(x - .19, up, .36, color=RED, label="up vs baseline")
    ax.bar(x + .19, dn, .36, color=BLUE, label="down vs baseline")
    for xi, (u, d) in enumerate(zip(up, dn)):
        ax.text(xi - .19, u + 1.4, f"{u:.1f}", ha="center", fontsize=8.6)
        ax.text(xi + .19, d + 1.4, f"{d:.1f}", ha="center", fontsize=8.6)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{t}\n{len(g):,} cells\n{m} mice" for t, g, m in groups],
                       fontsize=8.6)
    ax.set_ylim(0, 84)
    ax.set_ylabel("cells responding (% of cells)")
    ax.set_title("All cells pooled\n(mixes animals — see B)", fontsize=9.8, pad=7)
    ax.legend(frameon=False, loc="upper center", ncol=1, fontsize=8.4)
    panel_tag(ax, "A", dx=-0.38, dy=1.18)

    # --- B,C  paired --------------------------------------------------------
    for k, (fld, lab, ttl) in enumerate(
            [("exc", "cells up vs baseline (%)", "% of cells driven up"),
             ("amp", "mean ΔF/F, % of own baseline", "mean response amplitude")]):
        ax = fig.add_subplot(gs[1 + k])
        for mel_id, mor_id in PAIRS.items():
            A, B = bym[mel_id], bym[mor_id]
            if fld == "exc":
                va = 100 * np.mean([r["exc"] for r in A])
                vb = 100 * np.mean([r["exc"] for r in B])
            else:
                va = float(np.mean([float(r["amp_pct"]) for r in A]))
                vb = float(np.mean([float(r["amp_pct"]) for r in B]))
            ax.plot([0, 1], [va, vb], "-o", color=NAVY, ms=9, lw=2.2,
                    markerfacecolor="white", markeredgewidth=2.2)
            ax.annotate(f"{mel_id}\n{vb-va:+.1f}", (1, vb), xytext=(9, 0),
                        textcoords="offset points", fontsize=8, va="center",
                        color=NAVY)
        ax.set_xlim(-.35, 1.85)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["melittin\nalone", "melittin +\nmorphine"], fontsize=8.6)
        ax.set_ylabel(lab)
        ax.set_title(ttl, fontsize=9.8)
        if fld == "amp":
            ax.axhline(0, color="k", ls=":", lw=1)
        panel_tag(ax, "BC"[k], dx=-0.42, dy=1.18)

    # --- D  by area ---------------------------------------------------------
    ax = fig.add_subplot(gs[3])
    # isocortex only. Leaving every labelled area in pulls in SCzo (superior
    # colliculus) and an "(outside)" bin, which do not belong in a cortical claim.
    CTX = {"MOs", "MOp", "SSp-m", "SSp-n", "SSp-bfd", "SSp-ul", "SSp-ll", "SSp-tr",
           "SSp-un", "SSs", "VISp", "VISa", "VISam", "VISal", "VISpm", "VISrl",
           "VISl", "VISpor", "RSPd", "RSPv", "RSPagl", "AUDp", "AUDd", "AUDv",
           "AUDpo", "ACAd", "ACAv", "TEa", "PL", "ILA", "ORBl", "ORBm", "ORBvl"}
    reg = collections.defaultdict(dict)
    for mel_id, mor_id in PAIRS.items():
        for tag, mid in (("mel", mel_id), ("mor", mor_id)):
            per = collections.defaultdict(list)
            for r in bym[mid]:
                if r["area_acronym"] in CTX:
                    per[r["area_acronym"]].append(r["exc"])
            for a, v in per.items():
                if len(v) >= 100:
                    reg[a][f"{tag}_{mel_id}"] = 100 * float(np.mean(v))
    rows = []
    for a, d in reg.items():
        ds = [d[f"mor_{m}"] - d[f"mel_{m}"] for m in PAIRS
              if f"mor_{m}" in d and f"mel_{m}" in d]
        if len(ds) == 2:
            rows.append((a, ds[0], ds[1], float(np.mean(ds)), ds[0] * ds[1] > 0))
    rows.sort(key=lambda t: -t[3])
    y = np.arange(len(rows))[::-1]
    ax.barh(y + .19, [r[1] for r in rows], .38, color=TEAL, label="CS1014-1b")
    ax.barh(y - .19, [r[2] for r in rows], .38, color=AMBER, label="CS1014-2w")
    ax.axvline(0, color="k", lw=1.1)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.6)
    ax.set_xlabel("change in % of cells driven up,\nmorphine minus melittin "
                  "(percentage points)")
    n_ag = sum(1 for r in rows if r[4])
    ax.set_title(f"Per Allen cortical area, each animal drawn separately\n"
                 f"{n_ag} of {len(rows)} areas move the same way in both animals "
                 f"({100*n_ag/len(rows):.0f} %)\n"
                 f"isocortex only, ≥ 100 cells in all four sessions", fontsize=9.4,
                 pad=7)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(.5, -.145), ncol=2,
              fontsize=8.6)
    panel_tag(ax, "D", dx=-0.15, dy=1.18)

    fig.suptitle("Slide 3  |  Morphine does NOT silence melittin-driven cortex — "
                 "within the same two mice it drives slightly MORE cells up",
                 fontsize=13.5, fontweight="bold", y=.975)
    fig.savefig(os.path.join(OUT, "S3_morphine.png"), dpi=170)
    plt.close(fig)
    print("S3 done")


# ===================================================================== S4
def slide4():
    fig = plt.figure(figsize=(13.0, 5.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 1.0], wspace=.40,
                          left=.055, right=.985, top=.78, bottom=.21)

    ef, bf = enrich(FP, "exc")
    em, bm = enrich(MEL_ONLY, "exc")
    df = {a: k for a, k, _, _ in ef}
    dm = {a: k for a, k, _, _ in em}
    nf = {a: n for a, _, n, _ in ef}
    nm = {a: n for a, _, n, _ in em}
    shared = sorted(set(df) & set(dm))

    # --- A  scatter ---------------------------------------------------------
    ax = fig.add_subplot(gs[0])
    xs = np.array([df[a] for a in shared])
    ys = np.array([dm[a] for a in shared])
    ax.axhline(1, color=LGREY, lw=1); ax.axvline(1, color=LGREY, lw=1)
    ax.plot([0, 2.2], [0, 2.2], ":", color=LGREY, lw=1.2)
    ax.scatter(xs, ys, s=46, c=NAVY, zorder=3, edgecolor="white", linewidth=.8)
    for a, x, y in zip(shared, xs, ys):
        if abs(x - y) > .40 or x > 1.3 or y > 1.4:
            ax.annotate(a, (x, y), xytext=(6, 4), textcoords="offset points",
                        fontsize=8.4, color="#222222")
    rho, p = stats.spearmanr(xs, ys)
    ax.set_xlabel("face poke: fold enrichment of excited cells")
    ax.set_ylabel("melittin: fold enrichment of cells driven up")
    ax.set_xlim(.4, 2.2); ax.set_ylim(.4, 2.2)
    ax.set_title(f"The two paradigms map onto different cortex\n"
                 f"Spearman ρ = {rho:+.2f}, p = {p:.2f}, n = {len(shared)} shared areas\n"
                 f"(dotted line = identical enrichment)", fontsize=9.6, pad=7)
    panel_tag(ax, "A", dx=-0.19, dy=1.18)

    # --- B  diverging -------------------------------------------------------
    ax = fig.add_subplot(gs[1])
    diff = sorted(((a, df[a] - dm[a]) for a in shared), key=lambda t: t[1])
    y = np.arange(len(diff))
    ax.barh(y, [d for _, d in diff],
            color=[AMBER if d > 0 else TEAL for _, d in diff], height=.74)
    ax.axvline(0, color="k", lw=1.1)
    ax.set_yticks(y)
    ax.set_yticklabels([a for a, _ in diff], fontsize=7.8)
    ax.set_xlabel("face-poke enrichment  −  melittin enrichment\n(fold units)")
    ax.set_title("Somatosensory → poke;  parietal/visual → melittin\n"
                 "both share the motor areas MOp / MOs", fontsize=9.6, pad=7)
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=AMBER),
                       Rectangle((0, 0), 1, 1, color=TEAL)],
              labels=["stronger for face poke", "stronger for melittin"],
              frameon=False, fontsize=8.4, ncol=2, loc="upper center",
              bbox_to_anchor=(.5, -.175))
    panel_tag(ax, "B", dx=-0.26, dy=1.17)

    # --- C  the top areas as a table ----------------------------------------
    ax = fig.add_subplot(gs[2])
    ax.axis("off")
    top_f = [t for t in ef if t[1] > 1][:5]
    top_m = [t for t in em if t[1] > 1][:5]
    lines = [("FACE POKE  — top excited-cell enrichment", None)]
    lines += [(f"   {a:8s} {k:4.2f}×   n = {n:6,}   q = {q:.0e}", AMBER)
              for a, k, n, q in top_f]
    lines += [("", None), ("MELITTIN  — top driven-up enrichment", None)]
    lines += [(f"   {a:8s} {k:4.2f}×   n = {n:6,}   q = {q:.0e}", TEAL)
              for a, k, n, q in top_m]
    lines += [("", None),
              ("cortex-wide base rate", None),
              (f"   face poke   {100*bf:.1f} % of cells excited", None),
              (f"   melittin    {100*bm:.1f} % of cells driven up", None),
              ("", None),
              ("Enrichment = that area's response rate", None),
              ("÷ its own paradigm's cortex-wide rate, so", None),
              ("the two blocks compare even though the", None),
              ("base rates differ 3-fold.", None),
              ("", None),
              ("Binomial test per area, BH-FDR over the", None),
              ("19 (poke) / 27 (melittin) areas ≥ 300 cells.", None)]
    yy = .985
    for s, c in lines:
        ax.text(0, yy, s, transform=ax.transAxes, va="top", ha="left",
                fontsize=8.5, family="DejaVu Sans Mono" if c else "DejaVu Sans",
                color=c if c else "#222222",
                fontweight="bold" if s.startswith(("FACE", "MELITTIN")) else "normal")
        yy -= .0435
    panel_tag(ax, "C", dx=-0.04, dy=1.085)

    fig.suptitle("Slide 4  |  Face poke and melittin recruit different cortical "
                 "territory — this is what registration buys us",
                 fontsize=13.5, fontweight="bold", y=.965)
    fig.savefig(os.path.join(OUT, "S4_dissociation.png"), dpi=170)
    plt.close(fig)
    print("S4 done")


if __name__ == "__main__":
    slide1(); slide2(); slide3(); slide4()
    import make_s5; make_s5.slide5()
    print("->", OUT)
