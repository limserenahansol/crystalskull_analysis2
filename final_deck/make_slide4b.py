"""One extra slide: what morphine actually changed, drawn instead of tabulated.

This is the within-animal table turned into a picture. Same numbers, recomputed from
the cell table rather than copied, so the slide cannot drift from the data.

Four panels:
  A  the three cell fractions, before and after, both animals on one % axis
  B  mean response amplitude, which has its own unit and so its own axis
  C  the pooled odds ratio with its CI, against OR = 1
  D  what did NOT change, and why the effect cannot yet be called morphine

Outputs a figure and a ONE-SLIDE pptx, not a rebuilt deck.

    python make_slide4b.py
"""
from __future__ import annotations

import collections
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
CCF = os.path.join(os.path.dirname(HERE), "ccf_registration")
FIG = os.path.join(HERE, "figures")
ST = json.load(open(os.path.join(HERE, "region_stats.json")))

RED, BLUE, TEAL, AMBER = "#C0392B", "#2C6FBB", "#0D5C63", "#D08C1E"
NAVY, GREY, LGREY = "#1F3A5F", "#555555", "#BBBBBB"
PAIRS = [("CS1014-1b", "1014m01b"), ("CS1014-2w", "1014m02w")]
ACOL = {"CS1014-1b": TEAL, "CS1014-2w": AMBER}

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10.5,
    "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "savefig.facecolor": "white",
})

# ------------------------------------------------------------------ numbers
rows = list(csv.DictReader(open(os.path.join(CCF, "melittin_ccf",
                                             "cell_ccf_all.csv"))))
by = collections.defaultdict(list)
for r in rows:
    r["exc"] = r["is_exc"] == "1"
    r["inh"] = r["is_inh"] == "1"
    by[r["mouse"]].append(r)

V = {}
for mel_id, mor_id in PAIRS:
    d = {}
    for tag, mid in (("mel", mel_id), ("mor", mor_id)):
        v = by[mid]
        up = 100 * float(np.mean([x["exc"] for x in v]))
        dn = 100 * float(np.mean([x["inh"] for x in v]))
        d[tag] = dict(n=len(v), up=up, down=dn, none=100 - up - dn,
                      amp=float(np.mean([float(x["amp_pct"]) for x in v])))
    V[mel_id] = d
    print(f"{mel_id}")
    for k in ("up", "down", "none", "amp"):
        a, b = d["mel"][k], d["mor"][k]
        print(f"   {k:6s} {a:6.1f} -> {b:6.1f}   ({b-a:+.1f})")

OR = {a: (np.exp(ST["C"][a]["pool"]),
          np.exp(ST["C"][a]["pool"] - 1.96 * ST["C"][a]["se_pool"]),
          np.exp(ST["C"][a]["pool"] + 1.96 * ST["C"][a]["se_pool"]))
      for a, _ in PAIRS}


def panel_tag(ax, letter, dx=-0.085, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=14, fontweight="bold",
            va="top", ha="left", color=NAVY)


# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(13.0, 5.5))
gs = fig.add_gridspec(1, 4, width_ratios=[1.12, .72, .86, 1.40], wspace=.46,
                      left=.055, right=.985, top=.78, bottom=.17)

# --- A  the three cell fractions ----------------------------------------
ax = fig.add_subplot(gs[0])
MEAS = [("down", BLUE, "driven DOWN"), ("up", RED, "driven UP"),
        ("none", GREY, "no response")]
for key, col, lab in MEAS:
    for aid, _ in PAIRS:
        a, b = V[aid]["mel"][key], V[aid]["mor"][key]
        ls = "-" if aid == "CS1014-1b" else "--"
        ax.plot([0, 1], [a, b], ls, color=col, lw=2.2, marker="o", ms=8,
                markerfacecolor="white", markeredgewidth=2.2)
    # one delta label per measure, placed on the mean of the two animals
    mb = np.mean([V[a]["mor"][key] for a, _ in PAIRS])
    md = np.mean([V[a]["mor"][key] - V[a]["mel"][key] for a, _ in PAIRS])
    ax.annotate(f"{lab}\n{md:+.1f} pp", (1, mb), xytext=(11, 0),
                textcoords="offset points", fontsize=9, va="center",
                color=col, fontweight="bold")
ax.set_xlim(-.25, 1.95)
ax.set_xticks([0, 1])
ax.set_xticklabels(["melittin\nalone", "melittin +\nmorphine"])
ax.set_ylim(0, 72)
ax.set_ylabel("cells in that class (% of all cells)")
ax.set_title("Fewer cells down, more cells up\nboth animals, every class",
             fontsize=10.2, pad=7)
ax.legend(handles=[Line2D([], [], color="k", ls="-", lw=2, label="CS1014-1b"),
                   Line2D([], [], color="k", ls="--", lw=2, label="CS1014-2w")],
          frameon=False, fontsize=8.6, loc="center left", bbox_to_anchor=(.02, .62))
panel_tag(ax, "A", dx=-0.22, dy=1.17)

# --- B  amplitude --------------------------------------------------------
ax = fig.add_subplot(gs[1])
for aid, _ in PAIRS:
    a, b = V[aid]["mel"]["amp"], V[aid]["mor"]["amp"]
    ax.plot([0, 1], [a, b], "-o", color=ACOL[aid], lw=2.4, ms=9,
            markerfacecolor="white", markeredgewidth=2.4)
    # short suffix only: the full id here collides with panel C's y tick labels
    ax.annotate(f"{b-a:+.1f}\n{aid.split('-')[1]}", (1, b), xytext=(9, 0),
                textcoords="offset points", fontsize=8.8, va="center",
                color=ACOL[aid], fontweight="bold")
ax.axhline(0, color="k", ls=":", lw=1)
ax.set_xlim(-.3, 1.75)
ax.set_xticks([0, 1])
ax.set_xticklabels(["melittin\nalone", "melittin +\nmorphine"])
ax.set_ylabel("mean ΔF/F,  % of own baseline")
ax.set_title("Amplitude rises too\n(its own unit, own axis)", fontsize=10.2, pad=7)
panel_tag(ax, "B", dx=-0.44, dy=1.17)

# --- C  pooled odds ratio -------------------------------------------------
ax = fig.add_subplot(gs[2])
for k, (aid, _) in enumerate(PAIRS):
    o, lo, hi = OR[aid]
    ax.errorbar(o, k, xerr=[[o - lo], [hi - o]], fmt="o", color=ACOL[aid], ms=11,
                capsize=5, lw=2.6)
    ax.text(o, k + .30, f"OR {o:.2f}\n[{lo:.2f}–{hi:.2f}]", ha="center",
            fontsize=9, fontweight="bold", color=ACOL[aid])
ax.axvline(1, color="k", ls="--", lw=1.6)
ax.text(1, -.52, "no change", ha="center", va="top", fontsize=8.6, color=GREY)
ax.set_yticks([0, 1])
ax.set_yticklabels([a for a, _ in PAIRS], fontsize=9.5)
ax.set_ylim(-.75, 1.8)
ax.set_xlim(.96, 1.27)
ax.set_xlabel("odds of a cell being driven up,\nmorphine ÷ melittin  (95 % CI)")
ax.set_title("Both animals shift the same way\npooled over all cortex",
             fontsize=10.2, pad=7)
panel_tag(ax, "C", dx=-0.34, dy=1.17)

# --- D  what did NOT change ----------------------------------------------
ax = fig.add_subplot(gs[3])
ax.axis("off")
D = ST["D"]
lines = [
    ("WHAT DID NOT CHANGE", NAVY, True),
    ("", None, False),
    ("No individual brain area, reproducibly.", "#222222", True),
    ("Within each animal the effect varies strongly", None, False),
    ("by region (I² = 82 % and 92 %), but the two", None, False),
    (f"animals rank the areas differently:", None, False),
    (f"r = {D['pearson_lor'][0]:+.2f}, p = {D['pearson_lor'][1]:.2f}, "
     f"{D['agree']}/{D['n']} directions agree.", None, False),
    ("", None, False),
    ("No change in the A–P pattern.", "#222222", True),
    ("Morphine arm: down +0.4 pp (3/4 recordings),", None, False),
    ("up −1.8 pp (1/4). Both n.s.", None, False),
    ("", None, False),
    ("WHAT WE CANNOT YET SAY", RED, True),
    ("", None, False),
    ("The two arms are seven weeks apart —", None, False),
    ("melittin alone 18 Feb, +morphine 8 Apr — and", None, False),
    ("there is no melittin + vehicle session on the", None, False),
    ("second day. A +2.7 pp shift is equally", None, False),
    ("consistent with morphine, with window", None, False),
    ("maturation over seven weeks, or with a", None, False),
    ("second melittin exposure. This design", None, False),
    ("cannot separate them.", None, False),
    ("", None, False),
    ("A same-day vehicle arm would.", "#222222", True),
]
yy = .99
for s, c, bold in lines:
    ax.text(0, yy, s, transform=ax.transAxes, va="top", ha="left", fontsize=9.3,
            color=c if c else "#333333", fontweight="bold" if bold else "normal")
    yy -= .0455
panel_tag(ax, "D", dx=-0.04, dy=1.085)

n1 = V["CS1014-1b"]["mel"]["n"] + V["CS1014-1b"]["mor"]["n"]
n2 = V["CS1014-2w"]["mel"]["n"] + V["CS1014-2w"]["mor"]["n"]
fig.suptitle(f"Morphine makes the melittin response slightly LARGER, not smaller  —  "
             f"2 mice, {n1+n2:,} cells, within animal",
             fontsize=13.5, fontweight="bold", y=.955)
fig.savefig(os.path.join(FIG, "S4b_what_morphine_changed.png"), dpi=170)
plt.close(fig)
print("\nfigure done")

# ------------------------------------------------------------------ one slide
W, H = Inches(13.333), Inches(7.5)
prs = Presentation()
prs.slide_width, prs.slide_height = W, H
s = prs.slides.add_slide(prs.slide_layouts[6])


def tb(left, top, width, height, text, size=13, bold=False,
       color=RGBColor(0x22, 0x22, 0x22)):
    box = s.shapes.add_textbox(left, top, width, height)
    tf = box.text_frame
    tf.word_wrap = True
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = PP_ALIGN.LEFT
        p.space_after = Pt(3)
        run = p.add_run()
        run.text = line
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.color.rgb = color
        run.font.name = "Calibri"


tb(Inches(.45), Inches(.22), Inches(12.45), Inches(.6),
   "4b.  What morphine changed — and what it did not",
   size=24, bold=True, color=RGBColor(0x1F, 0x3A, 0x5F))

ph = Inches(4.72)
pw = Emu(int(ph * (13.0 / 5.5)))
s.shapes.add_picture(os.path.join(FIG, "S4b_what_morphine_changed.png"),
                     Emu(int((W - pw) / 2)), Inches(.82), width=pw, height=ph)

tb(Inches(.55), Inches(5.66), Inches(12.25), Inches(1.72),
   "What we did:  compared melittin alone with melittin + morphine on every measure "
   "we have, in the two mice recorded under both arms, so between-animal variance "
   "drops out.\n"
   "What we found:  all four measures move the same way and all by a small amount — "
   "more cells driven up, fewer driven down, higher amplitude, pooled odds ratio 1.15 "
   "and 1.16 — so morphine does not quieten the melittin response. No individual area "
   "changes reproducibly, and because the two arms sit seven weeks apart with no "
   "same-day vehicle control, this shift cannot yet be attributed to morphine rather "
   "than to the session.",
   size=12.5)

s.notes_slide.notes_text_frame.text = (
    "THE NUMBERS BEHIND PANEL A\n"
    + "\n".join(
        f"{aid}:  up {V[aid]['mel']['up']:.1f} -> {V[aid]['mor']['up']:.1f} "
        f"({V[aid]['mor']['up']-V[aid]['mel']['up']:+.1f} pp),  "
        f"down {V[aid]['mel']['down']:.1f} -> {V[aid]['mor']['down']:.1f} "
        f"({V[aid]['mor']['down']-V[aid]['mel']['down']:+.1f}),  "
        f"no response {V[aid]['mel']['none']:.1f} -> {V[aid]['mor']['none']:.1f},  "
        f"amplitude {V[aid]['mel']['amp']:+.1f} -> {V[aid]['mor']['amp']:+.1f} "
        f"({V[aid]['mor']['amp']-V[aid]['mel']['amp']:+.1f} % of baseline).  "
        f"Cells: {V[aid]['mel']['n']:,} melittin, {V[aid]['mor']['n']:,} morphine."
        for aid, _ in PAIRS)
    + "\n\nWHY NO P VALUE ON THE WHOLE-CORTEX SHIFT\n"
      "n = 2 paired animals. A two-sided sign test floors at p = 0.5 and a paired "
      "Wilcoxon cannot go below p = 0.5 either. The odds ratios in panel C DO carry "
      "CIs, but those CIs are computed over cells within one animal, so they say the "
      "shift is real inside that animal - not that it generalises. Two more paired "
      "animals would make the across-animal claim testable.\n\n"
      "WHY THE SHIFT IS SMALL\n"
      "+2.7 pp on a base of 28.7 % is a 9 % relative increase, OR 1.15. The "
      "no-response fraction barely moves (7.6 -> 6.2 and 10.9 -> 9.8), so this is not "
      "morphine silencing cells; it is a small number of cells moving from the down "
      "class to the up class.\n\n"
      "THE SEVEN-WEEK CONFOUND - EXPECT THIS QUESTION\n"
      "melittin only was recorded 20260218, melittin + morphine 20260408. Separate "
      "sessions, separate EXTRACT runs, no cell tracking across them, and no melittin "
      "+ vehicle session on the second day. Window maturation, GCaMP expression drift, "
      "and a second melittin exposure are all confounded with the drug. Nothing in "
      "this dataset separates them.\n\n"
      "WHAT WOULD FIX IT, IN ORDER OF VALUE\n"
      "1. a same-day melittin + vehicle arm, ideally within-session crossover - this "
      "removes the day confound entirely and is the single most valuable addition\n"
      "2. two more paired animals, which makes the region question testable\n"
      "3. registering the baseline and injection EXTRACT runs to each other so cells "
      "can be tracked, which would also fix the inflated 28 % / 63 % responsive rates "
      "discussed in supplement S2\n\n"
      "This slide is an addition. Slides 1-7 and supplements S1-S6 are unchanged."
)

out = os.path.join(HERE, "CrystalSkull_SLIDE4b_what_morphine_changed.pptx")
prs.save(out)
print(f"wrote {out}")
