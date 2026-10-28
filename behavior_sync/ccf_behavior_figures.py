"""Figures for the facial-behavior arm: inventory, method, and condition comparison.

  figB1  what data exists - 19 sessions, which drug arm, how long, calcium matched
  figB2  how automated grimace scoring works, drawn on one of our own frames
  figB3  baseline vs melittin vs melittin+morphine, real frames, within animal

figB2 is a METHOD schematic. The action-unit boxes are drawn at nominal positions to show
what the algorithm looks at - they are not model output, and no grimace score is computed
anywhere in this file. Scoring needs GrimACE or PainFace; see README_behavior.md.
"""
from __future__ import annotations

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
FRAMES = os.path.join(HERE, "frames")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

NAVY, TEAL, AMBER, RED, GREY = "#1F3A5F", "#0D5C63", "#D08C1E", "#C0392B", "#555555"
ARM_COL = {"melittin_only": TEAL, "melittin_morphine": AMBER}

rows = list(csv.DictReader(open(os.path.join(HERE, "behavior_calcium_manifest.csv"))))


def load(tag, k=0):
    d = os.path.join(FRAMES, tag)
    if not os.path.isdir(d):
        return None
    fs = sorted(f for f in os.listdir(d) if f.lower().endswith((".png", ".jpg")))
    if not fs:
        return None
    return plt.imread(os.path.join(d, fs[min(k, len(fs) - 1)]))


# ============================================================ figB1 inventory
fig, axes = plt.subplots(1, 2, figsize=(16, 7.5),
                         gridspec_kw=dict(width_ratios=[1.35, 1]))

ax = axes[0]
order = sorted(rows, key=lambda r: (r["mouse"], r["date"], r["condition"], r["repeat"]))
y = np.arange(len(order))
for i, r in enumerate(order):
    n = int(r["n_frames"]) if r["n_frames"] else 0
    mins = n / 35.0 / 60.0
    c = ARM_COL.get(r["arm"], GREY)
    hatch = "" if r["media"] == "mp4" else "//"
    ax.barh(i, mins, color=c, alpha=.85, edgecolor="white", hatch=hatch)
    ax.text(mins + .4, i, f"{n:,} fr", va="center", fontsize=7.5, color=GREY)
ax.set_yticks(y)
ax.set_yticklabels([f"{r['mouse']}  {r['date']}  {r['condition']}{r['repeat']}"
                    for r in order], fontsize=8.5)
ax.invert_yaxis()
ax.set_xlabel("session length (min, at 35 fps)")
ax.set_title("19 face-video sessions on the share", fontsize=12, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)
ax.legend(handles=[Rectangle((0, 0), 1, 1, color=ARM_COL["melittin_only"]),
                   Rectangle((0, 0), 1, 1, color=ARM_COL["melittin_morphine"]),
                   Rectangle((0, 0), 1, 1, fc="white", ec=GREY, hatch="//")],
          labels=["melittin only", "melittin + morphine", "loose .jpg (needs encoding)"],
          frameon=False, fontsize=9, loc="lower right")

ax = axes[1]
ax.axis("off")
mice = {}
for r in order:
    mice.setdefault(r["mouse"], set()).add(r["arm"])
both = [m for m, a in mice.items() if len(a) > 1]
n_mp4 = sum(1 for r in order if r["media"] == "mp4")
n_ts = sum(1 for r in order if r["has_timestamps"] == "yes")
n_cal = sum(1 for r in order if r["calcium_exists"] == "yes")
txt = "\n".join([
    "INVENTORY",
    "",
    f"  sessions              {len(order)}",
    f"  mice                  {len(mice)}  ({', '.join(sorted(mice))})",
    f"  already .mp4          {n_mp4}   (20260408)",
    f"  loose .jpg sequences  {len(order)-n_mp4}   (20260218, 20260219)",
    f"  per-frame timestamps  {n_ts} sessions have a _behavior_timestamp.csv;",
    "                        the .jpg sets carry the clock in the filename",
    f"  calcium folder found  {n_cal} / {len(order)}",
    "",
    "SYNC",
    "",
    "  Frame index x 1/35 s is NOT the wall clock. 1014m01b 20260408",
    "  Injection holds 1,504 s of frames spread over 3,761 s of wall",
    "  clock, so that recording is segmented. Always use the",
    "  timestamps, never the frame index.",
    "",
    "WORTH CONFIRMING WITH OMER",
    "",
    f"  {', '.join(both)} have face video under BOTH drug arms",
    "  (melittin only on 20260218, melittin + morphine on 20260408).",
    "  If the calcium is from the same animals, the melittin vs",
    "  melittin+morphine comparison is WITHIN animal, not the n = 4",
    "  vs 2 unpaired design the current analysis assumes.",
    "",
    "  Folder names are not proof of animal identity - ask before",
    "  re-running the statistics on that basis.",
])
ax.text(0, .99, txt, va="top", ha="left", fontsize=10.2, family="monospace",
        color="#333333", linespacing=1.5, transform=ax.transAxes)

fig.suptitle("Fig B1  |  Facial behavior recordings paired with the calcium sessions",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .94])
fig.savefig(os.path.join(OUT, "figB1_behavior_inventory.png"), dpi=160)
plt.close(fig)
print("figB1 done")

# ============================================================ figB2 method
base = load("1014m01b__20260408__Baseline", 3)
if base is not None:
    if base.ndim == 3:
        base = base.mean(axis=2)
    h, w = base.shape

    fig, axes = plt.subplots(1, 4, figsize=(19, 5.6))

    axes[0].imshow(base, cmap="gray")
    axes[0].set_title("1. Raw frame\n1024x768, 35 fps, IR-lit", fontsize=11, fontweight="bold")

    # frame-quality / selector step
    axes[1].imshow(base, cmap="gray")
    axes[1].set_title("2. Frame selection\nkeep frames where the face is front-on\n"
                      "(GrimACE: SLEAP nose/ear keypoints + a\ndensenet 'suitable for MGS' classifier)",
                      fontsize=10, fontweight="bold")

    # crop
    cy, cx = int(h * .42), int(w * .55)
    hh, ww = int(h * .46), int(w * .40)
    axes[2].imshow(base, cmap="gray")
    axes[2].add_patch(Rectangle((cx - ww, cy - hh), 2 * ww, 2 * hh,
                                fill=False, ec=AMBER, lw=3))
    axes[2].set_title("3. Face crop\nautomatic bounding box", fontsize=11, fontweight="bold")

    # action units
    axes[3].imshow(base, cmap="gray")
    # (fx, fy, half-width, half-height, label offset) read off this actual frame -
    # an earlier version used guessed fractions and the boxes landed on the light rod
    # and on empty background, which is worse than no boxes on a methods figure.
    aus = [("orbital tightening", .395, .105, .075, .075, (-.02, -.09)),
           ("nose bulge", .545, .280, .070, .075, (.12, -.02)),
           ("ear position", .900, .195, .060, .075, (-.13, -.09)),
           ("whisker change", .225, .350, .105, .090, (.02, .12)),
           ("cheek bulge", .700, .400, .085, .080, (.12, .04))]
    for name, fx, fy, hw, hh, (ox, oy) in aus:
        x, yy = fx * w, fy * h
        axes[3].add_patch(Rectangle((x - hw * w, yy - hh * h), 2 * hw * w, 2 * hh * h,
                                    fill=False, ec=RED, lw=2.2))
        lx, ly = (fx + ox) * w, (fy + oy) * h
        axes[3].annotate(name, (lx, ly), color=RED, fontsize=9, fontweight="bold",
                         ha="center", va="center",
                         bbox=dict(fc="black", alpha=.55, ec="none", pad=1.5))
    axes[3].set_title("4. Action units scored 0/1/2\nsummed -> MGS 0-8 (PainFace, 4 AUs)\n"
                      "or 0-10 (5 AUs)", fontsize=10, fontweight="bold")

    for a in axes:
        a.set_xticks([]); a.set_yticks([])
    for i in range(3):
        fig.add_artist(FancyArrowPatch(
            (.253 + i * .2385, .5), (.283 + i * .2385, .5),
            transform=fig.transFigure, arrowstyle="-|>", mutation_scale=22,
            color=NAVY, lw=2.5))

    fig.suptitle("Fig B2  |  How automated grimace scoring works, drawn on one of our own frames\n"
                 "SCHEMATIC: the boxes mark what the algorithm looks at. They are not model "
                 "output and no score is computed here.",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, .86])
    fig.savefig(os.path.join(OUT, "figB2_method_steps.png"), dpi=160)
    plt.close(fig)
    print("figB2 done")

# ============================================================ figB3 conditions
panels = [
    ("1014m01b", [("1014m01b__20260218__Baseline", "Baseline\n(20260218)"),
                  ("1014m01b__20260218__Injection1", "Melittin\n(20260218)"),
                  ("1014m01b__20260408__Baseline", "Baseline\n(20260408)"),
                  ("1014m01b__20260408__Injection", "Melittin + morphine\n(20260408)")]),
    ("1014m02w", [("1014m02w__20260218__Baseline", "Baseline\n(20260218)"),
                  ("1014m02w__20260218__Injection3", "Melittin\n(20260218)"),
                  ("1014m02w__20260408__Baseline", "Baseline\n(20260408)"),
                  ("1014m02w__20260408__Injection", "Melittin + morphine\n(20260408)")]),
    ("CS0204-1b", [("CS0204-1b__20260219__Baseline", "Baseline\n(20260219)"),
                   ("CS0204-1b__20260219__Injection1", "Melittin\n(20260219)"),
                   (None, ""), (None, "")]),
    ("CS0204-2w", [("CS0204-2w__20260219__Baseline", "Baseline\n(20260219)"),
                   ("CS0204-2w__20260219__Injection1", "Melittin\n(20260219)"),
                   (None, ""), (None, "")]),
]
fig, axes = plt.subplots(len(panels), 4, figsize=(15, 3.6 * len(panels)))
for r, (mouse, cells) in enumerate(panels):
    for c, (tag, lab) in enumerate(cells):
        ax = axes[r, c]
        ax.set_xticks([]); ax.set_yticks([])
        if tag is None:
            ax.axis("off")
            if c == 2:
                ax.text(.5, .5, "no morphine session\nfor this mouse", ha="center",
                        va="center", fontsize=10, color=GREY, style="italic",
                        transform=ax.transAxes)
            continue
        im = load(tag, 4)
        if im is None:
            ax.text(.5, .5, "missing", ha="center", va="center", transform=ax.transAxes)
            continue
        ax.imshow(im if im.ndim == 2 else im.mean(axis=2), cmap="gray")
        col = AMBER if "morphine" in lab.lower() else (TEAL if "Melittin" in lab else GREY)
        ax.set_title(lab, fontsize=10, fontweight="bold", color=col)
        for s in ax.spines.values():
            s.set_color(col); s.set_linewidth(2.5)
    axes[r, 0].set_ylabel(mouse, fontsize=12, fontweight="bold", color=NAVY)
fig.suptitle("Fig B3  |  Same animals, three conditions - raw face frames\n"
             "1014m01b and 1014m02w were recorded under BOTH melittin alone and "
             "melittin + morphine, so their comparison can be made within animal",
             fontsize=13.5, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .93])
fig.savefig(os.path.join(OUT, "figB3_conditions.png"), dpi=150)
plt.close(fig)
print("figB3 done")
print("->", OUT)
