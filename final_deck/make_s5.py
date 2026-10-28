"""Slide 5 figure: five visual panels only.

The written argument (what we ran, what happened, why, what to change) belongs in
the slide body and the speaker notes. The first version tried to carry it inside the
figure and the text ran off both edges.
"""
from __future__ import annotations

import csv
import glob
import os

import cv2
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

from make_final_figures import (AMBER, BEH, GREY, OUT, RED, TEAL, panel_tag)

GDIR = os.path.join(os.path.expanduser("~"), "Downloads", "grimace")


def slide5():
    fig = plt.figure(figsize=(13.0, 5.5))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.10, 1.05, 1.05],
                          height_ratios=[1.02, 1], wspace=.38, hspace=.66,
                          left=.145, right=.985, top=.80, bottom=.13)

    man = list(csv.DictReader(open(os.path.join(BEH,
                                                "behavior_calcium_manifest.csv"))))

    # --- A  inventory -------------------------------------------------------
    ax = fig.add_subplot(gs[:, 0])
    order = sorted(man, key=lambda r: (r["mouse"], r["date"], r["condition"],
                                       r["repeat"]))
    col = {"melittin_only": TEAL, "melittin_morphine": AMBER}
    for i, r in enumerate(order):
        n = int(r["n_frames"]) if r["n_frames"] else 0
        ax.barh(i, n / 35 / 60, color=col.get(r["arm"], GREY),
                hatch="" if r["media"] == "mp4" else "///",
                edgecolor="white", height=.76)
    ax.set_yticks(np.arange(len(order)))
    ax.set_yticklabels(
        [f"{r['mouse']}  {r['date'][4:6]}/{r['date'][6:]}  "
         f"{'base' if r['condition'].startswith('Base') else 'inj'}{r['repeat']}"
         for r in order], fontsize=7.0)
    ax.invert_yaxis()
    ax.set_xlabel("session length (min, at 35 fps)")
    n_cal = sum(1 for r in order if r["calcium_exists"] == "yes")
    ax.set_title(f"{len(order)} face-video sessions indexed\n"
                 f"{n_cal} have a matching calcium folder", fontsize=9.8, pad=7)
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=TEAL),
                       Rectangle((0, 0), 1, 1, color=AMBER),
                       Rectangle((0, 0), 1, 1, fc="white", ec=GREY, hatch="///")],
              labels=["melittin only", "melittin + morphine", "loose .jpg frames"],
              frameon=False, fontsize=7.4, loc="lower right")
    panel_tag(ax, "A", dx=-0.52, dy=1.11)

    # --- B,C  representative frames ----------------------------------------
    cap = cv2.VideoCapture(os.path.join(GDIR, "extract", "GrimACE_Live_Demo",
                                        "videos", "a3_c_0h_front_4.avi"))
    cap.set(cv2.CAP_PROP_POS_FRAMES, 120)
    ok, theirs = cap.read()
    cap.release()
    theirs = (cv2.cvtColor(theirs, cv2.COLOR_BGR2GRAY).astype(float) / 255
              if ok else None)

    ours = plt.imread(sorted(glob.glob(os.path.join(
        BEH, "frames", "1014m01b__20260408__Baseline", "*.png")))[2])
    if ours.ndim == 3:
        ours = ours[..., :3].mean(axis=2)

    for k, (im, ttl, c, tag, tdx) in enumerate([
            (theirs, "GrimACE rig  (their example video)", TEAL, "B", -.30),
            (ours, "our crystal-skull rig", AMBER, "C", -.16)]):
        ax = fig.add_subplot(gs[0, 1 + k])
        if im is not None:
            ax.imshow(im, cmap="gray", vmin=0, vmax=1)
            ax.set_title(f"{ttl}\npixel mean {im.mean():.2f},  SD {im.std():.2f}",
                         fontsize=9.2, color=c, fontweight="bold", pad=6)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color(c); sp.set_linewidth(2.5); sp.set_visible(True)
        panel_tag(ax, tag, dx=tdx, dy=1.26)

    # --- D  intensity histogram --------------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    bins = np.linspace(0, 1, 60)
    ax.hist(theirs.ravel(), bins=bins, density=True, color=TEAL, alpha=.60,
            label=f"GrimACE rig   SD {theirs.std():.2f}")
    ax.hist(ours.ravel(), bins=bins, density=True, color=AMBER, alpha=.60,
            label=f"our rig   SD {ours.std():.2f}")
    ax.axvline(.5175, color=RED, ls="--", lw=1.6)
    ax.annotate("scorer is hard-coded\nto mean 0.5175, SD 0.14",
                (.5175, ax.get_ylim()[1] * .62), xytext=(-8, 0),
                textcoords="offset points", fontsize=7.6, color=RED,
                ha="right", va="center")
    ax.set_xlabel("pixel intensity  (0 = black, 1 = white)")
    ax.set_ylabel("density")
    ax.set_title("Our frames are bimodal: dark background +\n"
                 f"bright mouse, SD {ours.std()/theirs.std():.1f}× theirs",
                 fontsize=9.4, pad=6)
    ax.legend(frameon=False, fontsize=7.8, loc="upper right")
    panel_tag(ax, "D", dx=-0.19, dy=1.22)

    # --- E  face-cropper acceptance ----------------------------------------
    ax = fig.add_subplot(gs[1, 2])
    gd = os.path.join(GDIR, "GrimACE_manuscript", "GrimaceRecorder")
    g_ctrl = list(csv.DictReader(open(os.path.join(
        gd, "a3_c_0h_front_4_grimace.csv"))))
    g_ours = list(csv.DictReader(open(os.path.join(gd, "ours_base_raw.csv"))))
    nc = sum(1 for r in g_ctrl if r["crop"] == "yolo")
    no = sum(1 for r in g_ours if r["crop"] == "yolo")
    ax.bar([0, 1], [100 * nc / len(g_ctrl), 100 * no / len(g_ours)],
           color=[TEAL, AMBER], width=.52)
    for x, (a, b) in enumerate([(nc, len(g_ctrl)), (no, len(g_ours))]):
        ax.text(x, 100 * a / b + 2.5, f"{a} / {b} frames", ha="center",
                fontsize=9.2, fontweight="bold")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["GrimACE rig\n(their video)",
                        "our rig\n(1014m01b baseline)"], fontsize=8.6)
    ax.set_ylim(0, 68)
    ax.set_ylabel("frames accepted by the cropper (%)")
    ax.set_title("Same weights, same code, same sampling:\n"
                 "the cropper never fires on our footage",
                 fontsize=9.4, pad=6)
    panel_tag(ax, "E", dx=-0.26, dy=1.22)

    fig.suptitle("Slide 5  |  Face video is indexed and synced; automated grimace "
                 "scoring needs a rig change, not a new model",
                 fontsize=13.5, fontweight="bold", y=.965)
    fig.savefig(os.path.join(OUT, "S5_behavior.png"), dpi=170)
    plt.close(fig)
    print("S5 done")


if __name__ == "__main__":
    slide5()
