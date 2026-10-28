"""Figures for the PI deck: the paper's theory, the transfer attempt, and the corrected results.

Numbers here are my own verified measurements unless the caption says otherwise.
"""
import glob
import os
import zipfile

import cv2
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
F = os.path.join(HERE, "figures")
os.makedirs(F, exist_ok=True)
NAVY, TEAL, AMBER, RED, GREY, OK = "#1F3A5F", "#0D5C63", "#D08C1E", "#C0392B", "#555555", "#1E8449"
DATA = r"C:/Users/hsollim/Downloads/grimace/Data.zip"

# ------------------------------------------------------------------ C1 theory
fig, ax = plt.subplots(figsize=(15, 7.2))
ax.axis("off")
ax.text(.5, .97, "Mouse Grimace Scale (Langford 2010) and how GrimACE automates it",
        ha="center", fontsize=15.5, fontweight="bold", color=NAVY, transform=ax.transAxes)
AUS = [("Orbital tightening", "eye narrowing / squeezing shut"),
       ("Nose bulge", "rounded bulge on the nose bridge"),
       ("Cheek bulge", "rounded bulge in the cheek"),
       ("Ear position", "ears rotate back and apart"),
       ("Whisker change", "whiskers pulled back or fanned forward")]
ax.text(.03, .86, "THEORY - five facial action units, each graded 0 / 1 / 2",
        fontsize=13, fontweight="bold", color=TEAL, transform=ax.transAxes)
for i, (n, d) in enumerate(AUS):
    y = .78 - i * .073
    ax.add_patch(Rectangle((.03, y - .020), .016, .040, color=TEAL, transform=ax.transAxes))
    ax.text(.055, y, n, fontsize=12, fontweight="bold", va="center", transform=ax.transAxes)
    ax.text(.26, y, d, fontsize=11, color=GREY, va="center", transform=ax.transAxes)
ax.text(.03, .37, "Sum of the five  ->  MGS 0-10 manual.   GrimACE grades each 0-3  ->  0-15.",
        fontsize=11.5, style="italic", color=GREY, transform=ax.transAxes)
ax.text(.03, .28,
        "Why automate: manual scoring needs a trained expert, takes hours per animal, and\n"
        "scorers disagree with each other. GrimACE replaces the human with four networks.",
        fontsize=11.5, color="#333333", transform=ax.transAxes)
STEPS = [("frame_quality", "mobilenetv3\n1 channel", "is this frame\nscorable at all?"),
         ("face_cropper", "YOLO\ndetector", "where is\nthe head?"),
         ("grimace_scorer", "ViT-B/16\n+ 5 heads", "grade each\naction unit"),
         ("top pose", "YOLO pose", "posture /\nmovement")]
ax.text(.03, .175, "GrimACE ARCHITECTURE", fontsize=13, fontweight="bold", color=TEAL,
        transform=ax.transAxes)
for i, (nm, arch, job) in enumerate(STEPS):
    x = .05 + i * .235
    ax.add_patch(Rectangle((x, .02), .19, .115, fc="#EAF2F4", ec=TEAL, lw=2, transform=ax.transAxes))
    ax.text(x + .095, .112, nm, ha="center", fontsize=11.5, fontweight="bold", color=NAVY,
            transform=ax.transAxes)
    ax.text(x + .095, .078, arch, ha="center", fontsize=9.5, color=TEAL, transform=ax.transAxes)
    ax.text(x + .095, .040, job, ha="center", fontsize=9, color=GREY, transform=ax.transAxes)
    if i < 3:
        ax.add_patch(FancyArrowPatch((x + .195, .078), (x + .232, .078), transform=ax.transAxes,
                                     arrowstyle="-|>", mutation_scale=18, color=NAVY, lw=2))
fig.savefig(f"{F}/figC1_theory.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("C1")

# --------------------------------------------------- C2 normalisation mismatch
fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
ax[0].bar([0, 1], [0.581, 0.502], yerr=[0.037, 0.024], capsize=6, color=[TEAL, RED], width=.55)
ax[0].axhline(0.5175, color="k", ls="--", lw=2)
ax[0].text(1.42, 0.5175, "scorer expects\n0.5175", fontsize=9.5, va="center")
ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["their training\ncrops (n=45)", "OUR crops\n(n=12)"])
ax[0].set_ylabel("mean intensity"); ax[0].set_ylim(0, .78)
ax[0].set_title("Mean - close enough", fontsize=12.5, fontweight="bold", color=OK)
ax[0].spines[["top", "right"]].set_visible(False)
ax[1].bar([0, 1], [0.137, 0.256], yerr=[0.014, 0.016], capsize=6, color=[TEAL, RED], width=.55)
ax[1].axhline(0.140, color="k", ls="--", lw=2)
ax[1].text(1.42, 0.140, "scorer expects\n0.140", fontsize=9.5, va="center")
ax[1].set_xticks([0, 1]); ax[1].set_xticklabels(["their training\ncrops", "OUR crops"])
ax[1].set_ylabel("intensity SD (contrast)"); ax[1].set_ylim(0, .33)
ax[1].set_title("Contrast - 1.87x TOO HIGH", fontsize=12.5, fontweight="bold", color=RED)
ax[1].spines[["top", "right"]].set_visible(False)
fig.suptitle("The scorer hard-codes  (x - 0.5175) / 0.140\n"
             "A brightness gain moves mean and SD together, so it can never hit both targets",
             fontsize=13.5, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .85])
fig.savefig(f"{F}/figC2_normalisation.png", dpi=150)
plt.close(fig)
print("C2")

# ------------------------------------------------------- C3 threshold artifact
fig, ax = plt.subplots(figsize=(11, 5.4))
c = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.415, 0.0, 0.445]
ax.bar(range(len(c)), c, color=[OK if v > .25 else RED for v in c], width=.7)
ax.axhline(0.9, color=RED, lw=2.5)
ax.text(9.6, 0.9, "hard-coded\ngate = 0.90", color=RED, fontsize=10, va="center", fontweight="bold")
ax.axhline(0.25, color=GREY, ls="--", lw=1.8)
ax.text(9.6, 0.25, "YOLO default\n0.25", color=GREY, fontsize=9.5, va="center")
ax.set_xlabel("our frames"); ax.set_ylabel("YOLO face confidence"); ax.set_ylim(0, 1)
ax.set_title("The cropper is not blind - it sits under a gate tuned for their rig\n"
             "2/10 of our frames clear YOLO's own default; none clear the 0.90 the app hard-codes\n"
             "the source carries the comment: TODO Decide on confidence threshold",
             fontsize=12, fontweight="bold")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{F}/figC3_threshold.png", dpi=150)
plt.close(fig)
print("C3")

# --------------------------------------------------- C4 standardisation effect
fig, ax = plt.subplots(figsize=(9.8, 5.6))
g = [0.8, 1.0, 1.2, 1.4, 1.6]
ax.plot(g, [6.5, 3.5, 1.5, 1.0, 1.0], "-o", color=RED, lw=2.5, ms=9,
        label="scalar gain, my first attempt - swing 5.5")
ax.plot(g, [2.5, 2.5, 2.17, 1.0, 0.5], "-o", color=AMBER, lw=2.5, ms=9,
        label="per-image standardise - swing 2.0")
ax.plot(g, [2.3, 2.3, 2.0, 3.2, 4.0], "-o", color=TEAL, lw=2.5, ms=9,
        label="their own data, reference - swing 2.0")
ax.set_xlabel("brightness gain applied before scoring"); ax.set_ylabel("mean MGS")
ax.set_title("Fixing the normalisation removes most of the instability\n"
             "but does not by itself make the measurement usable",
             fontsize=12.5, fontweight="bold")
ax.legend(frameon=False, fontsize=9.5)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(f"{F}/figC4_standardise.png", dpi=150)
plt.close(fig)
print("C4")

# ----------------------------------------------------------- C5 ears in frame
our = cv2.imread(sorted(glob.glob(f"{HERE}/frames/1014m01b__20260408__Baseline/*.png"))[3],
                 cv2.IMREAD_GRAYSCALE)
their = None
with zipfile.ZipFile(DATA) as z:
    mgs = [i for i in z.namelist() if i.startswith("Data/MGS_scores/images") and i.endswith(".png")]
    for i in mgs[::97]:
        with z.open(i) as fh:
            a = cv2.imdecode(np.frombuffer(fh.read(), np.uint8), cv2.IMREAD_GRAYSCALE)
        if a is not None and 300 < a.shape[0] < 520:
            their = a
            break
fig, ax = plt.subplots(1, 2, figsize=(13.5, 6.2))
ax[0].imshow(their, cmap="gray")
ax[0].set_title("THEIR training crop\nboth ears in frame, eyes clear of the edge",
                fontsize=12.5, fontweight="bold", color=OK)
ax[1].imshow(our, cmap="gray")
h, w = our.shape
ax[1].add_patch(Rectangle((0, 0), w, 58, fc=RED, alpha=.32))
ax[1].text(w / 2, 30, "EARS ARE ABOVE THIS EDGE - not in frame", ha="center", va="center",
           color="white", fontsize=11.5, fontweight="bold")
for x in (400, 715):
    ax[1].add_patch(Rectangle((x - 55, 18), 110, 82, fill=False, ec=AMBER, lw=2.5, ls="--"))
ax[1].text(557, 132, "eyes clipped at the top edge", ha="center", color=AMBER,
           fontsize=11, fontweight="bold")
ax[1].set_title("OUR frame\ncamera sits below-front", fontsize=12.5, fontweight="bold", color=RED)
for a in ax:
    a.set_xticks([]); a.set_yticks([])
fig.suptitle("The real mechanism: 2 of the 5 action units are physically unobservable in our footage\n"
             "the model grades ear position 3/3 on frames that contain no ears",
             fontsize=13.5, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .87])
fig.savefig(f"{F}/figC5_ears.png", dpi=150)
plt.close(fig)
print("C5")

# --------------------------------------------------------- C6 what to change
fig, ax = plt.subplots(figsize=(14, 6.4))
ax.axis("off")
ax.text(.5, .96, "What has to change for grimace scoring to work here",
        ha="center", fontsize=15, fontweight="bold", color=NAVY, transform=ax.transAxes)
rows = [("Camera angle", "below-front; ears out of frame, eyes clipped",
         "re-aim / widen so the whole head is in frame", "rig", OK),
        ("Illumination", "hard IR side-light, contrast SD 0.256",
         "soften and diffuse; target SD near 0.14", "rig", OK),
        ("Background", "dark, mouse is bright on black",
         "light background, as in both validated setups", "rig", OK),
        ("Cropper gate", "hard-coded 0.90, tuned for their rig",
         "lower it, or skip the cropper - our prep is head-fixed", "one line", OK),
        ("Normalisation", "scalar gain cannot fix mean and SD together",
         "per-image standardise to (0.5175, 0.140)", "one line", OK),
        ("Scorer weights", "trained on their animals and rig",
         "only retrain if the rig fixes above are not enough", "big job", RED)]
ax.text(.03, .855, "what", fontsize=11, fontweight="bold", color=GREY, transform=ax.transAxes)
ax.text(.17, .855, "problem now", fontsize=11, fontweight="bold", color=GREY, transform=ax.transAxes)
ax.text(.50, .855, "fix", fontsize=11, fontweight="bold", color=GREY, transform=ax.transAxes)
ax.text(.88, .855, "cost", fontsize=11, fontweight="bold", color=GREY, transform=ax.transAxes)
for i, (a_, b_, c_, d_, col) in enumerate(rows):
    y = .77 - i * .115
    ax.add_patch(Rectangle((.02, y - .045), .96, .095, fc="#F7F7F7" if i % 2 else "white",
                           ec="none", transform=ax.transAxes))
    ax.text(.03, y, a_, fontsize=11.5, fontweight="bold", va="center", transform=ax.transAxes)
    ax.text(.17, y, b_, fontsize=10.5, color=GREY, va="center", transform=ax.transAxes)
    ax.text(.50, y, c_, fontsize=10.5, color="#333333", va="center", transform=ax.transAxes)
    ax.text(.88, y, d_, fontsize=10.5, fontweight="bold", color=col, va="center",
            transform=ax.transAxes)
ax.text(.5, .04, "Five of the six are rig or one-line changes. Only the last is an ML project,\n"
        "and it is only needed if fixing the rig is not enough.",
        ha="center", fontsize=12, fontweight="bold", color=NAVY, transform=ax.transAxes)
fig.savefig(f"{F}/figC6_what_to_change.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("C6")
print("done")
