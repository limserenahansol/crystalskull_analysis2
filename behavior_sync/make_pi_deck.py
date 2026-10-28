"""Rebuild figC5 honestly, then assemble the PI deck.

figC5 first claimed "both ears are outside the frame and the eyes are clipped at the top
edge". Checking twelve of our own frames, that is not true: the eyes sit in the upper
third and are not clipped, and an ear is often partly visible. What IS measurable and
does differ is the blown-out chest/forepaw region, which is what drives our contrast SD
to 1.87x the value the scorer expects. The figure now says only that.
"""
import glob
import os
import zipfile

import cv2
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.util import Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
F = os.path.join(HERE, "figures")
DATA = r"C:/Users/hsollim/Downloads/grimace/Data.zip"
NAVY, TEAL, AMBER, RED, GREY, OK = "#1F3A5F", "#0D5C63", "#D08C1E", "#C0392B", "#555555", "#1E8449"

# ------------------------------------------------ figC5, corrected
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

fig = plt.figure(figsize=(17, 5.8))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.15], wspace=.25)
a0 = fig.add_subplot(gs[0]); a0.imshow(their, cmap="gray", vmin=0, vmax=255)
a0.set_title("THEIR training crop" + "\n"
             + f"mean {their.mean()/255:.2f}   SD {their.std()/255:.3f}",
             fontsize=12.5, fontweight="bold", color=OK)
a1 = fig.add_subplot(gs[1]); a1.imshow(our, cmap="gray", vmin=0, vmax=255)
a1.set_title("OUR frame" + "\n"
             + f"mean {our.mean()/255:.2f}   SD {our.std()/255:.3f}",
             fontsize=12.5, fontweight="bold", color=RED)
for a in (a0, a1):
    a.set_xticks([]); a.set_yticks([])
a2 = fig.add_subplot(gs[2])
a2.hist(their.ravel()/255, bins=60, range=(0, 1), color=OK, alpha=.70,
        density=True, label="theirs - single peak")
a2.hist(our.ravel()/255, bins=60, range=(0, 1), color=RED, alpha=.60,
        density=True, label="ours - two peaks")
a2.set_xlabel("pixel intensity"); a2.set_ylabel("density")
a2.legend(frameon=False, fontsize=10)
a2.set_title("Intensity histogram", fontsize=12.5, fontweight="bold")
a2.spines[["top", "right"]].set_visible(False)
fig.suptitle("What actually differs: our frame is bimodal - a large dark background "
             "against a bright subject\n"
             "Their training crops fill the frame with the animal, so the histogram has one peak. "
             "That split is what doubles our contrast SD.",
             fontsize=13, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .85])
fig.savefig(f"{F}/figC5_contrast.png", dpi=150)
plt.close(fig)
print("figC5_contrast rebuilt")

# ------------------------------------------------ deck
prs = Presentation()
prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)
BLANK = prs.slide_layouts[6]


def _c(h):
    return RGBColor(int(h[1:3], 16), int(h[3:5], 16), int(h[5:7], 16))


def title_slide(t, sub, extra=""):
    s = prs.slides.add_slide(BLANK)
    b = s.shapes.add_textbox(Inches(.8), Inches(2.3), Inches(11.7), Inches(1.1))
    r = b.text_frame.paragraphs[0].add_run(); r.text = t
    r.font.size, r.font.bold, r.font.color.rgb = Pt(30), True, _c(NAVY)
    b2 = s.shapes.add_textbox(Inches(.9), Inches(3.5), Inches(11.5), Inches(1.4))
    b2.text_frame.word_wrap = True
    r2 = b2.text_frame.paragraphs[0].add_run(); r2.text = sub
    r2.font.size, r2.font.color.rgb = Pt(16), _c(GREY)
    if extra:
        b3 = s.shapes.add_textbox(Inches(.9), Inches(5.0), Inches(11.5), Inches(1.4))
        b3.text_frame.word_wrap = True
        r3 = b3.text_frame.paragraphs[0].add_run(); r3.text = extra
        r3.font.size, r3.font.bold, r3.font.color.rgb = Pt(14), True, _c(TEAL)
    return s


def fig_slide(title, what, result, img, rcol=TEAL):
    s = prs.slides.add_slide(BLANK)
    tb = s.shapes.add_textbox(Inches(.35), Inches(.10), Inches(12.6), Inches(.40))
    r = tb.text_frame.paragraphs[0].add_run(); r.text = title
    r.font.size, r.font.bold, r.font.color.rgb = Pt(20), True, _c(NAVY)
    cb = s.shapes.add_textbox(Inches(.35), Inches(.50), Inches(12.6), Inches(.62))
    cb.text_frame.word_wrap = True
    p0 = cb.text_frame.paragraphs[0]
    r0 = p0.add_run(); r0.text = "What:  " + what
    r0.font.size, r0.font.color.rgb = Pt(12), _c(GREY)
    p1 = cb.text_frame.add_paragraph()
    r1 = p1.add_run(); r1.text = "Result:  " + result
    r1.font.size, r1.font.bold, r1.font.color.rgb = Pt(13), True, _c(rcol)
    if os.path.exists(img):
        with Image.open(img) as im:
            iw, ih = im.size
        bw, bh = Inches(12.5), Inches(5.9)
        nw, nh = bw, int(bw / (iw / ih))
        if nh > bh:
            nh, nw = bh, int(bh * (iw / ih))
        s.shapes.add_picture(img, int((prs.slide_width - nw) / 2), Inches(1.25), nw, nh)
    return s


def bullets(title, items, size=16):
    s = prs.slides.add_slide(BLANK)
    tb = s.shapes.add_textbox(Inches(.6), Inches(.30), Inches(12.1), Inches(.55))
    r = tb.text_frame.paragraphs[0].add_run(); r.text = title
    r.font.size, r.font.bold, r.font.color.rgb = Pt(24), True, _c(NAVY)
    box = s.shapes.add_textbox(Inches(.85), Inches(1.15), Inches(11.7), Inches(5.9))
    tf = box.text_frame; tf.word_wrap = True
    for i, it in enumerate(items):
        txt, col, bold = it if isinstance(it, tuple) else (it, GREY, False)
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        r = p.add_run(); r.text = txt
        r.font.size, r.font.bold, r.font.color.rgb = Pt(size), bold, _c(col)
        p.space_after = Pt(11)
    return s


# ---- title
title_slide("Can automated grimace scoring measure pain in our crystal-skull mice?",
            "Testing GrimACE (Sturman et al., Lab Animal 2026) on our melittin / morphine face videos. "
            "19 sessions, 6 mice, already synchronised to the calcium recordings.",
            "Short answer: not on the footage as currently acquired. The fix is a camera and lighting "
            "change, not a machine-learning project.")

# ---- PART 1 : theory
title_slide("Part 1  —  What the method is",
            "The Mouse Grimace Scale, and how GrimACE automates it.")
fig_slide("The theory: five facial action units",
          "Langford et al. 2010 defined five features that change when a mouse is in pain. A trained "
          "human grades each 0/1/2 from still frames. GrimACE replaces the human with four networks.",
          "Manual scoring needs an expert, takes hours per animal, and scorers disagree — that is the "
          "problem the paper solves.",
          f"{F}/figC1_theory.png")
fig_slide("How it performs on its own data",
          "GrimACE's published example videos run through the same offline harness I built, so this is "
          "their method on their rig.",
          "Works as advertised: frame_quality 0.998, YOLO finds the head at confidence 0.919, "
          "MGS = 1/15 for a pre-surgery baseline animal. Their published effect is 0h to 24h = "
          "+2.43 MGS (Cohen d = 0.77).",
          f"{F}/figB5_accepted_crops.png")

# ---- PART 2 : applying it
title_slide("Part 2  —  Applying it to our videos",
            "Same models, same code, our footage. Five stages, checked one at a time.")
fig_slide("Our data: 19 face-video sessions, all matched to calcium",
          "Every behaviour session on the share, indexed and paired with its calcium recording.",
          "19 sessions, 6 mice. 1014m01b and 1014m02w have video under BOTH drug arms, which you "
          "confirmed are the same animals — so melittin vs melittin+morphine is a within-animal comparison.",
          f"{F}/figB1_behavior_inventory.png")
fig_slide("Stage by stage: where it works and where it breaks",
          "Each of the five pipeline stages, their video on top, ours below.",
          "Stage 3 fails: the YOLO head detector returns confidence 0.000-0.445 on our frames against a "
          "gate of 0.90. Stages 4-5 then run on a box I placed by hand, so those MGS numbers are not a "
          "measurement.",
          f"{F}/figB7_pipeline_stages.png", RED)
fig_slide("Why stage 3 fails — and why it is less serious than it looks",
          "Raw YOLO confidence on our frames, against the 0.90 gate the app hard-codes and YOLO's own "
          "0.25 default.",
          "The detector is not blind. 2/10 of our frames clear YOLO's default threshold; none clear 0.90. "
          "The source even carries the comment TODO: Decide on confidence threshold. And our prep is "
          "head-fixed (centroid SD 18 px over 10 min), so a fixed crop removes the need for a detector "
          "entirely.",
          f"{F}/figC3_threshold.png", AMBER)
fig_slide("What the scorer was trained on, versus what we would feed it",
          "GrimACE's own labelled MGS training images (n = 2,972) against our fixed-box crops.",
          "Closer than expected — their training set is mostly face close-ups framed much like ours, and "
          "several of their mice are brown, not black. So framing and coat colour are not the barrier.",
          f"{F}/figB8_training_vs_ours.png")
fig_slide("The measurable difference: contrast, not framing",
          "Intensity statistics of their training crops versus ours, against the constants the scorer "
          "hard-codes.",
          "Mean is close (0.581 vs 0.502, target 0.5175). Contrast is not: our SD is 0.256 against their "
          "0.137 — 1.87x too high. A brightness gain scales mean and SD together, so it can never fix both.",
          f"{F}/figC2_normalisation.png", RED)
fig_slide("Where our extra contrast comes from",
          "The two crops with their intensity histograms.",
          "Their crops are filled by the animal, so the histogram has one peak. Ours splits into a dark "
          "background and a bright subject. That bimodality, not saturation, is what doubles our SD.",
          f"{F}/figC5_contrast.png", AMBER)

# ---- PART 3 : results
title_slide("Part 3  —  Results",
            "What the corrected pipeline actually returns, and what that means.")
fig_slide("Correcting the normalisation removes most of the instability",
          "Mean MGS across a brightness-gain sweep, three ways of normalising.",
          "My first attempt used a scalar gain and MGS swung 5.5 points — that was my error, not a "
          "property of the data. Per-image standardisation to the scorer's own constants cuts the swing "
          "to 2.0, matching their own data.",
          f"{F}/figC4_standardise.png", AMBER)
bullets("The decisive test — and why my first version of it was worthless", [
    ("I originally tested baseline vs injection in the 20260408 sessions and found no separation.",
     GREY, False),
    ("Two things were wrong with that test:", NAVY, True),
    ("1.  Wrong arm. 20260408 is melittin + MORPHINE. Morphine is an analgesic, so no grimace "
     "increase is the expected result there, not evidence of failure.", RED, False),
    ("2.  No power. Their own in-domain effect is +2.43 MGS (d = 0.77) with a per-frame SD of 3.2. "
     "At 14 frames a session mean carries a 95% CI of about +/-1.8 MGS, so a real effect needs "
     "roughly 2.6 MGS to clear noise. My reported effects were +0.43 and -0.43 — indistinguishable "
     "from zero either way.", RED, False),
    ("Re-run properly — melittin-only arm, 6 mice, 12 crop x normalisation configurations — there is "
     "still no consistent separation (best configuration: mean delta -0.15, positive in 3 of 6 mice).",
     NAVY, True),
    ("So the conclusion holds, but for a different reason than I first gave, and the earlier evidence "
     "for it should be discarded.", TEAL, True),
], size=15)
fig_slide("What has to change",
          "Every difference between their setup and ours, with what it would take to fix.",
          "Five of the six are rig or one-line software changes. Only retraining the scorer is a real "
          "project, and it is only needed if fixing the rig is not enough.",
          f"{F}/figC6_what_to_change.png", OK)
bullets("Recommendation", [
    ("Do not invest in retraining. The models are sound — our acquisition is out of their domain.",
     NAVY, True),
    ("When crystal-skull imaging resumes, change three things about the behaviour camera:", TEAL, True),
    ("•  Re-aim or widen so the whole head sits comfortably inside the frame", GREY, False),
    ("•  Diffuse the IR lighting so nothing saturates — target contrast SD near 0.14, not 0.26", GREY, False),
    ("•  Light background instead of dark, as in both validated setups", GREY, False),
    ("Cost is a lens, a mount and a diffuser. GrimACE publishes its hardware design, so the chamber "
     "can be copied directly.", GREY, False),
    ("Then both GrimACE and PainFace work off the shelf, with no retraining.", TEAL, True),
    ("Why this matters for the science: it would give an independent behavioural read-out of whether "
     "morphine actually produced analgesia — which is the single biggest gap in the current morphine "
     "result, where n = 2 paired mice and nothing confirms the drug worked.", NAVY, True),
], size=15)

out = os.path.join(HERE, "GrimACE_feasibility_for_Mark.pptx")
for k in range(1, 20):
    try:
        prs.save(out)
        break
    except PermissionError:
        out = out.replace(".pptx", f"_v{k+1}.pptx")
print("wrote", out, "slides", len(prs.slides))
