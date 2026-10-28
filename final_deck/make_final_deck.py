"""Build the 5-slide summary deck for Mark, plus the supplement slides behind it.

Five slides carry the story. Everything a question would reach for - the exact test,
the n, the colour-scale definition, the file paths, the caveats - sits in the speaker
notes and in supplement slides 6-10, so the five stay readable.

python-pptx only. Never PowerPoint COM: it closes the documents Hansol has open.

    python make_final_deck.py
"""
from __future__ import annotations

import os

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
OUT = os.path.join(HERE, "CrystalSkull_summary_5slides_for_Mark.pptx")

W, H = Inches(13.333), Inches(7.5)
NAVY = RGBColor(0x1F, 0x3A, 0x5F)
INK = RGBColor(0x22, 0x22, 0x22)
GREY = RGBColor(0x55, 0x55, 0x55)
RED = RGBColor(0xC0, 0x39, 0x2B)

prs = Presentation()
prs.slide_width, prs.slide_height = W, H
BLANK = prs.slide_layouts[6]


def tb(slide, left, top, width, height, text, size=13, bold=False, color=INK,
       align=PP_ALIGN.LEFT, space_after=4, font="Calibri"):
    box = slide.shapes.add_textbox(left, top, width, height)
    tf = box.text_frame
    tf.word_wrap = True
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(space_after)
        run = p.add_run()
        run.text = line
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.color.rgb = color
        run.font.name = font
    return box


def figure_slide(title, png, caption, notes):
    """Title on top, one figure, one or two sentences of did/found at the bottom."""
    s = prs.slides.add_slide(BLANK)
    tb(s, Inches(.45), Inches(.22), Inches(12.45), Inches(.6), title,
       size=24, bold=True, color=NAVY)

    # 13.0 x 5.5 in figures -> 2.364 aspect. Fit to 4.85 in tall, centred.
    ph = Inches(4.72)
    pw = Emu(int(ph * (13.0 / 5.5)))
    s.shapes.add_picture(os.path.join(FIG, png), Emu(int((W - pw) / 2)),
                         Inches(.82), width=pw, height=ph)

    tb(s, Inches(.55), Inches(5.66), Inches(12.25), Inches(1.72), caption,
       size=12.5, color=INK, space_after=3)
    s.notes_slide.notes_text_frame.text = notes
    return s


def text_slide(title, blocks, notes=""):
    """A supplement slide: a heading and two columns of plain text."""
    s = prs.slides.add_slide(BLANK)
    tb(s, Inches(.45), Inches(.22), Inches(12.45), Inches(.6), title,
       size=22, bold=True, color=NAVY)
    left_w = Inches(6.1)
    for k, body in enumerate(blocks):
        tb(s, Inches(.52 + k * 6.35), Inches(1.00), left_w, Inches(6.1), body,
           size=11.5, color=INK, space_after=3)
    if notes:
        s.notes_slide.notes_text_frame.text = notes
    return s


# =====================================================================  1
figure_slide(
    "1.  All 30 recordings registered to the Allen CCF — no manual landmarks",
    "S1_registration.png",
    "What we did:  registered all 30 FOVs (18 face-poke + 12 melittin) to the Allen "
    "CCFv3 dorsal-surface map by maximising the Dice overlap between the tissue mask "
    "and the atlas cortex mask, after detecting the sagittal midline from the image "
    "itself.\n"
    "What we found:  all 30 pass (Dice 0.912–0.942, median 0.931) and the fitted scale "
    "agrees to within 1.05× across every FOV, so a cell's Allen area label is now a "
    "usable variable in the analysis.",
    "METHOD\n"
    "Atlas: BrainGlobe allen_mouse_25um_v1.2, 'asr' orientation, 528 x 320 x 456. "
    "Projected top-down; the first labelled voxel along each dorsal ray gives the "
    "surface area, and cortical layers are collapsed to their parent area.\n\n"
    "Registration is a similarity transform (rotation, isotropic scale, translation). "
    "Rotation comes from the auto-detected midline plus the fixed 90 deg between the "
    "FOV frame (A-P horizontal) and the atlas frame (A-P vertical). Scale and "
    "translation are found by maximising Dice. A-P sense is resolved by running both "
    "and keeping the better Dice.\n\n"
    "WHY NOT LANDMARKS\n"
    "I first asked Hansol to click bregma and lambda. That was my error - the skull is "
    "removed in this preparation, so neither is visible. The clicked registrations gave "
    "a 1.69x spread in fitted scale and an implied FOV width of 24.6 mm, which is "
    "impossible. The mask method needs no clicking and gives 1.05x.\n\n"
    "WHAT 'DICE' MEANS HERE\n"
    "2 x |A and B| / (|A| + |B|) between the binarised FOV tissue mask and the atlas "
    "cortex mask, after transform. 1.0 is perfect; 0.912 is the worst of the 30.\n\n"
    "VALIDATION\n"
    "Cerebellum mislabels dropped from 25% of posterior cells to 0%. Three known "
    "stereotaxic coordinates each land in the correct Allen area (see panel D).\n\n"
    "FILES\n"
    "Transforms: C:\\Users\\hsollim\\behavior_task\\crystal_skull_analysis\\"
    "ccf_registration\\transforms\\*.mat (30 files, each with dice, um_per_px, "
    "midline_angle_deg, ap_flipped).\n"
    "QC table: ...\\final_deck\\reg_qc.csv")

# =====================================================================  2
figure_slide(
    "2.  Face poke: somatosensory excitation, posterior inhibition",
    "S2_facepoke.png",
    "What we did:  classified each of 43,156 cells (18 recordings, 4 mice) as excited "
    "or inhibited by comparing the 0–2 s window after each poke with the 2 s before it, "
    "trial-paired within cell, then asked which Allen areas are over-represented.\n"
    "What we found:  8.6 % of cells are excited and 27.6 % inhibited; excitation "
    "concentrates in SSs (2.01×), SSp-n (1.42×) and MOp (1.37×), and inhibition is "
    "stronger in the posterior half of the FOV in 18 of 18 recordings (+7.2 pp, "
    "Wilcoxon signed-rank p = 7.6e-06).",
    "RESPONSE CLASSIFICATION (face poke)\n"
    "Per cell, per poke: mean dF/F in the response window 0 to +2 s minus mean dF/F in "
    "the base window -2 to 0 s. Amplitude = the mean of that difference across trials. "
    "Two-sided t across trials, uncorrected p < 0.05; sign of the amplitude gives "
    "excited vs inhibited. This is a within-cell, trial-paired design.\n\n"
    "Note that an uncorrected alpha of 0.05 would call about 2.5% of cells excited and "
    "2.5% inhibited by chance alone. 8.6% and 27.6% are both above that floor, but the "
    "excited fraction is only ~3x chance, so I would not lean on the absolute number - "
    "the regional comparison is the robust part, because the same test is applied "
    "everywhere.\n\n"
    "PANEL B - ENRICHMENT\n"
    "Fold enrichment = that area's excited-cell rate divided by the cortex-wide rate "
    "of 8.6%. 19 areas have at least 300 cells. Two-sided binomial test of each area "
    "against the cortex-wide rate, Benjamini-Hochberg across the 19. Stars: "
    "*** q<0.001, ** q<0.01, * q<0.05. Red = enriched, blue = depleted.\n\n"
    "PANEL C - A-P GRADIENT\n"
    "Each recording is split at its own median AP coordinate, so the test does not "
    "depend on how far forward that FOV happened to sit. Wilcoxon signed-rank over the "
    "18 recordings, W = 0 (every single recording moves the same way), p = 7.6e-06.\n\n"
    "The 4 animals are 0004749770_m1, m200, m201 and m3, derived from the movie name. "
    "The 'mouse' column in the recording index holds the FOV label (m1f / m2c) and "
    "would have collapsed these to 3.")

# =====================================================================  3
figure_slide(
    "3.  Morphine does not silence the melittin response (within animal)",
    "S3_morphine.png",
    "What we did:  compared melittin alone with melittin + morphine in the two animals "
    "(CS1014-1b, CS1014-2w) that were recorded under both arms, which removes the "
    "between-animal variance that dominates the pooled comparison.\n"
    "What we found:  in both animals the fraction of cells driven up rose (+2.8 and "
    "+2.6 pp) and mean amplitude rose (+2.4 and +2.3 % of baseline). With n = 2 pairs "
    "no p value is meaningful, so none is quoted — the pairing sharpens the estimate, "
    "it does not supply the evidence.",
    "WHY PAIRED\n"
    "Hansol confirmed CS1014-1b = 1014m01b and CS1014-2w = 1014m02w. Every earlier "
    "version of this comparison treated the two arms as independent groups (n = 4 vs "
    "2, unpaired Mann-Whitney). For these two animals that was the wrong design. The "
    "unpaired estimate was roughly 3x larger, inflated by one outlier animal "
    "(CS0204-2w); the within-animal effect is +2.7 pp on average.\n\n"
    "WHY NO P VALUE\n"
    "A two-sided sign test on n = 2 has a floor of p = 0.5, and a paired Wilcoxon on "
    "n = 2 cannot go below p = 0.5 either. Quoting any p here would be theatre. Two "
    "more paired animals would make this testable.\n\n"
    "RESPONSE CLASSIFICATION (melittin) - IMPORTANT CAVEAT\n"
    "Baseline and injection are separate EXTRACT runs, so cells are NOT tracked across "
    "the injection. Each injection cell's amplitude is its mean dF/F minus the grand "
    "mean over all baseline cells, tested with a t across time bins. Because there are "
    "thousands of time bins, that test has enormous degrees of freedom and calls a very "
    "large fraction of cells responsive - which is why panel A shows 28% up and 63% "
    "down, far above the face-poke numbers. Those absolute percentages are NOT "
    "comparable to slide 2 and should not be read as biology.\n\n"
    "What IS valid is the comparison between the two arms and between areas, because "
    "the identical test is applied to both sides. That is what panels B, C and D show.\n\n"
    "PANEL D\n"
    "Isocortex only, areas with at least 100 cells in all four sessions. Leaving every "
    "labelled area in pulls in SCzo (superior colliculus) and an 'outside' bin, which "
    "do not belong in a cortical claim.")

# =====================================================================  4
figure_slide(
    "4.  Face poke and melittin recruit different cortical territory",
    "S4_dissociation.png",
    "What we did:  put the two paradigms on the same axis by expressing each area's "
    "response rate as a fold enrichment over its own paradigm's cortex-wide rate, so "
    "the 3-fold difference in base rate cancels.\n"
    "What we found:  face poke is strongest in somatosensory cortex (SSs 2.01×, SSp-n "
    "1.42×) while melittin is strongest in posterior-medial parietal / visual areas "
    "(VISam 1.74×, VISpm 1.46×, VISa 1.42×); the two area maps are uncorrelated "
    "(Spearman ρ = −0.22, p = 0.36, 19 shared areas), and only MOp / MOs are enriched "
    "in both.",
    "WHY FOLD ENRICHMENT AND NOT RAW RATES\n"
    "The two paradigms use different response tests (see slide 3 notes), so their "
    "cortex-wide rates differ 3-fold: 8.6% excited for face poke, 28.3% driven up for "
    "melittin. Dividing each area by its own paradigm's cortex-wide rate removes that "
    "offset and leaves the spatial pattern, which is the thing we want to compare.\n\n"
    "STATISTICS\n"
    "Two-sided binomial test per area against the cortex-wide rate of that paradigm, "
    "Benjamini-Hochberg FDR across the 19 (face poke) and 27 (melittin) areas that have "
    "at least 300 cells. Spearman rank correlation across the 19 areas present in both "
    "lists.\n\n"
    "HOW TO READ THE SPEARMAN RESULT\n"
    "rho = -0.22 with p = 0.36 means the two maps are UNCORRELATED, not "
    "anti-correlated. The honest statement is that knowing an area's face-poke "
    "enrichment tells you nothing about its melittin enrichment. With 19 areas this "
    "test has limited power, so a modest true correlation would not be detected.\n\n"
    "SANITY CHECK\n"
    "SSs at 2.01x for face poke comes from 422 cells, which is the smallest n in the "
    "enriched set. It survives FDR (q = 7e-08), but it is the one number on this slide "
    "I would most want a second cohort to confirm.\n\n"
    "The melittin side uses the melittin-only arm (81,834 cells, 4 mice), not the "
    "morphine arm, so the comparison is drug-free on both sides.")

# =====================================================================  5
figure_slide(
    "5.  Face video is indexed and synced; grimace scoring needs a rig change",
    "S5_behavior.png",
    "What we did:  indexed all 19 facial-video sessions against the calcium recordings, "
    "then ran the published GrimACE model chain (frame quality → face cropper → ViT "
    "grimace scorer, authors' own weights and code) first on GrimACE's own example "
    "video as a control, then on ours.\n"
    "What we found:  it works on their footage (6/12 frames accepted) and finds no face "
    "at all in ours (0/12). Our frames are bimodal — dark IR background, bright mouse, "
    "SD 1.5× theirs — against a scorer that hard-codes an input mean of 0.5175 and SD "
    "of 0.14, so this is a rig mismatch, not a software problem.",
    "WHAT WAS ACTUALLY RUN\n"
    "The GitHub repo ships no weights (only models/.gitkeep). The trained weights and "
    "the labelled dataset are on Zenodo record 15119195: grimace_scorer.pt 341.5 MB, "
    "top_pose_estimator_l.pt 85.4 MB, top_pose_estimator.pt 22.8 MB, face_cropper.pt "
    "21.5 MB, frame_quality.pt 21.1 MB, 493 MB in total, plus 2,970 human-scored PNG "
    "frames and 15 score tables from 3 independent scorers. Those weights were dropped "
    "into the cloned repo's models/ folder and driven from a GUI-free runner "
    "(run_grimace_offline.py) so no camera is needed.\n\n"
    "THE CONTROL RUN MATTERS\n"
    "Running on their own example video first is what makes the negative result "
    "interpretable. The install is correct; the footage is the difference.\n\n"
    "WHY A GAIN CANNOT FIX IT\n"
    "The scorer normalises with hard-coded constants (x - 0.5175) / 0.14. A scalar gain "
    "moves the mean and the SD together, so it cannot match both. An earlier gain sweep "
    "I ran looked like evidence of scorer instability; it was arithmetic, not evidence, "
    "and I am not showing it.\n\n"
    "RECOMMENDATION\n"
    "1. light background behind the mouse  2. diffuse IR with no specular hot spot  "
    "3. wider framing so both ears stay in frame. With those three, either GrimACE or "
    "PainFace should run with no retraining.\n\n"
    "PAINFACE\n"
    "PainFace (Zylka lab, UNC) is the alternative and is a web service - it needs an "
    "account that only Hansol can create, at painface.net. Converted clips are already "
    "staged at ...\\behavior_sync\\painface_upload\\.\n\n"
    "SYNC WARNING\n"
    "Frame index x 1/35 s is NOT wall clock. 1014m01b 20260408 Injection holds 1,504 s "
    "of frames spread over 3,761 s of wall clock, so that recording is segmented. "
    "Always use the per-frame timestamps.")

# ================================================================  SUPPLEMENT
text_slide(
    "S1.  Methods — registration to the Allen CCF",
    ["ATLAS\n"
     "BrainGlobe allen_mouse_25um_v1.2 annotation volume, 'asr' orientation, "
     "528 × 320 × 456 voxels at 25 µm.\n"
     "A top-down ray through each dorsal column takes the first labelled voxel as the "
     "surface area. Cortical layers are collapsed to their parent area, so SSp-bfd2/3 "
     "and SSp-bfd5 both become SSp-bfd.\n"
     "Bregma sits at atlas row 216, column 228.\n"
     "\n"
     "TRANSFORM\n"
     "Similarity: rotation, one isotropic scale, translation. Four degrees of freedom.\n"
     "• rotation = auto-detected sagittal midline + the fixed 90° between the FOV frame "
     "(A–P horizontal) and the atlas frame (A–P vertical). Detected angles ran −4.0° to "
     "+7.5° over the 30 FOVs.\n"
     "• scale + translation = maximise Dice between the binarised FOV tissue mask and "
     "the atlas cortex mask.\n"
     "• A–P sense = run both, keep the better Dice.\n"
     "\n"
     "No cell coordinates enter the fit. The transform is estimated from the FOV "
     "silhouette alone, then applied to the EXTRACT footprint centroids.",

     "WHY NOT LANDMARK CLICKING\n"
     "The skull is removed in this preparation, so bregma and lambda cannot be seen. "
     "An initial clicked-landmark pass gave a 1.69× spread in fitted scale and an "
     "implied FOV width of 24.6 mm. The mask method gives 1.05× and needs no input.\n"
     "\n"
     "QUALITY, ALL 30 RECORDINGS\n"
     "Dice   median 0.931,  range 0.912 – 0.942\n"
     "Scale  median 5.48 µm/px,  range 5.38 – 5.65,  spread 1.05×\n"
     "Midline angle  −4.0° to +7.5°\n"
     "\n"
     "VALIDATION\n"
     "• cerebellum mislabels among posterior cells: 25 % → 0 %\n"
     "• 3.8 mm posterior / 2.5 lateral → VISp\n"
     "• 1.5 mm posterior / 3.2 lateral → SSp-bfd\n"
     "• 1.0 mm anterior / 1.5 lateral → MOp\n"
     "\n"
     "A known failure that was fixed: the first midline detector succeeded on only "
     "12/30 and reported confidence 1.00 everywhere, because the score was clipped at 1 "
     "and the connected-component step split one hemisphere. The projection-based "
     "detector with an unclipped score succeeds on 30/30."],
    "This supplement exists so slide 1 does not have to carry it.")

text_slide(
    "S2.  Methods — how a cell is called excited or inhibited",
    ["The two paradigms use DIFFERENT tests, so their responsive-cell percentages are "
     "NOT comparable to each other. Read this slide before quoting either number.\n"
     "\n"
     "FACE POKE  —  within-cell, trial-paired\n"
     "For each cell and each poke:\n"
     "   amplitude = mean dF/F over 0 to +2 s   −   mean dF/F over −2 to 0 s\n"
     "The per-cell amplitude is the mean of that difference over trials. Two-sided t "
     "across trials, uncorrected p < 0.05. Sign of the amplitude gives excited vs "
     "inhibited.\n"
     "Peri-stimulus window −2 to +3 s; only pokes whose full window lies inside the "
     "recording are used.\n"
     "\n"
     "Result: 8.6 % excited, 27.6 % inhibited, of 43,156 cells.\n"
     "\n"
     "MELITTIN  —  across-session, cells NOT tracked\n"
     "Baseline and injection are separate EXTRACT runs, so there is no cell "
     "correspondence across the injection.\n"
     "   amplitude = that injection cell's mean dF/F  −  the grand mean over all "
     "baseline cells\n"
     "tested with a t across time bins, uncorrected p < 0.05.\n"
     "\n"
     "Result: 28.3 % up, 63.0 % down (melittin only).",

     "WHY THE MELITTIN PERCENTAGES ARE SO MUCH HIGHER\n"
     "The melittin test uses time bins as its sample, so the degrees of freedom are in "
     "the thousands and almost any consistent offset reaches p < 0.05. The face-poke "
     "test uses trials, which is both the correct unit and a far smaller n.\n"
     "\n"
     "So:  28 % up / 63 % down and 8.6 % / 27.6 % are NOT the same measurement and "
     "must not be compared to each other as biology.\n"
     "\n"
     "WHAT REMAINS VALID\n"
     "Every comparison in this deck is made between conditions or between areas that "
     "were run through the SAME test, so the test's bias cancels:\n"
     "• slide 2, areas within face poke\n"
     "• slide 3, melittin vs melittin+morphine, same test both sides\n"
     "• slide 4, each paradigm normalised by its own cortex-wide rate\n"
     "\n"
     "WHAT WOULD FIX IT\n"
     "Register the baseline and injection EXTRACT runs to each other so cells can be "
     "tracked across the injection, then run the melittin arm as a within-cell paired "
     "test like the face-poke arm. That is the single biggest methodological upgrade "
     "available on this dataset."],
    "If Mark asks one hard question about this deck, it will be this one.")

text_slide(
    "S3.  Every statistic in the deck, in one place",
    ["SLIDE 1 — registration\n"
     "Dice, 30 recordings: median 0.931, range 0.912–0.942. Descriptive, no test.\n"
     "Scale: median 5.48 µm/px, 1.05× spread. Descriptive.\n"
     "\n"
     "SLIDE 2 — face poke,  43,156 cells / 18 recordings / 4 mice\n"
     "Excited 8.6 %, inhibited 27.6 % (per-cell t across trials, p < 0.05 uncorrected).\n"
     "Area enrichment: two-sided binomial vs the 8.6 % cortex-wide rate, BH-FDR over "
     "19 areas with ≥ 300 cells.\n"
     "   SSs 2.01×  q = 7e-08       SSp-n 1.42×  q = 4e-02\n"
     "   MOp 1.37×  q = 5e-16       MOs  1.29×  q = 3e-09\n"
     "   SSp-ul 0.66×  q = 4e-25    RSPagl 0.60×  q = 7e-08\n"
     "A–P gradient: inhibited fraction, posterior half vs anterior half, each recording "
     "split at its own median AP. Wilcoxon signed-rank, W = 0, p = 7.6e-06, n = 18; "
     "18/18 recordings move the same way; mean +7.2 pp.\n"
     "\n"
     "SLIDE 3 — melittin,  144,367 cells total\n"
     "melittin only: 81,834 cells, 4 mice, 28.3 % up / 63.0 % down\n"
     "melittin+morphine: 62,533 cells, 2 mice, 31.4 % up / 60.6 % down\n"
     "Paired, within animal (n = 2 pairs, NO p value quoted):\n"
     "   CS1014-1b  up 28.5 → 31.3 (+2.8 pp)   amp −5.7 → −3.3 (+2.4)\n"
     "   CS1014-2w  up 28.9 → 31.6 (+2.6 pp)   amp +0.5 → +2.7 (+2.3)\n"
     "By area: 11 of 22 isocortical areas (≥100 cells in all four sessions) agree in "
     "direction between the two animals — 50 %, i.e. chance.",

     "SLIDE 4 — paradigm comparison\n"
     "Enrichment computed separately per paradigm against its own cortex-wide rate "
     "(face poke 8.6 % excited; melittin 28.3 % driven up), binomial + BH-FDR over 19 "
     "and 27 areas respectively.\n"
     "Spearman rank correlation across the 19 shared areas: ρ = −0.22, p = 0.36. "
     "Uncorrelated, not anti-correlated.\n"
     "\n"
     "SLIDE 5 — behaviour\n"
     "19 sessions indexed, 15 with a matching calcium folder.\n"
     "GrimACE face cropper, identical weights and code and sampling:\n"
     "   their example video   6 / 12 frames accepted\n"
     "   our 1014m01b baseline  0 / 12 frames accepted\n"
     "Frame statistics: GrimACE rig mean 0.74 SD 0.17; our rig mean 0.41 SD 0.25 "
     "(SD 1.5× theirs). Scorer hard-codes mean 0.5175, SD 0.14.\n"
     "Descriptive; no inferential test is run on the grimace arm, and none should be "
     "until the cropper works.\n"
     "\n"
     "MULTIPLE COMPARISONS\n"
     "Benjamini-Hochberg is applied within each enrichment family (the 19 or 27 areas "
     "of one paradigm). It is NOT applied across slides — the per-cell p < 0.05 that "
     "defines a responsive cell is uncorrected, by design, and is discussed in S2."],
    "")

text_slide(
    "S4.  Data inventory and where everything lives",
    ["CALCIUM — face poke\n"
     "18 recordings, 4 animals (0004749770_m1, m200, m201, m3), 43,156 cells.\n"
     "Cell table:\n"
     "C:\\Users\\hsollim\\behavior_task\\crystal_skull_analysis\\ccf_registration\\"
     "n18\\cell_ccf_all.csv\n"
     "\n"
     "CALCIUM — melittin\n"
     "12 recordings, 6 mice, 144,367 cells.\n"
     "melittin only  CS0204-1b, CS0204-2w, CS1014-1b, CS1014-2w\n"
     "melittin+morphine  1014m01b, 1014m02w   (= CS1014-1b, CS1014-2w)\n"
     "Cell table:\n"
     "C:\\Users\\hsollim\\behavior_task\\crystal_skull_analysis\\ccf_registration\\"
     "melittin_ccf\\cell_ccf_all.csv\n"
     "\n"
     "REGISTRATION\n"
     "C:\\Users\\hsollim\\behavior_task\\crystal_skull_analysis\\ccf_registration\\"
     "transforms\\   (30 .mat + 30 QC .png)\n"
     "Atlas map:  ...\\ccf_registration\\allen_dorsal_map.mat\n"
     "QC table:   ...\\final_deck\\reg_qc.csv",

     "BEHAVIOUR VIDEO\n"
     "19 sessions on  \\\\stg-tnr50a1.stanford.edu\\crystal_skull\\\n"
     "   20260408  5 × .mp4\n"
     "   20260218, 20260219  14 × loose .jpg sequences\n"
     "Manifest:\n"
     "C:\\Users\\hsollim\\behavior_task\\crystal_skull_analysis\\behavior_sync\\"
     "behavior_calcium_manifest.csv\n"
     "PainFace staging:  ...\\behavior_sync\\painface_upload\\\n"
     "\n"
     "THIS DECK\n"
     "C:\\Users\\hsollim\\behavior_task\\crystal_skull_analysis\\final_deck\\\n"
     "   CrystalSkull_summary_5slides_for_Mark.pptx\n"
     "   make_final_figures.py, make_s5.py   (figures)\n"
     "   compute_stats.py                    (every number on the slides)\n"
     "   stats.json, reg_qc.csv              (the computed values)\n"
     "   figures\\S1..S5*.png\n"
     "\n"
     "EARLIER DECKS, STILL ON DISK AND UNCHANGED\n"
     "...\\mark_brief\\CrystalSkull_facepoke_vs_melittin_Mark_with Allen CCF_"
     "hansol_v2.pptx  (29 slides)\n"
     "...\\combined_paradigm\\CrystalSkull_facepoke_vs_melittin_vs_morphine_v2.pptx  "
     "(24 slides)\n"
     "...\\behavior_sync\\GrimACE_feasibility_for_Mark.pptx  (16 slides)\n"
     "Plus figures_v2\\fig1–fig17 and behavior_sync\\figures\\figB1–B8, figC1–C6."],
    "")

text_slide(
    "S5.  What each claim still needs",
    ["STRONG AS IT STANDS\n"
     "• The registration. 30/30, Dice 0.931, 1.05× scale spread, three independent "
     "coordinate checks. Nothing here depends on a judgement call.\n"
     "• The face-poke A–P gradient. 18/18 recordings, W = 0. This does not depend on "
     "the responsive-cell threshold, because the same threshold is applied to both "
     "halves of the same recording.\n"
     "• The face-poke somatosensory enrichment. Survives FDR comfortably for SSs, MOp "
     "and MOs.\n"
     "\n"
     "NEEDS MORE ANIMALS\n"
     "• Morphine. n = 2 paired animals. The direction is consistent and the effect is "
     "small (+2.7 pp). Two more paired animals would make this testable; at n = 4 a "
     "sign test can reach p = 0.125, at n = 6 it can reach p = 0.031.\n"
     "• SSp-n at 1.42× rests on 335 cells, SSs at 2.01× on 422. Both survive FDR but "
     "are the thinnest numbers in the deck.",

     "NEEDS A METHOD CHANGE\n"
     "• The melittin responsive-cell percentages. Cells are not tracked across the "
     "injection, so the test compares each injection cell to the baseline population "
     "mean using time bins as the sample. Registering the two EXTRACT runs to each "
     "other would let the melittin arm use the same within-cell paired test as face "
     "poke, and would make its percentages directly comparable to slide 2.\n"
     "• Grimace scoring. Blocked on the rig, not on the model: light background, "
     "diffuse IR, wider framing.\n"
     "\n"
     "OPEN AND WAITING ON SOMEONE ELSE\n"
     "• PainFace needs an account at painface.net that only Hansol can create. The "
     "converted clips are already staged.\n"
     "• Whether the 4 face-poke animals overlap with any melittin animal has not been "
     "confirmed. If they do, slide 4 could be made within-animal too.\n"
     "\n"
     "NOT CLAIMED ANYWHERE IN THIS DECK\n"
     "• That any of these cells are nociceptive specifically. The face-poke paradigm "
     "does not separate pain from touch or from surprise; that separation needs the "
     "behavioural arm working."],
    "")

prs.save(OUT)
print(f"wrote {OUT}")
print(f"{len(prs.slides.__iter__.__self__._sldIdLst)} slides")
