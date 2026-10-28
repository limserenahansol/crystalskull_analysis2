"""Combined three-condition deck: face poke vs melittin vs melittin + morphine.

Answers Mark's four questions in order and ends on the paradigm decision.
Figures come from combined_paradigm/ (make_combined_paradigm_figs.py) plus the
existing per-paradigm example-neuron panels.
"""
from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

NAVY = RGBColor(0x1F, 0x3A, 0x5F)
TEAL = RGBColor(0x0D, 0x5C, 0x63)
GREY = RGBColor(0x55, 0x55, 0x55)
RED = RGBColor(0x9B, 0x1B, 0x30)
AMBER = RGBColor(0xB0, 0x71, 0x0B)

ROOT = Path(r"C:\Users\hsollim\behavior_task\crystal_skull_analysis")
CMB = ROOT / "combined_paradigm"
POKE = ROOT / "stim_triggered_figures" / "combined_0813_0818"
MEL = ROOT / "melittin_lab_meeting"
OUT = CMB / "CrystalSkull_facepoke_vs_melittin_vs_morphine.pptx"


def blank(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])


def set_run(p, text, size, bold, color, align=PP_ALIGN.LEFT):
    p.text = text
    p.font.size = Pt(size)
    p.font.bold = bold
    p.font.color.rgb = color
    p.font.name = "Calibri"
    p.alignment = align


def add_textbox(slide, l, t, w, h, text, size, bold, color, align=PP_ALIGN.LEFT):
    box = slide.shapes.add_textbox(Inches(l), Inches(t), Inches(w), Inches(h))
    tf = box.text_frame
    tf.word_wrap = True
    set_run(tf.paragraphs[0], text, size, bold, color, align)
    return tf


def picture_fit(slide, img_path, left, top, max_w, max_h):
    img_path = Path(img_path)
    if not img_path.is_file():
        add_textbox(slide, left, top + max_h / 2, max_w, 0.5,
                    f"Missing: {img_path.name}", 14, True, RED, PP_ALIGN.CENTER)
        return
    with Image.open(img_path) as im:
        w, h = im.size
    ar = w / float(h)
    width, height = max_w, max_w / ar
    if height > max_h:
        height = max_h
        width = height * ar
    slide.shapes.add_picture(str(img_path), Inches(left + (max_w - width) / 2.0),
                             Inches(top + (max_h - height) / 2.0), Inches(width), Inches(height))


def figure_slide(prs, title, what, result, img, result_color=TEAL):
    slide = blank(prs)
    add_textbox(slide, 0.4, 0.12, 12.5, 0.40, title, 20, True, NAVY)
    tf = add_textbox(slide, 0.4, 0.50, 12.5, 0.36, "What:  " + what, 12, False, GREY)
    p = tf.add_paragraph()
    set_run(p, "Result:  " + result, 13, True, result_color)
    picture_fit(slide, img, 0.30, 1.05, 12.7, 6.20)
    return slide


def two_fig_slide(prs, title, what, result, img_l, img_r, cap_l="", cap_r=""):
    slide = blank(prs)
    add_textbox(slide, 0.30, 0.10, 12.7, 0.38, title, 18, True, NAVY)
    tf = add_textbox(slide, 0.30, 0.46, 12.7, 0.34, "What:  " + what, 12, False, GREY)
    p = tf.add_paragraph()
    set_run(p, "Result:  " + result, 13, True, TEAL)
    picture_fit(slide, img_l, 0.20, 1.05, 6.4, 5.7)
    picture_fit(slide, img_r, 6.70, 1.05, 6.4, 5.7)
    if cap_l:
        add_textbox(slide, 0.20, 6.85, 6.4, 0.4, cap_l, 11, True, TEAL, PP_ALIGN.CENTER)
    if cap_r:
        add_textbox(slide, 6.70, 6.85, 6.4, 0.4, cap_r, 11, True, AMBER, PP_ALIGN.CENTER)
    return slide


def bullets(prs, title, items, sizes=17, color=GREY):
    slide = blank(prs)
    add_textbox(slide, 0.6, 0.32, 12.1, 0.55, title, 24, True, NAVY)
    box = slide.shapes.add_textbox(Inches(0.85), Inches(1.10), Inches(11.7), Inches(6.0))
    tf = box.text_frame
    tf.word_wrap = True
    for i, item in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        if isinstance(item, tuple):
            text, col, bold = item
        else:
            text, col, bold = item, color, False
        set_run(p, text, sizes, bold, col)
        p.space_after = Pt(11)
    return slide


def section_slide(prs, title, subtitle):
    slide = blank(prs)
    add_textbox(slide, 0.7, 2.5, 11.9, 0.9, title, 32, True, TEAL, PP_ALIGN.CENTER)
    add_textbox(slide, 1.0, 3.6, 11.3, 1.4, subtitle, 17, False, GREY, PP_ALIGN.CENTER)
    return slide


def table_slide(prs, title, header, rows, colw):
    slide = blank(prs)
    add_textbox(slide, 0.5, 0.30, 12.3, 0.5, title, 22, True, NAVY)
    n_r, n_c = len(rows) + 1, len(header)
    shape = slide.shapes.add_table(n_r, n_c, Inches(0.55), Inches(1.15),
                                   Inches(sum(colw)), Inches(0.5 + 0.52 * len(rows)))
    tbl = shape.table
    for j, w in enumerate(colw):
        tbl.columns[j].width = Inches(w)
    for j, h in enumerate(header):
        c = tbl.cell(0, j)
        c.text = h
        pr = c.text_frame.paragraphs[0]
        pr.font.size = Pt(13)
        pr.font.bold = True
        pr.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
        c.fill.solid()
        c.fill.fore_color.rgb = NAVY
    for i, row in enumerate(rows, start=1):
        for j, v in enumerate(row):
            c = tbl.cell(i, j)
            c.text = str(v)
            pr = c.text_frame.paragraphs[0]
            pr.font.size = Pt(12)
            pr.font.color.rgb = GREY
            pr.font.bold = (j == 0)
    return slide


def find_example(base, prefer):
    per = base / "example_neurons" / "per_recording"
    if not per.is_dir():
        return None, None
    dirs = sorted(d for d in per.iterdir() if d.is_dir())
    if not dirs:
        return None, None
    pick = next((d for d in dirs if prefer in d.name), dirs[0])
    return pick / "example_traces.png", pick / "example_maps.png"


def main():
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)

    # ---------------- title
    s = blank(prs)
    add_textbox(s, 0.8, 2.0, 11.7, 1.0,
                "Crystal skull: face poke vs melittin vs melittin + morphine", 30, True, NAVY, PP_ALIGN.CENTER)
    add_textbox(s, 1.0, 3.1, 11.3, 0.9,
                "One combined analysis of both candidate pain paradigms, so we can pick one for the SBI drug work",
                17, False, GREY, PP_ALIGN.CENTER)
    add_textbox(s, 1.0, 4.3, 11.3, 1.6,
                "Face poke   18 recordings · 4 mice · 43,156 cells   (13 + 18 Aug 2025)\n"
                "Melittin only   8 recordings · 4 mice · 81,834 cells\n"
                "Melittin + morphine   4 recordings · 2 mice · 62,533 cells",
                15, False, TEAL, PP_ALIGN.CENTER)
    add_textbox(s, 1.0, 6.3, 11.3, 0.5, "Hansol Lim · Schnitzer lab · 9 Sep 2026", 13, False, GREY, PP_ALIGN.CENTER)

    # ---------------- the ask
    bullets(prs, "Why this deck", [
        ("Mark asked four things:", NAVY, True),
        "1.  How many cells per mouse were excited vs inhibited?",
        "2.  Which Allen brain areas do they sit in?",
        "3.  How will we separate pain- / touch- / surprise-responsive cells?",
        "4.  How were the excited/inhibited cells affected by drug?",
        ("Question 4 is the one that decides the project.", NAVY, True),
        "Face poke has no drug arm — those 18 recordings were a feasibility test of the paradigm itself. "
        "Melittin does have a drug arm. So to answer 4 at all, both datasets are put through the same "
        "pipeline and the same figures here.",
        ("The deck ends on a decision that you and Greg need to make.", RED, True),
    ], sizes=16)

    # ---------------- methods caveat
    bullets(prs, "Read this before comparing the numbers", [
        ("The two paradigms do not measure the same thing, and the percentages are not interchangeable.",
         RED, True),
        ("Face poke — stimulus-locked, within one session.", TEAL, True),
        "A cell counts as excited or inhibited if its peri-poke response differs from its own pre-poke "
        "baseline, across ~55 pokes. Same cells throughout. This is a phasic, trial-averaged measure.",
        ("Melittin — session-level, baseline movie vs injection movie.", AMBER, True),
        "A cell counts as up or down if its mean dF/F over the injection recording differs from the "
        "baseline-session mean. Baseline and injection are separate movies with independent EXTRACT IDs, "
        "so these are NOT tracked neurons. Both sessions are cropped to the shorter of the two.",
        ("Consequence: the melittin fractions are much larger (~27% / ~64%) than the face-poke fractions "
         "(~9% / ~28%) largely because a sustained session-level shift is an easier thing to detect than "
         "a phasic response. Compare the two paradigms qualitatively, not number against number.", NAVY, True),
    ], sizes=14)

    # ---------------- Q1
    section_slide(prs, "Q1  |  Excited vs inhibited cells, per mouse",
                  "The direction of the effect is consistent in every animal, in all three conditions")

    figure_slide(prs, "Fig 1  |  Excited vs inhibited across the three conditions",
                 "Per-recording fractions, mean +/- SD, each dot one recording. Right: every recording plotted.",
                 "All three conditions are inhibition-dominated. Face poke ~9% up / ~28% down; melittin "
                 "~27% / ~64%; adding morphine does not move the split.",
                 CMB / "fig1_fractions_three_paradigms.png")

    figure_slide(prs, "Fig 2  |  Cells excited vs inhibited, per mouse  (Mark's question 1)",
                 "Cells pooled within each animal, expressed as % of that animal's cells.",
                 "Every mouse shows the same direction. Face poke: inhibition ~3x more common than "
                 "excitation in 4/4 animals. Melittin: down ~2.3x up in 4/4. No animal reverses.",
                 CMB / "fig2_per_mouse.png")

    table_slide(prs, "Q1 numbers — face poke, per animal (18 recordings)",
                ["Animal (day)", "Rec", "Cells", "Excited", "Inhibited", "Pop. sig."],
                [["m1 (13 Aug)", 4, "9,705", "896  (9.2%)", "2,406  (24.8%)", "4/4"],
                 ["m201 (13 Aug)", 6, "15,467", "1,151  (7.4%)", "4,496  (29.1%)", "5/6"],
                 ["m3 (13 Aug)", 2, "3,685", "475  (12.9%)", "1,113  (30.2%)", "2/2"],
                 ["m200 (18 Aug)", 2, "5,681", "553  (9.7%)", "1,105  (19.5%)", "0/2"],
                 ["m201 (18 Aug)", 4, "8,618", "644  (7.5%)", "2,788  (32.4%)", "3/4"],
                 ["All", 18, "43,156", "3,719  (8.6%)", "11,908  (27.6%)", "14/18"]],
                [2.5, 1.0, 1.7, 2.5, 2.6, 1.6])

    # ---------------- Q2
    section_slide(prs, "Q2  |  Which Allen areas?",
                  "Honest answer: not yet — the FOVs are not registered to the Allen CCF")

    figure_slide(prs, "Q2  |  What we can say today: a consistent anterior-posterior gradient",
                 "Correlation between each cell's A-P position and its poke response, one value per recording, "
                 "mirror2 flipped to align with mirror1.",
                 "Excitation biased anterior, inhibition biased posterior, in 18/18 recordings "
                 "(median rho = 0.13, range 0.03-0.23). But 'anterior' here is the image axis after the "
                 "mirror flip, NOT a bregma-referenced coordinate — so no Allen area names yet.",
                 POKE / "n18" / "takehome_anterior_posterior.png", result_color=AMBER)

    bullets(prs, "Q2  |  What it would take to name Allen areas", [
        ("Today the crystal-skull FOVs have no CCF registration at all.", RED, True),
        "Landmark-based warp of each FOV to the Allen dorsal-cortex map: bregma and lambda from the "
        "skull image, plus the midline and the transverse sinus as anchors.",
        "That gives per-cell area labels (SSp-bfd, SSp-m, MOp, MOs, RSP, VISp ...) and lets us report "
        "excited/inhibited counts per area, per mouse, and compare across animals.",
        "It also links this dataset to the TRAP whole-brain and GPCR atlases, which are already in CCF space.",
        ("I would rather build the registration once, on whichever paradigm we commit to, than build it "
         "twice. That is one more reason the paradigm decision comes first.", NAVY, True),
    ], sizes=16)

    # ---------------- Q3
    section_slide(prs, "Q3  |  Pain vs touch vs surprise",
                  "This dataset cannot separate them, and no reanalysis will fix that")

    bullets(prs, "Q3  |  Why the current face-poke data cannot answer this", [
        ("There is one stimulus, delivered one way. There is no contrast to dissociate the categories.",
         RED, True),
        "A face poke is simultaneously a touch, a noxious event and an unexpected event. Any cell that "
        "responds is consistent with all three labels. Splitting them post hoc would be curve-fitting.",
        ("What would actually answer it — one session, interleaved trials:", NAVY, True),
        "•  innocuous brush or air puff        -> touch",
        "•  calibrated von Frey / noxious pinch or radiant heat   -> pain, graded intensity",
        "•  unexpected tone or light flash     -> surprise, non-somatosensory",
        "•  stimulus-omission catch trials     -> expectation / prediction error",
        "•  same cells across all trial types, so each cell gets a response profile rather than one number",
        ("Then a cell is 'pain-responsive' if it scales with noxious intensity and not with brush, rather "
         "than because we called the stimulus painful. This is a different experiment from what we have.",
         NAVY, True),
    ], sizes=15)

    # ---------------- Q4
    section_slide(prs, "Q4  |  Drug effect",
                  "Face poke: no drug arm at all.  Melittin: morphine arm exists, and it points one way — but n = 2 mice")

    figure_slide(prs, "Fig 3  |  Melittin only vs melittin + morphine",
                 "Change in population mean dF/F, injection vs that animal's own baseline, as % of baseline. "
                 "Left: per recording. Right: per mouse.",
                 "Melittin alone lowers population activity in 7/8 recordings (mean -8.5%). With morphine "
                 "on board the shift is ~0 (-0.3%). But n = 4 vs 2 mice and it is not significant: "
                 "p = 0.15 per recording, p = 0.27 per mouse.",
                 CMB / "fig3_drug_effect.png", result_color=AMBER)

    figure_slide(prs, "Fig 4  |  Which way does morphine move the melittin response?",
                 "All four measures, melittin-only vs melittin + morphine, tested per recording and again "
                 "per animal (the two mirrors of a session are not independent samples).",
                 "Direction is consistent: every measure that moves, moves toward MORE activity with "
                 "morphine. Melittin alone suppresses the cortex; with morphine that suppression is gone. "
                 "But only fraction-UP reaches p < 0.05, and only at the recording level (p = 0.048); "
                 "per animal it is p = 0.13 and nothing else is close.",
                 CMB / "fig4_morphine_direction.png", result_color=AMBER)

    bullets(prs, "Q4  |  What we can and cannot say about morphine", [
        ("Can say — melittin is a strong, repeatable pain read-out.", TEAL, True),
        "~27% of cells up and ~64% down relative to baseline, in 4/4 animals, with a net suppression of "
        "population activity of -8.5% in 7/8 recordings.",
        ("Can say — the direction of the morphine effect is consistent.", TEAL, True),
        "With morphine on board there are more cells above baseline (27.2% -> 31.4%) and the net "
        "suppression is abolished (-8.5% -> -0.3%). Read physiologically: morphine reverses the "
        "melittin-induced cortical suppression. Nothing moves the other way.",
        ("Cannot say — that this is an established effect.", RED, True),
        "Only fraction-UP crosses p < 0.05, and only if the two mirrors of one session are counted as "
        "independent samples. At the animal level (n = 4 vs 2) nothing is significant: fraction UP "
        "p = 0.13, net population p = 0.27, UP amplitude p = 0.53.",
        ("Cannot say — anything about 'the neurons'.", RED, True),
        "Baseline and injection are separate movies with independent EXTRACT IDs, so no cell is followed "
        "across the drug. Every statement here is about populations, not identified neurons.",
        ("Fix that is available today: if 1014m01b / 1014m02w are the same animals as CS1014-1b / "
         "CS1014-2w on a later day, this becomes a paired within-animal test on the same data. I need "
         "someone to confirm the animal identities.", NAVY, True),
    ], sizes=13)

    # ---------------- examples
    section_slide(prs, "Same visual language, both paradigms",
                  "Example single-cell traces and the maps of where those cells sit in the FOV")

    p_tr, p_mp = find_example(POKE, "m1f_s00")
    m_tr, m_mp = find_example(MEL, "CS1014-1b")
    if p_tr:
        two_fig_slide(prs, "Face poke — example cells and their position in the FOV",
                      "Two excited and two inhibited cells per movie, EXTRACT activity with poke times, "
                      "and the same cells marked on the mesoscope FOV (P left, A right).",
                      "Poke lock is visible cell by cell: excited cells rise after each poke, inhibited "
                      "cells dip. E1/I1 more anterior, E2/I2 more posterior.",
                      p_tr, p_mp, "Traces", "FOV map")
    if m_tr:
        two_fig_slide(prs, "Melittin — example cells and their position in the FOV",
                      "Same panel design, but the comparison is injection vs baseline session rather than "
                      "peri-poke, and the cell IDs are injection-session IDs.",
                      "Up and down cells are intermingled across the FOV, as with face poke.",
                      m_tr, m_mp, "Traces", "FOV map")

    two_fig_slide(prs, "Population view, side by side",
                  "Left: face-poke peri-stimulus population response, pooled over 18 recordings. "
                  "Right: melittin baseline vs injection population dF/F.",
                  "Face poke gives a brief transient, a short dip, then a rebound, time-locked to the poke. "
                  "Melittin gives a slow sustained shift with no trial structure.",
                  POKE / "pooled_average" / "fig1_population_psth.png",
                  MEL / "cohort_average" / "fig1_population_psth.png",
                  "Face poke (trial-locked)", "Melittin (session-level)")

    # ---------------- decision
    section_slide(prs, "The decision", "Which paradigm do we build the SBI drug plan around?")

    table_slide(prs, "Face poke vs melittin — what each buys us",
                ["", "Face poke", "Melittin"],
                [["Pain type", "Phasic, mixed (touch + nox + surprise)", "Sustained inflammatory-type"],
                 ["Trials", "~55 pokes per session, trial-locked", "One injection, no trial structure"],
                 ["Cell tracking", "Same cells throughout session", "Independent IDs BL vs INJ"],
                 ["Statistics", "Clean per-cell trial statistics", "Session-level only"],
                 ["Drug arm", "None — no drug data at all", "Morphine exists (n = 2 mice, n.s.)"],
                 ["Throughput", "Fast, repeatable, same animal over days", "Slow, one shot per session"],
                 ["Best fit for", "Circuit mapping, stimulus panel designs", "Analgesic efficacy over time"]],
                [2.3, 4.9, 4.9])

    bullets(prs, "What I need from Mark and Greg", [
        ("Pick the paradigm. These are two different experimental plans and we should commit to one.",
         NAVY, True),
        "If melittin: the immediate job is a proper drug design — vehicle arm, within-animal morphine, "
        "enough mice to power it. The current morphine data is suggestive at best and I would not build "
        "on it as it stands.",
        "If face poke: the immediate job is the stimulus panel from Q3, so pain / touch / surprise is "
        "answerable by design, plus a drug arm which currently does not exist.",
        ("Either way I will add the Allen CCF registration once, on the chosen paradigm.", TEAL, True),
        ("I have paused crystal-skull imaging until this is settled — I would rather hold than collect "
         "another 18 recordings on the paradigm we end up not using.", RED, True),
    ], sizes=16)

    # ---------------- take home
    s = blank(prs)
    add_textbox(s, 0.7, 0.45, 11.9, 0.6, "Take home", 28, True, NAVY, PP_ALIGN.CENTER)
    add_textbox(s, 1.0, 1.5, 11.3, 4.8,
                "1.  Every mouse, in every condition, is inhibition-dominated. Face poke 8.6% excited vs "
                "27.6% inhibited across 43,156 cells; melittin 27% up vs 64% down. No animal reverses.\n\n"
                "2.  No Allen area names yet — the FOVs are not CCF-registered. What is solid is the "
                "anterior-excited / posterior-inhibited gradient, 18/18 recordings.\n\n"
                "3.  Pain vs touch vs surprise cannot be separated in this dataset. It needs an "
                "interleaved stimulus panel, not more analysis.\n\n"
                "4.  Morphine moves the melittin response toward more activity, not less: more cells "
                "above baseline (27.2% -> 31.4%) and the net cortical suppression abolished "
                "(-8.5% -> -0.3%). The direction is consistent across every measure, but at the animal "
                "level (n = 4 vs 2, no vehicle arm) it is not significant. It is a trend worth designing "
                "for, not a result to build on.\n\n"
                "5.  So the drug question is still open on both paradigms, and the decision on which one "
                "to use is the blocking step.",
                17, False, GREY)

    out = OUT
    for attempt in range(1, 20):
        try:
            prs.save(str(out))
            break
        except PermissionError:
            # file is open in PowerPoint; fall back to a versioned name
            out = OUT.with_name(f"{OUT.stem}_v{attempt + 1}{OUT.suffix}")
    else:
        raise RuntimeError("could not write the deck under any name")
    print("Wrote", out, "slides", len(prs.slides))


if __name__ == "__main__":
    main()
