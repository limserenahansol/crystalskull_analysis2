"""Update Hansol's Mark deck from slide 19 onward.

Two kinds of change:
  1. result sentences that overstate what the data support, rewritten
  2. slide 23's table image swapped for the regenerated one, and new slides appended
     covering the cell-count and cell-resolution analyses

Everything before slide 19 is left untouched.
"""
import os

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.util import Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(os.path.dirname(HERE), "ccf_registration", "figures_v2")
SRC = os.path.join(HERE, "CrystalSkull_facepoke_vs_melittin_Mark_with Allen CCF_hansol.pptx")

NAVY = RGBColor(0x1F, 0x3A, 0x5F)
TEAL = RGBColor(0x0D, 0x5C, 0x63)
GREY = RGBColor(0x55, 0x55, 0x55)
RED = RGBColor(0x9B, 0x1B, 0x30)

prs = Presentation(SRC)


def caption_box(slide):
    """The 'What: / Result:' textbox, identified by position."""
    for sh in slide.shapes:
        if sh.has_text_frame and sh.text_frame.text.strip().startswith("What:"):
            return sh
    return None


def set_caption(slide, what, result, result_color=TEAL):
    box = caption_box(slide)
    tf = box.text_frame
    tf.clear()
    p0 = tf.paragraphs[0]
    r0 = p0.add_run()
    r0.text = "What:  " + what
    r0.font.size = Pt(12)
    r0.font.color.rgb = GREY
    p1 = tf.add_paragraph()
    r1 = p1.add_run()
    r1.text = "Result:  " + result
    r1.font.size = Pt(13)
    r1.font.bold = True
    r1.font.color.rgb = result_color


def swap_picture(slide, new_path):
    """Replace the slide's picture, keeping position and box, preserving aspect."""
    old = next(sh for sh in slide.shapes if sh.shape_type == 13)
    l, t, w, h = old.left, old.top, old.width, old.height
    old._element.getparent().remove(old._element)
    with Image.open(new_path) as im:
        iw, ih = im.size
    ar = iw / ih
    bw, bh = w, h
    nw, nh = bw, int(bw / ar)
    if nh > bh:
        nh, nw = bh, int(bh * ar)
    slide.shapes.add_picture(new_path, l + (bw - nw) // 2, t + (bh - nh) // 2, nw, nh)


# ===================================================== slide 19
# "morphine reversed suppression by melittin" states a demonstrated effect. The direction
# is consistent, but at the animal level it is not significant (n = 4 vs 2, no vehicle arm).
set_caption(
    prs.slides[18],
    "Local % of injection cells up (red) or down (blue) vs that day's baseline, on Allen CCF.",
    "Morphine shifts the melittin response toward less suppression - net population change "
    "-8.5% to -0.3%, consistent in direction across every measure. Not significant at the "
    "animal level (n = 4 vs 2 mice, no vehicle arm), so a trend, not an established effect.",
    RED,
)

# ===================================================== slide 20  (was already correct)
set_caption(
    prs.slides[19],
    "Local % of poke-excited (red) and poke-inhibited (blue) cells on Allen CCF. "
    "Same cells, peri-poke vs pre-poke.",
    "Poke responses are sparse and inhibition-dominated (8.6% excited vs 27.6% inhibited, "
    "43,156 cells). Excitation is biased anterior and somatosensory (SSs, SSp-n, MOp, MOs); "
    "inhibition is richer in visual cortex.",
)

# ===================================================== slide 23  new table image
swap_picture(prs.slides[22], os.path.join(FIG, "fig_which_regions_table.png"))
set_caption(
    prs.slides[22],
    "Recording-level Wilcoxon of (area rate minus that recording's overall rate), "
    "BH-FDR within each assay x direction.",
    "Face-poke excitation: MOp. Face-poke inhibition: VISl, VISam, VISp, VISrl. "
    "Melittin up: VISam, MOp, VISa, MOs, SSp-tr. Melittin+morphine cannot be tested at "
    "4 recordings, so its areas are ranked by enrichment instead.",
)

# ===================================================== slide 24
# The old sentence said 0/24 areas survive FDR. The test was never run for this arm: it
# needs 8 recordings and there are 4, and at n = 4 the signed-rank floor is p = 0.125, so
# no effect size could ever reach q < 0.05. "Failed FDR" and "could not be tested" are
# different claims and only the second one is true.
set_caption(
    prs.slides[23],
    "Recording-level % up (melittin+morphine) minus % up (melittin), by Allen area.",
    "Morphine raises % up in almost every area (+3 to +11 points), most in VISam, VISa, "
    "RSPd. With 4 recordings the recording-level test cannot be run at all - the signed-rank "
    "floor is p = 0.125 at n = 4 - so this is ranked by effect size, not tested.",
    RED,
)

# ===================================================== new slides
blank = prs.slide_layouts[6]


def add_fig_slide(title, what, result, img, result_color=TEAL):
    s = prs.slides.add_slide(blank)
    tb = s.shapes.add_textbox(Inches(0.35), Inches(0.12), Inches(12.6), Inches(0.38))
    r = tb.text_frame.paragraphs[0].add_run()
    r.text = title
    r.font.size = Pt(20)
    r.font.bold = True
    r.font.color.rgb = NAVY
    cb = s.shapes.add_textbox(Inches(0.35), Inches(0.48), Inches(12.6), Inches(0.58))
    cb.text_frame.word_wrap = True
    p0 = cb.text_frame.paragraphs[0]
    r0 = p0.add_run()
    r0.text = "What:  " + what
    r0.font.size = Pt(12)
    r0.font.color.rgb = GREY
    p1 = cb.text_frame.add_paragraph()
    r1 = p1.add_run()
    r1.text = "Result:  " + result
    r1.font.size = Pt(13)
    r1.font.bold = True
    r1.font.color.rgb = result_color
    with Image.open(img) as im:
        iw, ih = im.size
    ar = iw / ih
    bw, bh = Inches(12.4), Inches(6.2)
    nw, nh = bw, int(bw / ar)
    if nh > bh:
        nh, nw = bh, int(bh * ar)
    s.shapes.add_picture(img, int((prs.slide_width - nw) / 2), Inches(1.15), nw, nh)
    return s


add_fig_slide(
    "Face poke  |  exact cell counts per Allen area",
    "Bar length is the number of cells, so it reflects how much of each area the FOV "
    "sampled. The % alongside is the fair between-area comparison.",
    "SSp: 2,509 inhibited vs 837 excited of 10,432 cells. VIS: 5,329 vs 1,217 of 16,937. "
    "MO: 2,924 vs 1,309 of 11,378. Every area is inhibition-dominated in every mouse.",
    os.path.join(FIG, "fig9_counts_facepoke.png"),
)

add_fig_slide(
    "Face poke  |  cell counts painted on the brain",
    "Each Allen area filled by its own count, no smoothing. Labels on the left hemisphere "
    "only; counts are pooled over both.",
    "The excited-minus-inhibited map is negative in every single area - strongest in VISp "
    "(-2,155), MOp (-905), MOs (-710). No area has net positive drive by cell count.",
    os.path.join(FIG, "fig13_brain_counts_facepoke.png"),
)

add_fig_slide(
    "Face poke  |  median amplitude at cell resolution",
    "Every cell drawn at its own CCF position. A median cannot be drawn per cell, so the "
    "right brain panel shows each cell's sign - the median is negative exactly when most "
    "cells are.",
    "44% of cells are above zero, 56% below; median -0.025. Only SSs, SSp-n, SSp-m, SSp-ul "
    "and SSp-un are net positive - the somatosensory face and limb areas. The A-P gradient "
    "survives: 17/18 recordings, p = 1.5e-05.",
    os.path.join(FIG, "fig7b_facepoke_median_cellres.png"),
)

add_fig_slide(
    "Face poke  |  mean vs median - why the maps differ",
    "The same cells summarised two ways. Mean is the quantity in the original take-home "
    "map; median is what the typical cell does.",
    "Excited cells average +0.44 but are only 8.6%; inhibited cells average -0.19 and are "
    "27.6%. The 2.4x magnitude asymmetry makes the MEAN positive anteriorly while the "
    "MEDIAN stays negative. Report the median alongside the mean - 8 of 19 areas flip sign "
    "between them.",
    os.path.join(FIG, "fig7_facepoke_amplitude_ccf.png"),
    RED,
)

add_fig_slide(
    "Melittin  |  mean amplitude in Allen coordinates",
    "Percent change from each animal's own baseline session. Same quantity type as the "
    "face-poke amplitude map, different measurement.",
    "Melittin produces a cortex-wide DROP (-11.2%) that face poke does not, with an "
    "anterior > posterior gradient on top of it (8/8 recordings, p = 0.008) in the same "
    "direction as face poke. The gradient is not specific to face poke; the offset is "
    "specific to melittin.",
    os.path.join(FIG, "fig8_melittin_amplitude_ccf.png"),
)

out = os.path.join(HERE, "CrystalSkull_facepoke_vs_melittin_Mark_with Allen CCF_hansol_v2.pptx")
for attempt in range(1, 20):
    try:
        prs.save(out)
        break
    except PermissionError:
        out = out.replace(".pptx", f"_{attempt + 1}.pptx")
print("wrote", out, "slides", len(prs.slides))
