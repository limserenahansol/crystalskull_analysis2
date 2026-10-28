"""Six slides: one figure each, plus a summary and a conclusion."""
from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

FIG = Path(__file__).resolve().parent / "grimace_output" / "figures"
OUT = Path(__file__).resolve().parent / "GrimACE_crystal_skull_pain_assay.pptx"

NAVY = RGBColor(0x1F, 0x3A, 0x5F)
TEAL = RGBColor(0x0D, 0x5C, 0x63)
RED = RGBColor(0x9B, 0x1B, 0x30)
BLACK = RGBColor(0x22, 0x22, 0x22)

SLIDES = [
    {
        "title": "How GrimACE was applied to our face video",
        "image": "fig01_algorithm_steps.png",
        "summary": "A frame is cropped, resized to 224 px, and scored by the published ViT. Five action units each return 0, 1, 2, or not rateable. MGS is the sum of the rateable units (0–10).",
        "conclusion": "The published weights do run on our infrared faces. On this example the score is 1, and that 1 point is the ear alone.",
        "conclusion_color": TEAL,
    },
    {
        "title": "What 0, 1, and 3 points look like on our mice",
        "image": "fig08_what_one_point_means.png",
        "summary": "0 means the feature is absent, 1 means moderate, 2 means obvious. In 192 sampled frames from our cameras, no action unit was ever scored 2.",
        "conclusion": "An MGS of 1 here is not a pain face. The MGS 4 frame is a crop that missed the eyes, so its nose and cheek points are not a real face score.",
        "conclusion_color": RED,
    },
    {
        "title": "Representative faces: baseline, melittin, melittin + morphine",
        "image": "fig02_representative_faces.png",
        "summary": "Rows are mice. Columns are melittin-day baseline, melittin, morphine-day baseline, and melittin + morphine. Gray cells have no video for that condition.",
        "conclusion": "The camera is a tight front view, which is the kind of image GrimACE was built for. CS0204, mouse1f, and mouse2c have no morphine-day face video.",
        "conclusion_color": TEAL,
    },
    {
        "title": "Mouse-level score across the pain conditions",
        "image": "fig03_mouse_mgs_conditions.png",
        "summary": "Each point is one mouse: the mean of sampled frames with at least 3 rateable units. Lines connect the same mouse. Every mean is below 3 on a 0–10 scale.",
        "conclusion": "Melittin and melittin + morphine do not raise MGS in this sample. That is not evidence that pain was absent, or that morphine reduced pain.",
        "conclusion_color": RED,
    },
    {
        "title": "Pros and cons: which features actually carry the score",
        "image": "fig04_action_units.png",
        "summary": "Pro: five standard MGS units, automated, and already matched to EXTRACT. Con: the ear sits near 1 in every condition, nose and cheek stay near 0, and whisker is not rateable on about 60–75% of frames.",
        "conclusion": "Only orbital tightening moves, and it moves down after melittin + morphine. Three of the five pain features are not informative on this camera.",
        "conclusion_color": RED,
    },
    {
        "title": "Can GrimACE be the crystal-skull pain assay?",
        "image": "fig05_delta_inj_minus_bl.png",
        "summary": "Within each mouse, injection minus baseline is between −0.63 and +0.14 MGS points. The black bar is the mean of mice. n is mice, not frames.",
        "conclusion": "Not yet. Keep it as an exploratory face score. A trained person has to score the same crops before any pain claim. These shifts are smaller than 1 point on a 0–10 scale.",
        "conclusion_color": RED,
    },
]


def _add_text(slide, left, top, width, height, text, size, color, bold=False, label=""):
    box = slide.shapes.add_textbox(left, top, width, height)
    tf = box.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    if label:
        lab = p.add_run()
        lab.text = label
        lab.font.size = Pt(size)
        lab.font.bold = True
        lab.font.color.rgb = color
        lab.font.name = "Calibri"
    run = p.add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.color.rgb = color
    run.font.name = "Calibri"
    return box


def _fit(path: Path, max_w: int, max_h: int) -> tuple[int, int]:
    with Image.open(path) as im:
        w, h = im.size
    scale = min(max_w / w, max_h / h)
    return int(w * scale), int(h * scale)


def main() -> None:
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    blank = prs.slide_layouts[6]
    max_w, max_h = Inches(12.4), Inches(4.55)

    for spec in SLIDES:
        slide = prs.slides.add_slide(blank)
        _add_text(slide, Inches(0.4), Inches(0.18), Inches(12.5), Inches(0.42),
                  spec["title"], 26, NAVY, bold=True)

        img = FIG / spec["image"]
        if not img.is_file():
            raise FileNotFoundError(img)
        nw, nh = _fit(img, max_w, max_h)
        left = int((prs.slide_width - nw) / 2)
        slide.shapes.add_picture(str(img), Emu(left), Inches(0.68), Emu(nw), Emu(nh))

        _add_text(
            slide, Inches(0.45), Inches(5.40), Inches(12.4), Inches(0.85),
            spec["summary"], 15, BLACK, label="Summary.  ",
        )
        _add_text(
            slide, Inches(0.45), Inches(6.30), Inches(12.4), Inches(0.95),
            spec["conclusion"], 16, spec["conclusion_color"], label="Conclusion.  ",
        )

    prs.save(OUT)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
