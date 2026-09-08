"""Lab-meeting PPT: method + 1-sentence 'what' + 1-sentence result per figure.

Builds both peak_idx and peak_pre_idx decks.

Run:  python build_crystal_skull_lab_pptx.py
"""

from pathlib import Path
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from PIL import Image

NAVY = RGBColor(0x1F, 0x3A, 0x5F)
TEAL = RGBColor(0x0D, 0x5C, 0x63)
RED = RGBColor(0x9B, 0x1B, 0x30)
GREY = RGBColor(0x55, 0x55, 0x55)

RECS = [
    ("m1f_s00_pwr8.38_mirror1", "mouse1f  s00  8.38 mW/cm2  mirror1"),
    ("m1f_s00_pwr8.38_mirror2", "mouse1f  s00  8.38 mW/cm2  mirror2"),
    ("m1f_s01_pwr20_mirror1", "mouse1f  s01  20 mW/cm2  mirror1"),
    ("m1f_s01_pwr20_mirror2", "mouse1f  s01  20 mW/cm2  mirror2"),
    ("m2c_s00_pwr8.38_mirror1", "mouse2c  s00  8.38 mW/cm2  mirror1"),
    ("m2c_s00_pwr8.38_mirror2", "mouse2c  s00  8.38 mW/cm2  mirror2"),
    ("m2c_s01_pwr20_mirror1", "mouse2c  s01  20 mW/cm2  mirror1"),
    ("m2c_s01_pwr20_mirror2", "mouse2c  s01  20 mW/cm2  mirror2"),
    ("m2c_s02_pwrNA_mirror1", "mouse2c  s02  power TBD  mirror1"),
    ("m2c_s02_pwrNA_mirror2", "mouse2c  s02  power TBD  mirror2"),
    ("m3_m1f_s00_pwr8.38_mirror1", "m3 (LabVIEW m1f)  s00  mirror1"),
    ("m3_m1f_s00_pwr8.38_mirror2", "m3 (LabVIEW m1f)  s00  mirror2"),
]

DECKS = {
    "pre": {
        "lm": Path(r"C:\Users\hsollim\behavior_task\crystal_skull_analysis\stim_triggered_figures\20250813_lab_meeting_pre\lab_meeting"),
        "out_name": "CrystalSkull_face_poke_labmeeting_pre.pptx",
        "event": "pre",
        "t0_short": "peak_pre_idx (onset candidate)",
        "title_sub": "2025-08-13   |   t = 0 at peak_pre_idx (onset candidate)   |   window 0.4–1.2 s vs -2–0 s",
        "takehome": "Take-home: activity dips at the pre-minimum then rises ~1.5 s later; red/blue maps show where those cells sit in each FOV.",
        "sync_body": "No behavior camera. A blocker on the arm changes mean-frame intensity. t = 0 is peak_pre_idx: the intensity minimum just before the positive brightness peak (candidate stimulus onset).",
        "next_body": "Compare this onset-candidate map with the peak_idx deck (stronger lock, 0.4–1.2 s window was chosen for that peak). FOVs are not aligned across mice or mirrors.",
        "figs": [
            {
                "file": "fig1_population_psth.png",
                "title": "Fig1  Population PSTH",
                "what": "Mean population activity (z-score) and fraction of active cells, aligned to each poke (t = 0 is peak_pre_idx).",
                "cohort_result": "Activity dips at poke onset, then rises and peaks about 1.5 s later — a delayed, poke-locked population response.",
            },
            {
                "file": "fig2_neuron_sta.png",
                "title": "Fig2  Neuron STA",
                "what": "Each neuron’s average trace around the poke, sorted by response size; right panel is the mean of excited vs inhibited cells.",
                "cohort_result": "A subset of cells go up after the poke and another subset go down; the all-cell mean is small because they partly cancel.",
            },
            {
                "file": "fig3_responsive_cells.png",
                "title": "Fig3  Excited vs inhibited cells",
                "what": "Per neuron: mean z in 0.4–1.2 s after t=0 vs -2–0 s before; signed-rank, p<0.05, sign of the difference.",
                "cohort_result": "About 5–13% of cells are excited and 17–24% inhibited; only 5 of 12 movies pass the population test with this window.",
            },
            {
                "file": "fig4_spatial_response.png",
                "title": "Fig4  Spatial poke-response map",
                "what": "Each neuron’s post−pre amplitude is painted onto its EXTRACT footprint. Red = excited, blue = inhibited. Left: anatomy. Middle: all cells. Right: p<0.05 only.",
                "cohort_result": "Excited and inhibited cells are intermingled across dorsal cortex; there is no single on-or-off region in these unaligned FOVs.",
            },
            {
                "file": "fig5_trial_timing.png",
                "title": "Fig5  Trial raster and poke interval",
                "what": "One row per poke (population z over time) and the histogram of time between detected pokes.",
                "cohort_result": "About 48–57 pokes per movie; median interval is 5.3 s (set to 3 s, with mechanical jitter). Long gaps are missed peaks.",
            },
        ],
    },
    "peak": {
        "lm": Path(r"C:\Users\hsollim\behavior_task\crystal_skull_analysis\stim_triggered_figures\20250813_lab_meeting\lab_meeting"),
        "out_name": "CrystalSkull_face_poke_labmeeting_peak.pptx",
        "event": "peak",
        "t0_short": "peak_idx (brightness peak)",
        "title_sub": "2025-08-13   |   t = 0 at peak_idx (positive brightness peak)   |   window 0.4–1.2 s vs -2–0 s",
        "takehome": "Take-home: cortex is tightly time-locked to the brightness peak (11/12 movies); red/blue maps show excited vs inhibited cells on each FOV.",
        "sync_body": "No behavior camera. A blocker on the arm changes mean-frame intensity. t = 0 is peak_idx: the positive brightness peak. This is the strongest neural alignment, but it is not proven to be the moment of skin contact.",
        "next_body": "Compare these maps with the peak_pre_idx deck (better physical candidate for skin-contact onset). FOVs are not aligned across mice or mirrors.",
        "figs": [
            {
                "file": "fig1_population_psth.png",
                "title": "Fig1  Population PSTH",
                "what": "Mean population activity (z-score) and fraction of active cells, aligned to each poke (t = 0 is peak_idx).",
                "cohort_result": "A sharp transient at t = 0, a brief dip, then a rebound peaking ~0.8 s — the same shape in all 12 movies.",
            },
            {
                "file": "fig2_neuron_sta.png",
                "title": "Fig2  Neuron STA",
                "what": "Each neuron’s average trace around the poke, sorted by response size; right panel is the mean of excited vs inhibited cells.",
                "cohort_result": "Excited cells rise after the peak; inhibited cells go down; the all-cell mean is small because the two groups partly cancel.",
            },
            {
                "file": "fig3_responsive_cells.png",
                "title": "Fig3  Excited vs inhibited cells",
                "what": "Per neuron: mean z in 0.4–1.2 s after t=0 vs -2–0 s before; signed-rank, p<0.05, sign of the difference.",
                "cohort_result": "About 6–14% of cells are excited and 23–31% inhibited; 11 of 12 movies have a significant population response.",
            },
            {
                "file": "fig4_spatial_response.png",
                "title": "Fig4  Spatial poke-response map",
                "what": "Each neuron’s post−pre amplitude is painted onto its EXTRACT footprint. Red = excited, blue = inhibited. Left: anatomy. Middle: all cells. Right: p<0.05 only.",
                "cohort_result": "Excited and inhibited cells are intermingled across dorsal cortex; there is no single on-or-off region in these unaligned FOVs.",
            },
            {
                "file": "fig5_trial_timing.png",
                "title": "Fig5  Trial raster and poke interval",
                "what": "One row per poke (population z over time) and the histogram of time between detected pokes.",
                "cohort_result": "About 48–57 pokes per movie; median interval is 5.3 s (set to 3 s, with mechanical jitter). Long gaps are missed peaks.",
            },
        ],
    },
    "post": {
        "lm": Path(r"C:\Users\hsollim\behavior_task\crystal_skull_analysis\stim_triggered_figures\20250813_lab_meeting_post\lab_meeting"),
        "out_name": "CrystalSkull_face_poke_labmeeting_post.pptx",
        "event": "post",
        "t0_short": "peak_post_idx (after the brightness peak)",
        "title_sub": "2025-08-13   |   t = 0 at peak_post_idx (minimum after the brightness peak)   |   window 0.4–1.2 s vs -2–0 s",
        "takehome": "Take-home: the neural response happens before this t = 0, so post is a poor onset; the maps still show where significant cells sit.",
        "sync_body": "No behavior camera. A blocker on the arm changes mean-frame intensity. t = 0 is peak_post_idx: the intensity minimum just after the positive brightness peak (candidate offset of the light change).",
        "next_body": "Compare pre (onset candidate) and peak (strongest neural lock) decks; post is mainly a negative control for timing. FOVs are not aligned across mice or mirrors.",
        "figs": [
            {
                "file": "fig1_population_psth.png",
                "title": "Fig1  Population PSTH",
                "what": "Mean population activity (z-score) and fraction of active cells, aligned to each poke (t = 0 is peak_post_idx).",
                "cohort_result": "The main rise is before t = 0 (~0.75 s earlier, at the brightness peak); at t = 0 activity dips, then slowly returns — this marker is after the neural event.",
            },
            {
                "file": "fig2_neuron_sta.png",
                "title": "Fig2  Neuron STA",
                "what": "Each neuron’s average trace around the poke, sorted by response size; right panel is the mean of excited vs inhibited cells.",
                "cohort_result": "Excited and inhibited groups still split, but their peaks sit before t = 0, so this alignment smears the response into the baseline window.",
            },
            {
                "file": "fig3_responsive_cells.png",
                "title": "Fig3  Excited vs inhibited cells",
                "what": "Per neuron: mean z in 0.4–1.2 s after t=0 vs -2–0 s before; signed-rank, p<0.05, sign of the difference.",
                "cohort_result": "Only 2 of 12 movies pass the population test; excited fraction is low (~2–6%) because the real response is mostly before t = 0.",
            },
            {
                "file": "fig4_spatial_response.png",
                "title": "Fig4  Spatial poke-response map",
                "what": "Each neuron’s post−pre amplitude is painted onto its EXTRACT footprint. Red = excited, blue = inhibited. Left: anatomy. Middle: all cells. Right: p<0.05 only.",
                "cohort_result": "Excited and inhibited cells are intermingled across dorsal cortex; there is no single on-or-off region in these unaligned FOVs.",
            },
            {
                "file": "fig5_trial_timing.png",
                "title": "Fig5  Trial raster and poke interval",
                "what": "One row per poke (population z over time) and the histogram of time between detected pokes.",
                "cohort_result": "About 48–57 pokes per movie; median interval is 5.3 s (set to 3 s, with mechanical jitter). Long gaps are missed peaks.",
            },
        ],
    },
}


def fmt_p(p):
    p = float(p)
    if p < 0.001:
        return f"p = {p:.1e}"
    return f"p = {p:.3f}"


def load_index(lm):
    rows = {}
    lines = (lm / "recording_index.csv").read_text(encoding="utf-8").strip().splitlines()
    hdr = lines[0].split(",")
    for line in lines[1:]:
        d = dict(zip(hdr, line.split(",")))
        rows[d["label"]] = d
    return rows


def rec_results(d, fig_file):
    n = int(d["n_neurons"])
    pk = int(d["n_pokes"])
    fe = 100 * float(d["frac_exc"])
    fi = 100 * float(d["frac_inh"])
    p = fmt_p(d["p_pop"])
    sig = "significant" if float(d["p_pop"]) < 0.05 else "not significant"
    if fig_file.endswith("psth.png"):
        return f"This FOV has {n} neurons and {pk} pokes; the population pre vs post test is {sig} ({p})."
    if fig_file.endswith("sta.png"):
        return f"Excited {fe:.1f}% and inhibited {fi:.1f}% of {n} neurons; see the STA split on the right."
    if fig_file.endswith("cells.png"):
        return f"{fe:.1f}% excited, {fi:.1f}% inhibited out of {n} neurons; population test {sig} ({p})."
    if "spatial_response" in fig_file or "NOT_response" in fig_file or "anatomy" in fig_file:
        return f"Red = excited, blue = inhibited on this FOV ({fe:.1f}% exc, {fi:.1f}% inh of {n} neurons)."
    if "timing" in fig_file:
        return f"{pk} pokes detected in this movie; intervals near 5.3 s, with occasional missed peaks."
    return ""


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


def two_line_caption(slide, what, result):
    tf = add_textbox(slide, 0.4, 0.72, 12.5, 0.42, "What:  " + what, 14, False, GREY)
    p = tf.add_paragraph()
    set_run(p, "Result:  " + result, 15, True, TEAL)


def picture_fit(slide, img_path, left, top, max_w, max_h):
    with Image.open(img_path) as im:
        w, h = im.size
    ar = w / float(h)
    width, height = max_w, max_w / ar
    if height > max_h:
        height = max_h
        width = height * ar
    slide.shapes.add_picture(
        str(img_path),
        Inches(left + (max_w - width) / 2.0),
        Inches(top + (max_h - height) / 2.0),
        Inches(width),
        Inches(height),
    )


def method_slide(prs, title, body):
    slide = blank(prs)
    add_textbox(slide, 0.7, 2.1, 11.9, 0.7, title, 22, True, NAVY, PP_ALIGN.CENTER)
    add_textbox(slide, 1.2, 3.0, 10.9, 2.2, body, 22, False, GREY, PP_ALIGN.CENTER)


def section_slide(prs, title, subtitle):
    slide = blank(prs)
    add_textbox(slide, 0.7, 2.4, 11.9, 0.9, title, 32, True, TEAL, PP_ALIGN.CENTER)
    add_textbox(slide, 1.0, 3.5, 11.3, 1.2, subtitle, 18, False, GREY, PP_ALIGN.CENTER)


def figure_slide(prs, title, what, result, img_path):
    slide = blank(prs)
    add_textbox(slide, 0.4, 0.16, 12.5, 0.5, title, 24, True, NAVY)
    two_line_caption(slide, what, result)
    if img_path is not None and Path(img_path).is_file():
        picture_fit(slide, img_path, 0.3, 1.55, 12.7, 5.7)
    else:
        add_textbox(slide, 1, 3.5, 11, 0.6, "Figure file missing", 20, True, RED, PP_ALIGN.CENTER)


BY_EVENT = Path(
    r"C:\Users\hsollim\behavior_task\crystal_skull_analysis\stim_triggered_figures\20250813_spatial_by_event"
)
WORK = Path(
    r"C:\Users\hsollim\behavior_task\crystal_skull_analysis\stim_triggered_figures\20250813_lab_meeting_pre\lab_meeting"
)

FIG4_WHAT = (
    "Each neuron's post-pre amplitude is painted onto its EXTRACT footprint. "
    "Red = excited, blue = inhibited. Left: anatomy. Middle: all cells. Right: p<0.05 only."
)
FIG4_RESULT_COHORT = (
    "Excited and inhibited cells are intermingled across dorsal cortex; "
    "there is no single on-or-off region in these unaligned FOVs."
)
BY_EVENT_WHAT = (
    "Same neurons and FOV, three t=0 definitions: pre (onset candidate), "
    "peak (brightness peak), post (after the peak). Shared color scale."
)
BY_EVENT_RESULT = (
    "The spatial pattern is similar across events: mixed red and blue cells. "
    "Peak shows the strongest lock; post is weaker because the response is mostly before t=0."
)


def save_pptx(prs, path):
    path = Path(path)
    try:
        prs.save(str(path))
        print("Wrote", path, "slides", len(prs.slides))
        return path
    except PermissionError:
        alt = path.with_name(path.stem + "_updated.pptx")
        prs.save(str(alt))
        print("Locked, wrote", alt, "slides", len(prs.slides))
        return alt


def copy_if_possible(src, dst):
    src = Path(src)
    dst = Path(dst)
    if src.resolve() == dst.resolve():
        return
    try:
        dst.write_bytes(src.read_bytes())
        print("Copied", dst)
    except PermissionError:
        alt = dst.with_name(dst.stem + "_updated.pptx")
        alt.write_bytes(src.read_bytes())
        print("Locked, copied", alt)


def build_deck(cfg):
    lm = cfg["lm"]
    per = lm / "per_recording"
    avg = lm / "cohort_average"
    out = lm / cfg["out_name"]
    idx = load_index(lm)
    figs = cfg["figs"]

    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)

    slide = blank(prs)
    add_textbox(slide, 0.8, 1.8, 11.7, 1.2, "Face poking on crystal skull", 36, True, NAVY, PP_ALIGN.CENTER)
    add_textbox(slide, 1.2, 3.2, 10.9, 1.4, cfg["title_sub"], 18, False, GREY, PP_ALIGN.CENTER)
    add_textbox(slide, 1.2, 4.6, 10.9, 1.2, cfg["takehome"], 18, True, TEAL, PP_ALIGN.CENTER)

    method_slide(
        prs, "Experiment",
        "Whole dorsal cortex Ca2+ imaging while a motor arm pokes the face. Every movie on this day is a poking session — there is no no-poke baseline.",
    )
    method_slide(
        prs, "Two recordings of the same experiment",
        "s00 = 8.38 mW/cm2, s01 = 20 mW/cm2. mirror1 and mirror2 are two optical paths (two FOVs at the same time).",
    )
    method_slide(prs, "How poke time is defined", cfg["sync_body"])
    method_slide(
        prs, "How traces are analyzed",
        "EXTRACT traces at 33 Hz. Cut -2 to +3 s around each poke. A cell is excited or inhibited if 0.4–1.2 s differs from -2–0 s (Wilcoxon signed-rank, p<0.05).",
    )
    method_slide(
        prs, "How the spatial map is made",
        "EXTRACT spatial footprints were reconstructed from the file (ndSparse CSC). Each cell's post-pre z is painted on its footprint. Red = excited, blue = inhibited.",
    )

    section_slide(prs, "All recordings together", "12 movies  |  t = 0 is " + cfg["t0_short"])
    for fig in figs:
        figure_slide(prs, fig["title"] + "   |   all movies", fig["what"], fig["cohort_result"], avg / fig["file"])

    figure_slide(
        prs,
        "Fig4  Spatial map by event   |   all movies",
        BY_EVENT_WHAT,
        BY_EVENT_RESULT,
        BY_EVENT / "cohort_spatial_by_event.png",
    )

    for folder, rec_title in RECS:
        d = idx.get(folder)
        for fig in figs:
            result = rec_results(d, fig["file"]) if d else fig["cohort_result"]
            figure_slide(prs, fig["title"] + "   |   " + rec_title, fig["what"], result, per / folder / fig["file"])
        figure_slide(
            prs,
            "Fig4  Spatial map by event   |   " + rec_title,
            BY_EVENT_WHAT,
            BY_EVENT_RESULT,
            BY_EVENT / "per_recording" / (folder + "_spatial_by_event.png"),
        )

    method_slide(prs, "Next", cfg["next_body"])
    written = save_pptx(prs, out)
    copy_if_possible(written, WORK / cfg["out_name"])
    return written


def build_by_event_deck():
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    slide = blank(prs)
    add_textbox(slide, 0.8, 1.8, 11.7, 1.2, "Spatial poke maps by event", 36, True, NAVY, PP_ALIGN.CENTER)
    add_textbox(
        slide, 1.2, 3.2, 10.9, 1.6,
        "2025-08-13   |   same cells, three t = 0 markers: pre / peak / post",
        18, False, GREY, PP_ALIGN.CENTER,
    )
    add_textbox(
        slide, 1.2, 4.8, 10.9, 1.2,
        "Take-home: excited and inhibited cells are mixed in every FOV; the pattern does not flip to one region when the event changes.",
        18, True, TEAL, PP_ALIGN.CENTER,
    )
    method_slide(
        prs, "How to read these maps",
        "Red = excited, blue = inhibited after the poke. Columns are pre (onset candidate), peak (brightness peak), and post (after the peak). Color scale is shared within each movie.",
    )
    figure_slide(
        prs, "All movies  |  pre vs peak vs post",
        BY_EVENT_WHAT, BY_EVENT_RESULT, BY_EVENT / "cohort_spatial_by_event.png",
    )
    for folder, rec_title in RECS:
        figure_slide(
            prs, rec_title,
            BY_EVENT_WHAT, BY_EVENT_RESULT,
            BY_EVENT / "per_recording" / (folder + "_spatial_by_event.png"),
        )
    out = BY_EVENT / "CrystalSkull_face_poke_labmeeting_spatial_by_event.pptx"
    written = save_pptx(prs, out)
    copy_if_possible(written, WORK / written.name)
    copy_if_possible(written, DECKS["peak"]["lm"] / written.name)
    return written


def replace_first_picture(slide, img_path):
    from pptx.enum.shapes import MSO_SHAPE_TYPE

    pic = None
    for shape in slide.shapes:
        if shape.shape_type == MSO_SHAPE_TYPE.PICTURE:
            pic = shape
            break
    if pic is None or img_path is None or not Path(img_path).is_file():
        return False
    left, top, width, height = pic.left, pic.top, pic.width, pic.height
    pic._element.getparent().remove(pic._element)
    slide.shapes.add_picture(str(img_path), left, top, width, height)
    return True


def replace_text(slide, pairs):
    for shape in slide.shapes:
        if not shape.has_text_frame:
            continue
        for p in shape.text_frame.paragraphs:
            full = "".join(run.text for run in p.runs)
            new = full
            for old, repl in pairs:
                new = new.replace(old, repl)
            if new != full and p.runs:
                p.runs[0].text = new
                for run in p.runs[1:]:
                    run.text = ""


def rec_folder_from_title(title):
    for folder, rec_title in RECS:
        if rec_title in title:
            return folder
    return None


def patch_existing_pptx(path, event_lm):
    path = Path(path)
    if not path.is_file():
        return None
    prs = Presentation(str(path))
    n_rep = 0
    for slide in prs.slides:
        title = ""
        for shape in slide.shapes:
            if shape.has_text_frame and shape.text_frame.text.strip().startswith("Fig4"):
                title = shape.text_frame.text.strip()
                break
        if not title:
            continue
        if "by event" in title.lower():
            folder = rec_folder_from_title(title)
            if folder is None and "all movies" in title:
                img = BY_EVENT / "cohort_spatial_by_event.png"
            elif folder is not None:
                img = BY_EVENT / "per_recording" / (folder + "_spatial_by_event.png")
            else:
                continue
        else:
            folder = rec_folder_from_title(title)
            if folder is None and "all movies" in title:
                img = event_lm / "cohort_average" / "fig4_spatial_response.png"
            elif folder is not None:
                img = event_lm / "per_recording" / folder / "fig4_spatial_response.png"
            else:
                continue
        if replace_first_picture(slide, img):
            n_rep += 1
        replace_text(
            slide,
            [
                ("Fig4  FOV anatomy  (not a response map)", "Fig4  Spatial poke-response map"),
                ("EXTRACT session-mean image: average brightness over the whole movie, not pre vs post poke.", FIG4_WHAT),
                ("Cells are dense across dorsal cortex, but this picture cannot say which regions turned on or off after poking.", FIG4_RESULT_COHORT),
                ("This is only where cells sit in this FOV — not which cells responded to the poke.", FIG4_RESULT_COHORT),
            ],
        )
    written = save_pptx(prs, path)
    print("  patched", n_rep, "Fig4 slides in", written)
    return written


def main():
    written = {}
    for key in ("peak", "pre", "post"):
        written[key] = build_deck(DECKS[key])
    written["by_event"] = build_by_event_deck()

    # Keep the user's working copies (including hansol) on the new spatial fig4.
    patch_jobs = [
        (WORK / "CrystalSkull_face_poke_labmeeting_peak_hansol.pptx", DECKS["peak"]["lm"]),
        (WORK / "CrystalSkull_face_poke_labmeeting_peak.pptx", DECKS["peak"]["lm"]),
        (WORK / "CrystalSkull_face_poke_labmeeting_pre.pptx", DECKS["pre"]["lm"]),
        (WORK / "CrystalSkull_face_poke_labmeeting_post.pptx", DECKS["post"]["lm"]),
        (DECKS["peak"]["lm"] / "CrystalSkull_face_poke_labmeeting.pptx", DECKS["peak"]["lm"]),
    ]
    for ppt, lm in patch_jobs:
        try:
            patch_existing_pptx(ppt, lm)
        except PermissionError:
            alt = ppt.with_name(ppt.stem + "_spatial.pptx")
            print("Could not patch locked file; rebuilding beside it:", alt)
        except Exception as exc:
            print("Patch skipped", ppt, exc)

    copy_if_possible(written["peak"], WORK / "CrystalSkull_face_poke_labmeeting_peak_hansol_spatial.pptx")


if __name__ == "__main__":
    main()
