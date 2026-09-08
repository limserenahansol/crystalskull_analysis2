"""Combined 20250813 vs 20250818 lab-meeting PPT (comparison first)."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

NAVY = RGBColor(0x1F, 0x3A, 0x5F)
TEAL = RGBColor(0x0D, 0x5C, 0x63)
RED = RGBColor(0x9B, 0x1B, 0x30)
GREY = RGBColor(0x55, 0x55, 0x55)

ROOT = Path(r"C:\Users\hsollim\behavior_task\crystal_skull_analysis\stim_triggered_figures")
LM13 = ROOT / "20250813_lab_meeting" / "lab_meeting"
LM18 = ROOT / "20250818_lab_meeting" / "lab_meeting"
OUT = ROOT / "combined_0813_0818"
WORK = ROOT / "20250813_lab_meeting_pre" / "lab_meeting"
OUT_NAME = "CrystalSkull_face_poke_labmeeting_0813_0818.pptx"


def load_index(lm):
    rows = []
    lines = (lm / "recording_index.csv").read_text(encoding="utf-8").strip().splitlines()
    hdr = lines[0].split(",")
    for line in lines[1:]:
        rows.append(dict(zip(hdr, line.split(","))))
    return rows


def load_cohort_csv(path):
    lines = Path(path).read_text(encoding="utf-8").strip().splitlines()
    hdr = lines[0].split(",")
    return [dict(zip(hdr, line.split(","))) for line in lines[1:]]


def stats(rows, cohort):
    p = np.array([float(r["p_pop"]) for r in rows])
    fe = np.array([float(r["frac_exc"]) for r in rows]) * 100
    fi = np.array([float(r["frac_inh"]) for r in rows]) * 100
    pk = np.array([int(r["n_pokes"]) for r in rows])
    nn = np.array([int(r["n_neurons"]) for r in rows])
    iti = np.array([float(r["iti_median"]) for r in cohort])
    amp = np.array([float(r["mean_amp"]) for r in cohort])
    return {
        "n": len(rows),
        "n_sig": int((p < 0.05).sum()),
        "fe_mean": fe.mean(),
        "fe_min": fe.min(),
        "fe_max": fe.max(),
        "fi_mean": fi.mean(),
        "fi_min": fi.min(),
        "fi_max": fi.max(),
        "pk_mean": pk.mean(),
        "pk_min": pk.min(),
        "pk_max": pk.max(),
        "nn_mean": nn.mean(),
        "iti": np.median(iti),
        "amp_mean": amp.mean(),
        "fe": fe,
        "fi": fi,
        "p": p,
    }


def write_summary_figs(s13, s18, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(11.5, 4.6), dpi=180)
    ax.axis("off")
    ax.set_title("Same analysis, two poke-timing definitions", fontsize=14, fontweight="bold", color="#1F3A5F", pad=12)
    col_labels = ["", "13 Aug  (peak_idx)", "18 Aug  (intensity DIP)"]
    cells = [
        ["t = 0", "Positive brightness peak", "Single-frame intensity dip"],
        ["Movies / mice", "12 movies, 3 mice (m1f, m2c, m3)", "6 movies, 2 mice (m1f=m200, m2c=m201)"],
        ["Population p<0.05", f"{s13['n_sig']} / {s13['n']}", f"{s18['n_sig']} / {s18['n']}"],
        ["Excited cells", f"{s13['fe_mean']:.1f}%  ({s13['fe_min']:.1f}–{s13['fe_max']:.1f})",
         f"{s18['fe_mean']:.1f}%  ({s18['fe_min']:.1f}–{s18['fe_max']:.1f})"],
        ["Inhibited cells", f"{s13['fi_mean']:.1f}%  ({s13['fi_min']:.1f}–{s13['fi_max']:.1f})",
         f"{s18['fi_mean']:.1f}%  ({s18['fi_min']:.1f}–{s18['fi_max']:.1f})"],
        ["Pokes / movie", f"{s13['pk_min']:.0f}–{s13['pk_max']:.0f}  (mean {s13['pk_mean']:.0f})",
         f"{s18['pk_min']:.0f}–{s18['pk_max']:.0f}  (mean {s18['pk_mean']:.0f})"],
        ["Median ITI", f"{s13['iti']:.2f} s", f"{s18['iti']:.2f} s"],
        ["Mean pop. amp (z)", f"{s13['amp_mean']:+.3f}", f"{s18['amp_mean']:+.3f}"],
    ]
    table = ax.table(cellText=cells, colLabels=col_labels, loc="center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.0, 1.55)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#DDDDDD")
        if r == 0:
            cell.set_facecolor("#1F3A5F")
            cell.set_text_props(color="white", weight="bold")
        elif c == 0:
            cell.set_facecolor("#EEF3F6")
            cell.set_text_props(weight="bold", color="#1F3A5F")
        else:
            cell.set_facecolor("white")
            cell.set_text_props(color="#333333")
    fig.tight_layout()
    table_path = out_dir / "comparison_table.png"
    fig.savefig(table_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    rng = np.random.default_rng(0)
    fig, ax = plt.subplots(figsize=(8.4, 4.8), dpi=180)
    x = np.array([0, 1, 3, 4], dtype=float)
    means = [s13["fe_mean"], s13["fi_mean"], s18["fe_mean"], s18["fi_mean"]]
    colors = ["#D94A3A", "#2F5FBF", "#D94A3A", "#2F5FBF"]
    ax.bar(x, means, width=0.7, color=colors, edgecolor="k", linewidth=0.6, zorder=2)
    j13 = rng.uniform(-0.12, 0.12, size=len(s13["fe"]))
    j18 = rng.uniform(-0.12, 0.12, size=len(s18["fe"]))
    ax.scatter(0 + j13, s13["fe"], c="k", s=22, zorder=3)
    ax.scatter(1 + j13, s13["fi"], c="k", s=22, zorder=3)
    ax.scatter(3 + j18, s18["fe"], c="k", s=22, zorder=3)
    ax.scatter(4 + j18, s18["fi"], c="k", s=22, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(["0813 exc", "0813 inh", "0818 exc", "0818 inh"])
    ax.set_ylabel("% of neurons")
    ax.set_ylim(0, 45)
    ax.set_title("Excited vs inhibited fractions (one dot = one movie)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.axvline(2, color="#CCCCCC", lw=1)
    fig.tight_layout()
    bars_path = out_dir / "comparison_frac_bars.png"
    fig.savefig(bars_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return table_path, bars_path


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


def method_slide(prs, title, body):
    slide = blank(prs)
    add_textbox(slide, 0.7, 2.0, 11.9, 0.7, title, 22, True, NAVY, PP_ALIGN.CENTER)
    add_textbox(slide, 1.1, 2.8, 11.1, 3.4, body, 20, False, GREY, PP_ALIGN.CENTER)


def section_slide(prs, title, subtitle):
    slide = blank(prs)
    add_textbox(slide, 0.7, 2.4, 11.9, 0.9, title, 32, True, TEAL, PP_ALIGN.CENTER)
    add_textbox(slide, 1.0, 3.5, 11.3, 1.4, subtitle, 18, False, GREY, PP_ALIGN.CENTER)


def picture_fit(slide, img_path, left, top, max_w, max_h):
    img_path = Path(img_path)
    if not img_path.is_file():
        add_textbox(slide, left, top + max_h / 2, max_w, 0.5, "Figure missing", 16, True, RED, PP_ALIGN.CENTER)
        return
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


def two_fig_slide(prs, title, what, result, img_l, cap_l, img_r, cap_r):
    slide = blank(prs)
    add_textbox(slide, 0.35, 0.12, 12.6, 0.42, title, 22, True, NAVY)
    tf = add_textbox(slide, 0.35, 0.52, 12.6, 0.38, "What:  " + what, 13, False, GREY)
    p = tf.add_paragraph()
    set_run(p, "Result:  " + result, 14, True, TEAL)
    add_textbox(slide, 0.35, 0.98, 6.2, 0.32, cap_l, 14, True, NAVY, PP_ALIGN.CENTER)
    add_textbox(slide, 6.8, 0.98, 6.2, 0.32, cap_r, 14, True, NAVY, PP_ALIGN.CENTER)
    picture_fit(slide, img_l, 0.25, 1.32, 6.4, 5.9)
    picture_fit(slide, img_r, 6.7, 1.32, 6.4, 5.9)


def figure_slide(prs, title, what, result, img_path):
    slide = blank(prs)
    add_textbox(slide, 0.4, 0.14, 12.5, 0.42, title, 22, True, NAVY)
    tf = add_textbox(slide, 0.4, 0.56, 12.5, 0.4, "What:  " + what, 13, False, GREY)
    p = tf.add_paragraph()
    set_run(p, "Result:  " + result, 14, True, TEAL)
    picture_fit(slide, img_path, 0.35, 1.35, 12.6, 5.9)


def bullets(prs, title, items):
    slide = blank(prs)
    add_textbox(slide, 0.6, 0.35, 12.1, 0.55, title, 24, True, NAVY)
    box = slide.shapes.add_textbox(Inches(0.9), Inches(1.15), Inches(11.5), Inches(5.8))
    tf = box.text_frame
    tf.word_wrap = True
    for i, item in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        set_run(p, "•  " + item, 18, False, GREY)
        p.space_after = Pt(14)
    return slide


def write_n18_two_slide_deck():
    n18 = OUT / "n18"
    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)
    figure_slide(
        prs,
        "Population peri-poke response   (n = 18 recordings)",
        "Each gray line is one recording. Red / black = mean ± SEM across recordings. Peri-window -2 to +3 s around the poke.",
        "The mean trace is poke-locked: a transient at t = 0, a brief dip, then a rebound. Recordings share the same time course.",
        n18 / "fig1_n18_psth.png",
    )
    figure_slide(
        prs,
        "Poke-responsive cells   (n = 18 recordings)",
        "Left: mean STA of excited, inhibited, and all cells. Right: % of neurons in each recording (one dot = one recording). A cell is excited or inhibited if 0.4–1.2 s differs from -2–0 s (signed-rank, p<0.05).",
        "Mean 8.8% excited and 28.2% inhibited. The all-cell STA stays near zero because the two groups cancel. 14 of 18 recordings p<0.05.",
        n18 / "fig2_n18_cells.png",
    )
    paths = [
        OUT / "CrystalSkull_face_poke_n18.pptx",
        WORK / "CrystalSkull_face_poke_n18.pptx",
    ]
    blob = None
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        try:
            if blob is None:
                prs.save(str(path))
                blob = path.read_bytes()
            else:
                path.write_bytes(blob)
            print("Wrote", path, "slides", len(prs.slides))
        except PermissionError:
            alt = path.with_name(path.stem + "_copy.pptx")
            if blob is None:
                prs.save(str(alt))
                blob = alt.read_bytes()
            else:
                alt.write_bytes(blob)
            print("Locked, wrote", alt)


def append_n18_to_existing():
    n18 = OUT / "n18"
    candidates = [
        OUT / "CrystalSkull_face_poke_labmeeting_0813_0818_pooled.pptx",
        OUT / "CrystalSkull_face_poke_labmeeting_0813_0818.pptx",
        WORK / "CrystalSkull_face_poke_labmeeting_0813_0818_pooled.pptx",
        WORK / "CrystalSkull_face_poke_labmeeting_0813_0818.pptx",
    ]
    for src in candidates:
        if not src.is_file():
            continue
        try:
            prs = Presentation(str(src))
        except Exception:
            continue
        figure_slide(
            prs,
            "Population peri-poke response   (n = 18 recordings)",
            "Each gray line is one recording. Red / black = mean ± SEM across recordings. Peri-window -2 to +3 s around the poke.",
            "The mean trace is poke-locked: a transient at t = 0, a brief dip, then a rebound. Recordings share the same time course.",
            n18 / "fig1_n18_psth.png",
        )
        figure_slide(
            prs,
            "Poke-responsive cells   (n = 18 recordings)",
            "Left: mean STA of excited, inhibited, and all cells. Right: % of neurons in each recording (one dot = one recording). A cell is excited or inhibited if 0.4–1.2 s differs from -2–0 s (signed-rank, p<0.05).",
            "Mean 8.8% excited and 28.2% inhibited. The all-cell STA stays near zero because the two groups cancel. 14 of 18 recordings p<0.05.",
            n18 / "fig2_n18_cells.png",
        )
        try:
            prs.save(str(src))
            print("Appended n=18 slides to", src, "now", len(prs.slides), "slides")
        except PermissionError:
            alt = src.with_name(src.stem + "_n18slides.pptx")
            prs.save(str(alt))
            print("Locked, wrote", alt, "slides", len(prs.slides))


def main():
    write_n18_two_slide_deck()
    OUT.mkdir(parents=True, exist_ok=True)
    idx13 = load_index(LM13)
    idx18 = load_index(LM18)
    coh13 = load_cohort_csv(ROOT / "20250813_lab_meeting" / "cohort_stim_summary.csv")
    coh18 = load_cohort_csv(ROOT / "20250818_lab_meeting" / "cohort_stim_summary.csv")
    s13 = stats(idx13, coh13)
    s18 = stats(idx18, coh18)
    table_path, bars_path = write_summary_figs(s13, s18, OUT)
    csv_path = OUT / "pooled_stim_summary.csv"
    with csv_path.open("w", encoding="utf-8") as f:
        f.write("day,label,movie,n_neurons,n_pokes,frac_exc,frac_inh,p_pop\n")
        for day, rows in (("0813", idx13), ("0818", idx18)):
            for r in rows:
                f.write(
                    f"{day},{r['label']},{r['movie']},{r['n_neurons']},{r['n_pokes']},"
                    f"{r['frac_exc']},{r['frac_inh']},{r['p_pop']}\n"
                )

    avg13 = LM13 / "cohort_average"
    avg18 = LM18 / "cohort_average"
    map13 = LM13 / "aligned_per_recording" / "cohort_aligned_sig.png"
    map18 = LM18 / "aligned_per_recording" / "cohort_aligned_sig.png"
    ap13 = avg13 / "fig4b_anterior_posterior_bias.png"
    ap_both = OUT / "ap_bias_both_days.png"

    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)

    slide = blank(prs)
    add_textbox(slide, 0.7, 1.55, 11.9, 1.1, "Face poking on crystal skull", 36, True, NAVY, PP_ALIGN.CENTER)
    add_textbox(
        slide, 1.1, 2.85, 11.1, 1.0,
        "Pooled analysis  |  18 movies, 43,156 neurons, 1,098 pokes  |  13 Aug peak + 18 Aug dip",
        18, False, GREY, PP_ALIGN.CENTER,
    )
    add_textbox(
        slide, 1.2, 4.05, 10.9, 2.0,
        "Take-home: pooling all movies, 14/18 are significant; 8.6% of cells are excited and 27.6% inhibited. The all-cell mean stays near zero because the two groups cancel. The pooled PSTH is a mix of two t=0 definitions, so day means are shown with the grand mean.",
        18, True, TEAL, PP_ALIGN.CENTER,
    )

    method_slide(
        prs, "Experiment",
        "Whole dorsal cortex Ca2+ imaging while a motor arm pokes the face. Both days are poking sessions (no no-poke baseline). m1f and m2c were imaged on both days; m3 only on 13 Aug. mirror1 / mirror2 are two optical paths.",
    )

    n18 = OUT / "n18"
    figure_slide(
        prs,
        "Population peri-poke response   (n = 18 recordings)",
        "Each gray line is one recording. Red / black = mean ± SEM across recordings. Peri-window -2 to +3 s around the poke.",
        "The mean trace is poke-locked: a transient at t = 0, a brief dip, then a rebound. Recordings share the same time course.",
        n18 / "fig1_n18_psth.png",
    )
    figure_slide(
        prs,
        "Poke-responsive cells   (n = 18 recordings)",
        "Left: mean STA of excited, inhibited, and all cells. Right: % of neurons in each recording (one dot = one recording). A cell is excited or inhibited if 0.4–1.2 s differs from -2–0 s (signed-rank, p<0.05).",
        "Mean 8.8% excited and 28.2% inhibited. The all-cell STA stays near zero because the two groups cancel. 14 of 18 recordings p<0.05.",
        n18 / "fig2_n18_cells.png",
    )

    method_slide(
        prs, "The only intended difference is t = 0",
        "No behavior camera. A blocker on the arm changes mean-frame intensity. 13 Aug peak_idx = positive brightness peak. 18 Aug peak_idx = single-frame intensity DIP (Omer). Analysis windows are identical: -2 to +3 s peri, 0.4–1.2 s vs -2–0 s, Wilcoxon signed-rank.",
    )
    method_slide(
        prs, "How traces and maps are made",
        "EXTRACT temporal_weights at 33 Hz. Spatial maps use Omer's EXTRACT spatial_weights (footprints reconstructed from the ndSparse CSC in each *_ext.mat). Red = excited, blue = inhibited. mirror2 is flipped so anterior is on the RIGHT (P left, A right).",
    )

    pooled = OUT / "pooled_average"
    section_slide(prs, "All 18 movies together", "Equal-movie grand mean  |  12 from 13 Aug + 6 from 18 Aug")
    method_slide(
        prs, "How pooling is done",
        "Each movie contributes one PSTH / STA. The black line is the mean of those 18 traces (equal movie weight). Neuron-weighted and trial-weighted means overlay it. Cells are also pooled: 8.6% excited and 27.6% inhibited of 43,156 neurons. FOVs cannot be averaged, so Fig4 is the full set of aligned maps.",
    )
    figure_slide(
        prs,
        "Pooled Fig1  |  Average trace of all 18 movies",
        "Gray/tint = each movie. Orange = 13 Aug mean. Teal = 18 Aug mean. Black = mean of all 18.",
        "The pooled black line has a small peak near t=0, a dip, and a rebound. That shape is a blend: 13 Aug locks to the peak, 18 Aug rises before the dip. It is not a third biological condition.",
        pooled / "fig1_population_psth.png",
    )
    figure_slide(
        prs,
        "Pooled Fig1b  |  Weighting check",
        "Same 18 movies averaged with equal movie weight, neuron weight, or trial weight.",
        "The three pooled curves sit on top of each other. The grand mean is not an artifact of one large movie.",
        pooled / "fig1b_weighting.png",
    )
    figure_slide(
        prs,
        "Pooled Fig2  |  Average STA of all 18 movies",
        "Mean of per-movie excited / inhibited / all-cell STAs. Dotted = 13 Aug, dashed = 18 Aug.",
        "Pooled excited cells rise after t=0; pooled inhibited cells fall. All-cells (black) stay near zero because the groups cancel.",
        pooled / "fig2_neuron_sta.png",
    )
    figure_slide(
        prs,
        "Pooled Fig3  |  All cells from all movies",
        "Bars = mean ± SD across 18 movies. Yellow diamond = all 43,156 neurons pooled. Orange circles = 13 Aug, teal squares = 18 Aug.",
        "Movie-mean 8.8% exc / 28.2% inh; neuron-pooled 8.6% / 27.6%. 14/18 movies p<0.05.",
        pooled / "fig3_responsive_cells.png",
    )
    figure_slide(
        prs,
        "Pooled Fig4  |  All FOVs (cannot pixel-average)",
        "Aligned EXTRACT footprint maps for every recording in the pool. P left, A right.",
        "The repeating pattern is mixed red/blue with anterior more excited. There is no single average cortex image because FOVs are not registered.",
        pooled / "fig4_spatial_aligned.png",
    )
    figure_slide(
        prs,
        "Pooled Fig5  |  All pokes together",
        "Orange bars = 13 Aug movies, teal = 18 Aug. Histogram = every inter-poke interval from all 18 movies.",
        "1,098 detected pokes. Pooled ITI median 5.24 s (13 Aug 5.30 s dominates the median; 18 Aug is 3.36 s).",
        pooled / "fig5_trial_timing.png",
    )

    section_slide(prs, "The same numbers, split by day", "Because t = 0 is not the same physical marker")

    figure_slide(
        prs,
        "Numbers at a glance",
        "Same cell-classification windows on both days. 18 Aug has fewer movies, a shorter ITI, and a dip-based t = 0.",
        f"Population lock is tighter on 13 Aug ({s13['n_sig']}/{s13['n']} movies) than 18 Aug ({s18['n_sig']}/{s18['n']}), but the exc < inh split is the same direction both days.",
        table_path,
    )
    figure_slide(
        prs,
        "Excited vs inhibited fractions",
        "Each black dot is one movie. Red bars = mean % excited; blue bars = mean % inhibited.",
        "Inhibition dominates on both days (~6–14% excited vs ~19–37% inhibited). 18 Aug m201 is especially inhibition-heavy.",
        bars_path,
    )

    two_fig_slide(
        prs,
        "Fig1  Population PSTH",
        "Mean population z and fraction active, aligned to each poke.",
        "13 Aug: sharp transient at t=0, dip, rebound ~0.8 s. 18 Aug: activity already rising before t=0, then a crash at the dip.",
        avg13 / "fig1_population_psth.png", "13 Aug  |  t = 0 is the brightness peak",
        avg18 / "fig1_population_psth.png", "18 Aug  |  t = 0 is the intensity dip",
    )
    two_fig_slide(
        prs,
        "Fig2  Neuron STA",
        "Cohort-mean stimulus-triggered averages of excited, inhibited, and all cells.",
        "Both days split into an up group and a down group. The all-cell mean is modest because they cancel — more completely on 18 Aug.",
        avg13 / "fig2_neuron_sta.png", "13 Aug",
        avg18 / "fig2_neuron_sta.png", "18 Aug",
    )
    two_fig_slide(
        prs,
        "Fig3  Excited vs inhibited cells",
        "Per neuron: 0.4–1.2 s after t=0 vs -2–0 s before; signed-rank p<0.05.",
        f"13 Aug: {s13['n_sig']}/{s13['n']} movies significant. 18 Aug: {s18['n_sig']}/{s18['n']} (both m201 s00, plus m201 s01 mirror1). m200 on 18 Aug is not significant (missed dips).",
        avg13 / "fig3_responsive_cells.png", "13 Aug",
        avg18 / "fig3_responsive_cells.png", "18 Aug",
    )
    two_fig_slide(
        prs,
        "Fig4  Spatial maps  (aligned, significant cells)",
        "EXTRACT footprints x poke amplitude. mirror2 flipped. P = posterior (left), A = anterior (right).",
        "Cells are locally intermingled on both days. After the same AP alignment, anterior tends to be more excited and posterior more inhibited.",
        map13, "13 Aug  |  12 FOVs",
        map18, "18 Aug  |  6 FOVs",
    )
    two_fig_slide(
        prs,
        "Fig5  Poke counts and interval",
        "Detected events per movie and the histogram of time between events.",
        "13 Aug median ITI 5.30 s (~48–57 pokes). 18 Aug median 3.36 s; m201 is dense (~83–91), m200 has ~20 s gaps from missed dips.",
        avg13 / "fig5_trial_timing.png", "13 Aug",
        avg18 / "fig5_trial_timing.png", "18 Aug",
    )

    if ap13.is_file():
        figure_slide(
            prs,
            "Anterior–posterior bias  |  13 Aug",
            "Spearman ρ between each cell's LR position and poke amplitude. After flipping mirror2, anterior is to the right in every movie.",
            "12/12 movies ρ > 0 (median 0.13, signrank p = 4.9e-4): anterior more excited, posterior more inhibited. Not a CCF registration — AP is the LR axis along the hemispheres.",
            ap13,
        )
    if ap_both.is_file():
        figure_slide(
            prs,
            "Anterior–posterior bias  |  both days",
            "Same centroid vs amplitude Spearman ρ after the mirror2 flip.",
            "If 18 Aug points also sit above zero, the AP bias is not special to the brightness-peak alignment.",
            ap_both,
        )

    bullets(
        prs,
        "Interpretation",
        [
            "Pooled across 18 movies / 43,156 cells: 8.6% excited, 27.6% inhibited, 14/18 movies significant. This is the combined biology.",
            "The pooled black PSTH is the average of two clocks (peak vs dip). Read orange vs teal to see timing; read the yellow-diamond cell counts to see the shared exc/inh split.",
            "All-cell STA near zero is cancellation, not absence of a poke response. Excited and inhibited groups are both locked; they have opposite sign.",
            "FOVs cannot be averaged. The pooled spatial result is that the same mixed AP pattern appears in the full set of 18 maps.",
            "Skin-contact time is still not proven. Do not treat the pooled t=0 as a single physical onset.",
        ],
    )
    bullets(
        prs,
        "What we should not over-claim",
        [
            "FOV is not CCF/bregma registered. 'Anterior' means the right side after aligning mirror2 to mirror1, along the hemispheres.",
            "Power is labeled only on 13 Aug (s00 = 8.38, s01 = 20 mW/cm2). Do not copy those values onto 18 Aug.",
            "m200 on 18 Aug failing the population test is consistent with missed dips, not evidence that this mouse has no poke map (its spatial map still shows mixed red/blue).",
            "We reconstructed EXTRACT spatial_weights from the file (CSC). That is Omer's footprint matrix, not a new segmentation.",
        ],
    )

    section_slide(prs, "13 Aug  in one place", "t = 0 is peak_idx (positive brightness peak)  |  12 movies")
    figure_slide(prs, "13 Aug  |  Fig1  Population PSTH", "Aligned to the brightness peak.", "Transient at t=0, brief dip, rebound ~0.8 s in all 12 movies.", avg13 / "fig1_population_psth.png")
    figure_slide(prs, "13 Aug  |  Fig2  Neuron STA", "Excited vs inhibited vs all cells.", "The two groups partly cancel in the population mean.", avg13 / "fig2_neuron_sta.png")
    figure_slide(prs, "13 Aug  |  Fig3  Responsive cells", "0.4–1.2 s vs -2–0 s, signed-rank.", f"{s13['n_sig']}/{s13['n']} movies p<0.05; ~{s13['fe_mean']:.0f}% exc, ~{s13['fi_mean']:.0f}% inh.", avg13 / "fig3_responsive_cells.png")
    figure_slide(prs, "13 Aug  |  Fig4  Aligned spatial maps", "P left, A right. Red = excited, blue = inhibited.", "Intermingled cells; anterior-excited / posterior-inhibited after the mirror2 flip.", map13)
    figure_slide(prs, "13 Aug  |  Fig5  Trial timing", "Pokes per movie and ITI.", "Median interval 5.30 s; about 48–57 pokes per movie.", avg13 / "fig5_trial_timing.png")

    section_slide(prs, "18 Aug  in one place", "t = 0 is peak_idx (intensity DIP)  |  6 movies")
    figure_slide(prs, "18 Aug  |  Fig1  Population PSTH", "Aligned to the intensity dip.", "Rise into t=0, sharp drop at the dip, slow recovery — same shape in all 6 movies.", avg18 / "fig1_population_psth.png")
    figure_slide(prs, "18 Aug  |  Fig2  Neuron STA", "Excited vs inhibited vs all cells.", "Same up/down split; all-cell mean is close to zero because inhibition wins on m201.", avg18 / "fig2_neuron_sta.png")
    figure_slide(prs, "18 Aug  |  Fig3  Responsive cells", "Same 0.4–1.2 s window as 13 Aug.", f"{s18['n_sig']}/{s18['n']} movies p<0.05; m201 carries the effect, m200 does not.", avg18 / "fig3_responsive_cells.png")
    figure_slide(prs, "18 Aug  |  Fig4  Aligned spatial maps", "P left, A right. Same EXTRACT footprint method as 13 Aug.", "Still mixed red/blue, with the same AP direction after flipping mirror2.", map18)
    figure_slide(prs, "18 Aug  |  Fig5  Trial timing", "Detected dips and ITI.", "Median 3.36 s. m200 ~43–45 dips with ~20 s gaps; m201 ~83–91.", avg18 / "fig5_trial_timing.png")

    method_slide(
        prs, "Next",
        "The pooled black line is the biology that repeats (exc/inh split, mixed maps). Keep day means when discussing timing (peak vs dip). A shared physical marker would be needed before treating all 18 t=0 as one onset.",
    )

    out = OUT / OUT_NAME
    try:
        prs.save(str(out))
        print("Wrote", out, "slides", len(prs.slides))
        written = out
    except PermissionError:
        written = OUT / "CrystalSkull_face_poke_labmeeting_0813_0818_with_n18.pptx"
        prs.save(str(written))
        print("Locked original, wrote", written, "slides", len(prs.slides))
    WORK.mkdir(parents=True, exist_ok=True)
    dest = WORK / written.name
    try:
        dest.write_bytes(written.read_bytes())
        print("Copied", dest)
    except PermissionError:
        alt = WORK / (written.stem + "_copy.pptx")
        alt.write_bytes(written.read_bytes())
        print("Locked, copied", alt)


if __name__ == "__main__":
    main()
