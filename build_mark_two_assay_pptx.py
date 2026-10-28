"""Short Mark deck: face poke (per-mouse stats) + melittin / morphine (same visual language)."""

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
GREY = RGBColor(0x55, 0x55, 0x55)
RED = RGBColor(0x9B, 0x1B, 0x30)

ROOT = Path(r"C:\Users\hsollim\behavior_task\crystal_skull_analysis")
POKE = ROOT / "stim_triggered_figures" / "combined_0813_0818"
MEL = ROOT / "melittin_lab_meeting"
OUT_DIR = ROOT / "mark_brief"
OUT = OUT_DIR / "CrystalSkull_facepoke_vs_melittin_Mark.pptx"

IDX13 = ROOT / "stim_triggered_figures" / "20250813_lab_meeting" / "lab_meeting" / "recording_index.csv"
IDX18 = ROOT / "stim_triggered_figures" / "20250818_lab_meeting" / "lab_meeting" / "recording_index.csv"
MEL_IDX = MEL / "recording_index.csv"


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


def method_slide(prs, title, body):
    slide = blank(prs)
    add_textbox(slide, 0.7, 1.9, 11.9, 0.7, title, 22, True, NAVY, PP_ALIGN.CENTER)
    add_textbox(slide, 1.1, 2.7, 11.1, 3.6, body, 18, False, GREY, PP_ALIGN.CENTER)


def section_slide(prs, title, subtitle):
    slide = blank(prs)
    add_textbox(slide, 0.7, 2.5, 11.9, 0.9, title, 32, True, TEAL, PP_ALIGN.CENTER)
    add_textbox(slide, 1.0, 3.5, 11.3, 1.5, subtitle, 18, False, GREY, PP_ALIGN.CENTER)


def figure_slide(prs, title, what, result, img_path):
    slide = blank(prs)
    add_textbox(slide, 0.4, 0.12, 12.5, 0.38, title, 20, True, NAVY)
    tf = add_textbox(slide, 0.4, 0.48, 12.5, 0.34, "What:  " + what, 12, False, GREY)
    p = tf.add_paragraph()
    set_run(p, "Result:  " + result, 13, True, TEAL)
    picture_fit(slide, img_path, 0.30, 0.95, 12.7, 6.3)


def two_fig_slide(prs, title, what, result, img_l, img_r):
    slide = blank(prs)
    add_textbox(slide, 0.30, 0.10, 12.7, 0.36, title, 18, True, NAVY)
    tf = add_textbox(slide, 0.30, 0.44, 12.7, 0.32, "What:  " + what, 12, False, GREY)
    p = tf.add_paragraph()
    set_run(p, "Result:  " + result, 13, True, TEAL)
    picture_fit(slide, img_l, 0.20, 0.90, 6.5, 6.35)
    picture_fit(slide, img_r, 6.65, 0.90, 6.5, 6.35)


def bullets(prs, title, items):
    slide = blank(prs)
    add_textbox(slide, 0.6, 0.35, 12.1, 0.55, title, 24, True, NAVY)
    box = slide.shapes.add_textbox(Inches(0.9), Inches(1.15), Inches(11.5), Inches(5.8))
    tf = box.text_frame
    tf.word_wrap = True
    for i, item in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        set_run(p, "•  " + item, 18, False, GREY)
        p.space_after = Pt(12)


def load_csv(path):
    if not Path(path).is_file():
        return []
    lines = Path(path).read_text(encoding="utf-8").strip().splitlines()
    hdr = lines[0].split(",")
    return [dict(zip(hdr, line.split(","))) for line in lines[1:]]


def mouse_name(raw):
    if raw in ("m3_m1f", "m3"):
        return "m3"
    return raw


def poke_per_mouse(rows):
    mice = {}
    for r in rows:
        m = mouse_name(r["mouse"])
        n = int(r["n_neurons"])
        n_exc = int(round(float(r["frac_exc"]) * n))
        n_inh = int(round(float(r["frac_inh"]) * n))
        d = mice.setdefault(m, {"n": 0, "exc": 0, "inh": 0, "movies": 0, "sig": 0})
        d["n"] += n
        d["exc"] += n_exc
        d["inh"] += n_inh
        d["movies"] += 1
        if float(r["p_pop"]) < 0.05:
            d["sig"] += 1
    order = [m for m in ("m1f", "m2c", "m3") if m in mice]
    return [(m, mice[m]) for m in order]


def write_poke_mouse_fig(pairs, out_path):
    names = [p[0] for p in pairs]
    fe = [100.0 * p[1]["exc"] / p[1]["n"] for p in pairs]
    fi = [100.0 * p[1]["inh"] / p[1]["n"] for p in pairs]
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.2), dpi=180, gridspec_kw={"width_ratios": [1.05, 1.25]})
    ax = axes[0]
    x = np.arange(len(names))
    ax.bar(x - 0.18, fe, 0.36, color="#D94A3A", edgecolor="k", linewidth=0.6, label="Excited")
    ax.bar(x + 0.18, fi, 0.36, color="#2F5FBF", edgecolor="k", linewidth=0.6, label="Inhibited")
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.set_ylabel("% of neurons (pooled within mouse)")
    ax.set_ylim(0, 45)
    ax.set_title("Excited vs inhibited by mouse")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax2 = axes[1]
    ax2.axis("off")
    cells = [["Mouse", "Movies", "Neurons", "Excited", "Inhibited"]]
    for m, d in pairs:
        cells.append(
            [
                m,
                f"{d['movies']}  ({d['sig']} p<0.05)",
                f"{d['n']:,}",
                f"{d['exc']:,}   ({100 * d['exc'] / d['n']:.1f}%)",
                f"{d['inh']:,}   ({100 * d['inh'] / d['n']:.1f}%)",
            ]
        )
    tot_n = sum(d["n"] for _, d in pairs)
    tot_e = sum(d["exc"] for _, d in pairs)
    tot_i = sum(d["inh"] for _, d in pairs)
    tot_m = sum(d["movies"] for _, d in pairs)
    tot_s = sum(d["sig"] for _, d in pairs)
    cells.append(
        [
            "All",
            f"{tot_m}  ({tot_s} p<0.05)",
            f"{tot_n:,}",
            f"{tot_e:,}   ({100 * tot_e / tot_n:.1f}%)",
            f"{tot_i:,}   ({100 * tot_i / tot_n:.1f}%)",
        ]
    )
    table = ax2.table(cellText=cells, loc="center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.05, 1.85)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#DDDDDD")
        if r == 0:
            cell.set_facecolor("#1F3A5F")
            cell.set_text_props(color="white", weight="bold")
        elif r == len(cells) - 1:
            cell.set_facecolor("#EEF3F6")
            cell.set_text_props(weight="bold", color="#1F3A5F")
        else:
            cell.set_facecolor("white")
    ax2.set_title("Pooled cell counts  |  0.4–1.2 s vs −2–0 s, p<0.05", fontsize=12, color="#1F3A5F", pad=8)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return tot_n, tot_e, tot_i, tot_m, tot_s


def write_melittin_change_fig(rows, out_path):
    only = [r for r in rows if r.get("cohort") == "melittin_only"]
    morph = [r for r in rows if r.get("cohort") == "melittin_morphine"]
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.4), dpi=180, sharey=True)
    specs = [
        (axes[0], only, "Melittin only", "#0D8A8A", "8 recordings  |  CS1014 + CS0204"),
        (axes[1], morph, "Melittin + morphine", "#E07A2F", "4 recordings  |  paired CS1014 mice"),
    ]
    for ax, recs, title, color, sub in specs:
        yb, yi = [], []
        for i, r in enumerate(recs):
            b = float(r["mean_dff_bl"])
            inj = float(r["mean_dff_inj"])
            yb.append(b)
            yi.append(inj)
            ax.plot([0, 1], [b, inj], "-", color="#BBBBBB", lw=1.2, zorder=2)
            ax.scatter([0, 1], [b, inj], c=[color, color], s=36, zorder=3, edgecolors="k", linewidths=0.4)
        if yb:
            ax.plot([0, 1], [np.mean(yb), np.mean(yi)], "-", color=color, lw=3.0, zorder=4)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Baseline", "Injection"])
        ax.set_title(title + "\n" + sub, fontsize=12)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        n_neg = sum(float(r["d_pop_dff"]) < 0 for r in recs)
        ax.text(0.5, 0.03, f"{n_neg}/{len(recs)} recordings  Δ < 0", transform=ax.transAxes, ha="center", fontsize=10, color="#555555")
    axes[0].set_ylabel("Population mean EXTRACT dF/F-like")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_melittin_mouse_fig(rows, out_path):
    mice = {}
    for r in rows:
        m = r["mouse"]
        d = mice.setdefault(
            m,
            {
                "cohort": r["cohort"],
                "n": 0,
                "up": 0,
                "down": 0,
                "d": [],
                "movies": 0,
                "sig": 0,
            },
        )
        n = int(r["n_neurons_inj"])
        d["n"] += n
        d["up"] += int(r["n_up"])
        d["down"] += int(r["n_down"])
        d["d"].append(float(r["d_pop_dff"]))
        d["movies"] += 1
        if float(r["p_pop"]) < 0.05:
            d["sig"] += 1
    order = ["CS1014-1b", "CS1014-2w", "CS0204-1b", "CS0204-2w", "1014m01b", "1014m02w"]
    order = [m for m in order if m in mice]
    fig, ax = plt.subplots(figsize=(12.4, 5.2), dpi=180)
    ax.axis("off")
    cells = [["Mouse", "Assay", "Planes", "INJ neurons", "Up vs BL mean", "Down vs BL mean", "Mean Δ dF/F"]]
    for m in order:
        d = mice[m]
        assay = "melittin only" if d["cohort"] == "melittin_only" else "melittin + morphine"
        cells.append(
            [
                m,
                assay,
                str(d["movies"]),
                f"{d['n']:,}",
                f"{d['up']:,}  ({100 * d['up'] / d['n']:.1f}%)",
                f"{d['down']:,}  ({100 * d['down'] / d['n']:.1f}%)",
                f"{np.mean(d['d']):+.2f}",
            ]
        )
    table = ax.table(cellText=cells, loc="center", cellLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.05, 1.7)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#DDDDDD")
        if r == 0:
            cell.set_facecolor("#1F3A5F")
            cell.set_text_props(color="white", weight="bold")
        elif r >= 1 and order[r - 1] in ("1014m01b", "1014m02w"):
            cell.set_facecolor("#FFF4E8")
        else:
            cell.set_facecolor("white")
    ax.set_title(
        "Melittin  |  injection cells vs baseline session mean (not the same cell IDs).  Δ = INJ − BL population mean.",
        fontsize=11,
        color="#1F3A5F",
        pad=10,
    )
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    poke_rows = load_csv(IDX13) + load_csv(IDX18)
    for r in poke_rows:
        r["mouse"] = mouse_name(r["mouse"])
    pairs = poke_per_mouse(poke_rows)
    poke_fig = OUT_DIR / "facepoke_per_mouse.png"
    tot_n, tot_e, tot_i, tot_m, tot_s = write_poke_mouse_fig(pairs, poke_fig)

    mel_rows = load_csv(MEL_IDX)
    change_fig = OUT_DIR / "melittin_bl_vs_inj.png"
    mouse_fig = OUT_DIR / "melittin_per_mouse.png"
    if mel_rows:
        write_melittin_change_fig(mel_rows, change_fig)
        write_melittin_mouse_fig(mel_rows, mouse_fig)

    only = [r for r in mel_rows if r.get("cohort") == "melittin_only"]
    morph = [r for r in mel_rows if r.get("cohort") == "melittin_morphine"]
    n_only_neg = sum(float(r["d_pop_dff"]) < 0 for r in only)
    n_morph_neg = sum(float(r["d_pop_dff"]) < 0 for r in morph)

    n18 = POKE / "n18"
    pooled = POKE / "pooled_average"
    avg = MEL / "cohort_average"
    fig1_only = avg / "fig1_population_psth_melittin_only.png"
    fig1_morph = avg / "fig1_population_psth_morphine.png"
    fig3_only = avg / "fig3_responsive_cells_melittin_only.png"
    fig3_morph = avg / "fig3_responsive_cells_morphine.png"
    if not fig1_only.is_file():
        fig1_only = avg / "fig1_population_psth.png"
    if not fig1_morph.is_file():
        fig1_morph = avg / "fig1_population_psth.png"
    if not fig3_only.is_file():
        fig3_only = avg / "fig3_responsive_cells.png"
    if not fig3_morph.is_file():
        fig3_morph = avg / "fig3_responsive_cells.png"

    prs = Presentation()
    prs.slide_width = Inches(13.333)
    prs.slide_height = Inches(7.5)

    slide = blank(prs)
    add_textbox(slide, 0.7, 1.5, 11.9, 1.0, "Two crystal-skull assays", 34, True, NAVY, PP_ALIGN.CENTER)
    add_textbox(
        slide, 1.1, 2.7, 11.1, 0.9,
        "Face poking  vs  melittin  (± morphine)",
        20, False, GREY, PP_ALIGN.CENTER,
    )
    add_textbox(
        slide, 1.2, 3.9, 10.9, 2.2,
        "Same visual language: population traces, up vs down cells, FOV maps. Face poke is peri-event. Melittin is baseline vs injection (not time-locked, cells not matched). For SBI: which paradigm?",
        18, True, TEAL, PP_ALIGN.CENTER,
    )

    method_slide(
        prs,
        "How to read the two assays",
        "Face poke (13+18 Aug): 18 movies, 3 mice. t=0 is the poke. Excited/inhibited = 0.4–1.2 s vs −2–0 s, signed-rank p<0.05. No drug. Melittin: 12 BL/INJ recordings. Both arms cropped to the shorter duration. Up/down = injection cell vs baseline session mean (independent EXTRACT IDs). Melittin+morphine uses the same CS1014 mice on a different day.",
    )

    section_slide(prs, "1.  Face poking", "No drug  |  poke-locked  |  18 movies")

    figure_slide(
        prs,
        "Face poke  |  population peri-poke (n = 18)",
        "Gray = one movie. Red/black = mean ± SEM. Peri-window −2 to +3 s.",
        "A brief transient, a dip, then a rebound. 14/18 movies p<0.05. The all-cell mean is modest because excited and inhibited cells cancel.",
        n18 / "fig1_n18_psth.png",
    )
    figure_slide(
        prs,
        "Face poke  |  excited vs inhibited cells",
        "Left: mean STA. Right: % of neurons per movie. Same 0.4–1.2 s vs −2–0 s rule.",
        f"Movie-mean ~8.8% excited / 28.2% inhibited. Pooled: {tot_e:,} / {tot_n:,} excited ({100 * tot_e / tot_n:.1f}%) and {tot_i:,} inhibited ({100 * tot_i / tot_n:.1f}%).",
        n18 / "fig2_n18_cells.png",
    )
    figure_slide(
        prs,
        "Face poke  |  how many cells per mouse",
        "Pooled within mouse across movies (both days). Excited/inhibited as above.",
        f"{tot_s}/{tot_m} movies p<0.05. Every mouse has more inhibited than excited cells.",
        poke_fig,
    )
    two_fig_slide(
        prs,
        "Face poke  |  where they sit",
        "Left: all FOVs (P left, A right). Right: AP bias of poke amplitude. Not CCF/Allen registered.",
        "Excited and inhibited cells are intermingled. 18/18 recordings: more excitation anterior, more inhibition posterior.",
        pooled / "fig4_spatial_aligned.png",
        n18 / "takehome_anterior_posterior.png",
    )

    section_slide(prs, "2.  Melittin  ± morphine", "Not poke-locked  |  BL vs injection  |  matched to the shorter duration")

    two_fig_slide(
        prs,
        "Melittin only  |  baseline vs injection  (n = 8)",
        "Same layout as face-poke Fig1/Fig3. Time is the recording, not a poke. Gray-tint = each recording.",
        f"Population mean activity falls after injection in {n_only_neg}/8 recordings. Down cells outnumber up cells, as in face poke.",
        fig1_only,
        fig3_only,
    )
    two_fig_slide(
        prs,
        "Melittin + morphine  |  baseline vs injection  (n = 4)",
        "Same CS1014 mice as the melittin-only pair, different day. Same crop-to-shorter-duration rule.",
        f"Direction splits: {n_morph_neg}/4 recordings Δ < 0. This is exploratory (small n).",
        fig1_morph,
        fig3_morph,
    )
    figure_slide(
        prs,
        "Three cases  |  baseline → injection",
        "Each line is one recording. Left: melittin only. Right: melittin + morphine. Y = population mean EXTRACT activity (8th-pct dF/F-like).",
        f"Melittin only: {n_only_neg}/8 go down. Morphine: {n_morph_neg}/4. The thick line is the cohort mean.",
        change_fig,
    )
    two_fig_slide(
        prs,
        "Melittin  |  per mouse, and paired CS1014",
        "Left: injection-cell up/down vs BL session mean (not the same cells). Right: paired Δ (INJ − BL) for CS1014-1b/2w vs 1014m01b/02w.",
        "Melittin-only Δ is negative in every plane. In 3/4 paired CS1014 mouse×mirror lines, morphine is less suppressed or reversed.",
        mouse_fig,
        avg / "fig_takehome_paired.png",
    )
    two_fig_slide(
        prs,
        "Melittin  |  FOV maps  (same red/blue language)",
        "EXTRACT footprints × amp. Red = up vs BL mean, blue = down. P left, A right. Cells are not matched to baseline.",
        "Local mix of up and down, as in face poke. Do not read this as a poke map.",
        avg / "fig4_spatial_response_melittin_only.png" if (avg / "fig4_spatial_response_melittin_only.png").is_file() else avg / "fig4_spatial_response.png",
        avg / "fig4_spatial_response_morphine.png" if (avg / "fig4_spatial_response_morphine.png").is_file() else avg / "fig4_spatial_response.png",
    )

    bullets(
        prs,
        "What this does and does not say",
        [
            "Face poke: every mouse has ~7–13% excited and ~23–30% inhibited poke-locked cells. No Allen/CCF areas yet — only AP on the skull.",
            "Face poke cannot split pain vs touch vs surprise. There is no drug on those movies.",
            "Melittin is a different assay. Baseline and injection EXTRACT IDs are independent, so we cannot follow the same cell into morphine.",
            "Session-level: melittin alone lowers population activity in 7/8 recordings (one CS1014-2w plane goes up after matching durations). Morphine pairing: 3/4 mouse×mirror lines less suppressed or reversed. Small n.",
            "Please confirm which paradigm to take forward for SBI (face poke vs melittin). Imaging is paused until then.",
        ],
    )

    try:
        prs.save(str(OUT))
        print("Wrote", OUT, "slides", len(prs.slides))
    except PermissionError:
        alt = OUT_DIR / "CrystalSkull_facepoke_vs_melittin_Mark_copy.pptx"
        prs.save(str(alt))
        print("Locked, wrote", alt, "slides", len(prs.slides))


if __name__ == "__main__":
    main()
