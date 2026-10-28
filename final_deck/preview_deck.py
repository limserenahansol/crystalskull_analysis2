"""Render a layout proof of the deck so overflow can actually be seen.

There is no LibreOffice on this machine, so the .pptx cannot be rasterised directly.
This walks the shapes with python-pptx and redraws each slide in matplotlib at the
real geometry and the real font sizes. It is not pixel-identical to PowerPoint - it
exists to answer one question: does any text run past its box or off the slide?

Boxes that overflow are outlined in red.

    python preview_deck.py   ->  final_deck/preview/slide_NN.png
"""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pptx import Presentation
from pptx.util import Emu

HERE = os.path.dirname(os.path.abspath(__file__))
PPTX = os.path.join(HERE, "CrystalSkull_summary_5slides_for_Mark.pptx")
OUT = os.path.join(HERE, "preview")
os.makedirs(OUT, exist_ok=True)

DPI = 110
prs = Presentation(PPTX)
SW, SH = Emu(prs.slide_width).inches, Emu(prs.slide_height).inches

# PowerPoint wraps at the box width; matplotlib does not wrap for us, so wrap by
# hand using an estimated character width. 0.50 em is close for Calibri at these
# sizes - deliberately a little wide, so the proof errs toward reporting overflow.
EM = 0.50


def wrap(line, width_in, pt):
    cpl = max(int(width_in / (pt / 72.0 * EM)), 8)
    if len(line) <= cpl:
        return [line]
    out, cur = [], ""
    for w in line.split(" "):
        if len(cur) + len(w) + 1 <= cpl:
            cur = f"{cur} {w}".strip()
        else:
            out.append(cur)
            cur = w
    if cur:
        out.append(cur)
    return out


for idx, slide in enumerate(prs.slides, start=1):
    fig = plt.figure(figsize=(SW, SH), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, SW); ax.set_ylim(SH, 0)
    ax.axis("off")
    ax.add_patch(Rectangle((0, 0), SW, SH, fc="white", ec="#DDDDDD"))
    problems = []

    for sh in slide.shapes:
        L, T = Emu(sh.left).inches, Emu(sh.top).inches
        Wd, Ht = Emu(sh.width).inches, Emu(sh.height).inches

        if sh.shape_type == 13 or sh.has_text_frame is False:   # picture
            ax.add_patch(Rectangle((L, T), Wd, Ht, fc="#EEF2F7", ec="#8899AA"))
            ax.text(L + Wd / 2, T + Ht / 2, "[figure]", ha="center", va="center",
                    color="#667788", fontsize=11)
            if L < 0 or T < 0 or L + Wd > SW or T + Ht > SH:
                problems.append(f"picture runs off the slide")
            continue

        y = T + 0.06
        used = 0.06
        for p in sh.text_frame.paragraphs:
            txt = "".join(r.text for r in p.runs)
            pt = 12.0
            bold = False
            for r in p.runs:
                if r.font.size:
                    pt = r.font.size.pt
                bold = bold or bool(r.font.bold)
            lh = pt / 72.0 * 1.22 + (p.space_after.pt / 72.0 if p.space_after else 0)
            for seg in (wrap(txt, Wd, pt) if txt else [""]):
                ax.text(L + 0.04, y, seg, ha="left", va="top", fontsize=pt * .92,
                        fontweight="bold" if bold else "normal", color="#222222")
                y += lh
                used += lh
        if y > SH - 0.03:
            problems.append(f"text box at top={T:.2f}in overflows the slide bottom "
                            f"(ends at {y:.2f} of {SH:.2f})")
            ax.add_patch(Rectangle((L, T), Wd, max(used, Ht), fc="none",
                                   ec="red", lw=2))
        elif used > Ht + 0.25:
            problems.append(f"text box at top={T:.2f}in is taller than declared "
                            f"({used:.2f} vs {Ht:.2f} in) - fine if nothing is below")
            ax.add_patch(Rectangle((L, T), Wd, used, fc="none", ec="orange", lw=1.5))
        else:
            ax.add_patch(Rectangle((L, T), Wd, Ht, fc="none", ec="#CCCCCC", lw=.8))

    fig.savefig(os.path.join(OUT, f"slide_{idx:02d}.png"), dpi=DPI)
    plt.close(fig)
    flag = "  <-- CHECK" if any("overflows" in p for p in problems) else ""
    print(f"slide {idx:2d}{flag}")
    for p in problems:
        print(f"     {p}")

print(f"\n-> {OUT}")
