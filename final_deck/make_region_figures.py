"""Two more figures: melittin by region, and whether morphine is region-specific.

S6 is built to sit beside S2 (face poke) - same panel order, same colour meanings,
same enrichment recipe - so the two paradigms can be read off one against the other.

S7 answers the question in two layers, because one layer alone is misleading:
  within an animal, does the morphine effect differ by region?   (Cochran's Q)
  do the two animals agree on WHICH regions?                     (cross-animal r)
The first is strongly yes and the second is flatly no, and a slide that showed only
the first would claim a region-specific morphine effect this dataset cannot support.

    python make_region_figures.py
"""
from __future__ import annotations

import collections
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy import stats

from make_final_figures import (AMBER, BLUE, BORDER_CMAP, BOUND, GREY, LABEL, LGREY,
                                MEL_ONLY, NAVY, OUT, RED, TEAL, panel_tag)

HERE = os.path.dirname(os.path.abspath(__file__))
ST = json.load(open(os.path.join(HERE, "region_stats.json")))


# ===================================================================== S6
def slide6():
    fig = plt.figure(figsize=(13.0, 5.5))
    gs = fig.add_gridspec(1, 4, width_ratios=[.92, 1.20, .92, .92], wspace=.50,
                          left=.05, right=.985, top=.75, bottom=.25)

    # --- A  spatial map -----------------------------------------------------
    ax = fig.add_subplot(gs[0])
    ax.imshow(np.where(LABEL > 0, .94, 1.0), cmap="gray", vmin=0, vmax=1)
    ax.imshow(np.ma.masked_where(~BOUND, BOUND), cmap=BORDER_CMAP, alpha=.45,
              interpolation="nearest")
    up = [r for r in MEL_ONLY if r["exc"]]
    dn = [r for r in MEL_ONLY if r["inh"]]
    pts = ([(r["col"], r["row"], RED) for r in up[::14]] +
           [(r["col"], r["row"], BLUE) for r in dn[::31]])
    rng = np.random.default_rng(0)
    pts = [pts[i] for i in rng.permutation(len(pts))]   # draw order, not colour, decides
    ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=1.6,
               c=[p[2] for p in pts], alpha=.45, linewidths=0, rasterized=True)
    ax.set_xlim(60, 400); ax.set_ylim(460, 60)
    ax.set_xticks([]); ax.set_yticks([])
    ax.plot([90, 90 + 2000 / 25], [440, 440], "k-", lw=3)
    ax.text(90 + 1000 / 25, 432, "2 mm", ha="center", fontsize=9)
    ax.annotate("anterior", (.98, .975), xycoords="axes fraction", ha="right",
                va="top", fontsize=8.4, color=GREY)
    ax.annotate("posterior", (.98, .025), xycoords="axes fraction", ha="right",
                va="bottom", fontsize=8.4, color=GREY)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=RED, ms=7,
                              label=f"up  {len(up):,} cells  (1 in 14 drawn)"),
                       Line2D([], [], marker="o", ls="", color=BLUE, ms=7,
                              label=f"down  {len(dn):,} cells  (1 in 31 drawn)")],
              frameon=False, loc="upper center", bbox_to_anchor=(.5, -.015),
              fontsize=8.4, handletextpad=.4)
    ax.set_title("Where the responding cells are", fontsize=10.5,
                 fontweight="bold", pad=7)
    panel_tag(ax, "A", dx=-0.07, dy=1.17)

    # --- B  enrichment ------------------------------------------------------
    ax = fig.add_subplot(gs[1])
    up_e = [d for d in ST["A"]["up"] if d["area"] not in ("(outside)", "SCzo")]
    sel = up_e[:7] + up_e[-4:]
    y = np.arange(len(sel))[::-1]
    ax.barh(y, [d["k"] for d in sel],
            color=[RED if d["k"] > 1 else BLUE for d in sel], height=.72)
    ax.axvline(1, color="k", lw=1.2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{d['area']}  (n={d['n']:,})" for d in sel], fontsize=8.8)
    for yy, d in zip(y, sel):
        q = d["q"]
        star = "***" if q < 1e-3 else ("**" if q < .01 else ("*" if q < .05 else "n.s."))
        ax.text(d["k"] + .03, yy, f"{d['k']:.2f}×  {star}", va="center",
                fontsize=8.3, color="#222222")
    ax.set_xlim(0, 2.25)
    ax.set_xlabel("fold enrichment of cells driven up\n"
                  f"(area rate ÷ cortex-wide rate, {100*ST['A']['base_up']:.1f} %)")
    ax.set_title("Melittin drives posterior-medial parietal cortex\n"
                 "27 areas with ≥ 300 cells; binomial test vs cortex-wide\n"
                 "rate, BH-FDR.  *** q<0.001  ** q<0.01  * q<0.05",
                 fontsize=9.4, pad=7)
    panel_tag(ax, "B", dx=-0.40, dy=1.18)

    # --- C  up and down are one map, not two --------------------------------
    ax = fig.add_subplot(gs[2])
    du = {d["area"]: d["k"] for d in ST["A"]["up"]}
    dd = {d["area"]: d["k"] for d in ST["A"]["down"]}
    sh = sorted(set(du) & set(dd))
    xs = np.array([du[a] for a in sh]); ys = np.array([dd[a] for a in sh])
    ax.scatter(xs, ys, s=38, c=NAVY, zorder=3, edgecolor="white", linewidth=.8)
    for a, x, yv in zip(sh, xs, ys):
        if x > 1.3 or x < .3 or yv > 1.35:
            ax.annotate(a, (x, yv), xytext=(5, 4), textcoords="offset points",
                        fontsize=7.8)
    ax.axhline(1, color=LGREY, lw=1); ax.axvline(1, color=LGREY, lw=1)
    ax.set_xlabel("enrichment of cells driven UP")
    ax.set_ylabel("enrichment of cells driven DOWN")
    ax.set_title(f"'Up' and 'down' are one map, not two\n"
                 f"Spearman ρ = {ST['A']['rho_updown']:+.2f}, "
                 f"p = {ST['A']['p_updown']:.0e}, n = {len(sh)} areas\n"
                 f"{100*ST['A']['base_up']:.0f} % + "
                 f"{100*ST['A']['base_down']:.0f} % = "
                 f"{100*(ST['A']['base_up']+ST['A']['base_down']):.0f} % of cells "
                 f"are called responsive",
                 fontsize=9.2, pad=7)
    panel_tag(ax, "C", dx=-0.30, dy=1.23)

    # --- D  A-P gradient, melittin vs face poke -----------------------------
    ax = fig.add_subplot(gs[3])
    B = ST["B"]
    specs = [("face poke\ninhibited", B["facepoke_inh"], TEAL),
             ("melittin\ndown", B["mel_inh_4"], AMBER)]
    for k, (lab, d, c) in enumerate(specs):
        a = np.array(d["ant"]); p_ = np.array(d["pos"])
        x0, x1 = k * 1.25, k * 1.25 + .55
        for va, vp in zip(a, p_):
            ax.plot([x0, x1], [va, vp], "-", color=LGREY, lw=1.0, zorder=1)
        ax.scatter(np.full_like(a, x0), a, s=26, color=c, zorder=3,
                   edgecolor="white", linewidth=.7)
        ax.scatter(np.full_like(p_, x1), p_, s=26, color=c, zorder=3,
                   edgecolor="white", linewidth=.7)
        ax.plot([x0, x1], [a.mean(), p_.mean()], "-o", color="k", lw=2.2, ms=7,
                zorder=4)
        # the stats go BELOW the axis. Put them inside and they land on the melittin
        # points, which sit at 56-76 % - the whole top of the panel is occupied.
        unit = "rec" if d["n"] > 4 else "mice"
        ax.text((x0 + x1) / 2, -.16,
                f"{lab.replace(chr(10), ' ')}\n"
                f"{d['delta']:+.1f} pp\n"
                f"{d['same']}/{d['n']} {unit}\n"
                f"p = {d['p']:.2g}",
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=7.8, color=c, fontweight="bold")
    ax.set_xticks([0, .55, 1.25, 1.80])
    ax.set_xticklabels(["ant", "post", "ant", "post"], fontsize=8.6)
    ax.set_xlim(-.3, 2.1)
    ax.set_ylim(10, 80)
    ax.set_ylabel("cells in that half (%)")
    ax.set_title("Face poke's posterior gradient\nis not demonstrable in melittin\n"
                 "Wilcoxon signed-rank, each recording\nsplit at its own median AP",
                 fontsize=9.2, pad=7)
    panel_tag(ax, "D", dx=-0.30, dy=1.23)

    n_sess = len({r["mouse"] for r in MEL_ONLY})
    fig.suptitle(f"Slide 6  |  Melittin: {ST['A']['n_cells']:,} cells, "
                 f"{n_sess} mice (8 recordings)  —  "
                 f"{100*ST['A']['base_up']:.1f} % driven up, "
                 f"{100*ST['A']['base_down']:.1f} % driven down",
                 fontsize=13.5, fontweight="bold", y=.965)
    fig.savefig(os.path.join(OUT, "S6_melittin_regions.png"), dpi=170)
    plt.close(fig)
    print("S6 done")


# ===================================================================== S7
def slide7():
    fig = plt.figure(figsize=(13.0, 5.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.22, 1.00, 1.02], wspace=.40,
                          left=.055, right=.995, top=.78, bottom=.21)

    C = ST["C"]; D = ST["D"]
    A1, A2 = "CS1014-1b", "CS1014-2w"
    m1 = {r["area"]: r for r in C[A1]["areas"]}
    m2 = {r["area"]: r for r in C[A2]["areas"]}
    areas = D["areas"]
    order = sorted(areas, key=lambda a: -(m1[a]["lor"] + m2[a]["lor"]) / 2)

    # --- A  forest ----------------------------------------------------------
    ax = fig.add_subplot(gs[0])
    # VISpor in CS1014-2w has a CI running to +/- 4. Letting the axis follow it
    # squashes the other 21 areas into a stripe, so clip and mark what was clipped.
    XL = 1.6
    y = np.arange(len(order))[::-1]
    clipped = 0
    for yy, a in zip(y, order):
        for off, mm, c in ((+.20, m1, TEAL), (-.20, m2, AMBER)):
            r = mm[a]
            lo, hi = r["lor"] - 1.96 * r["se"], r["lor"] + 1.96 * r["se"]
            ax.plot([max(lo, -XL), min(hi, XL)], [yy + off] * 2, "-", color=c,
                    lw=1.6, alpha=.85)
            for v, d in ((lo, -1), (hi, 1)):
                if abs(v) > XL:
                    clipped += 1
                    ax.plot(d * XL, yy + off, marker=("<" if d < 0 else ">"),
                            color=c, ms=4.5)
            if abs(r["lor"]) <= XL:
                ax.plot(r["lor"], yy + off, "o", color=c, ms=4.2)
    ax.axvline(0, color="k", lw=1.2)
    for mm, c, key in ((m1, TEAL, A1), (m2, AMBER, A2)):
        ax.axvline(C[key]["pool"], color=c, ls="--", lw=1.4, alpha=.8)
    ax.set_xlim(-XL, XL)
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=7.6)
    ax.set_xlabel("log odds ratio for being driven up,\n"
                  "melittin + morphine  vs  melittin alone\n"
                  f"dot = estimate,  bar = 95 % CI,  dashed = that animal's pooled\n"
                  f"value.  Arrowheads mark {clipped} CI ends clipped at ±{XL}.",
                  fontsize=8.6)
    ax.set_title(f"Within each animal the morphine effect differs by region\n"
                 f"{A1}:  Q = {C[A1]['Q']:.0f}, df = {C[A1]['df']}, "
                 f"p = {C[A1]['pQ']:.0e}, I² = {C[A1]['I2']:.0f} %\n"
                 f"{A2}:  Q = {C[A2]['Q']:.0f}, df = {C[A2]['df']}, "
                 f"p = {C[A2]['pQ']:.0e}, I² = {C[A2]['I2']:.0f} %",
                 fontsize=9.2, pad=7)
    ax.legend(handles=[Line2D([], [], color=TEAL, lw=3, label=A1),
                       Line2D([], [], color=AMBER, lw=3, label=A2)],
              frameon=True, framealpha=.9, edgecolor="none", fontsize=8.2,
              loc="lower left")
    panel_tag(ax, "A", dx=-0.17, dy=1.18)

    # --- B  do the animals agree ---------------------------------------------
    ax = fig.add_subplot(gs[1])
    x = np.array(D["lor1"]); yv = np.array(D["lor2"])
    lim = 1.1 * max(np.abs(np.r_[x, yv]).max(), .8)
    ax.axhline(0, color=LGREY, lw=1); ax.axvline(0, color=LGREY, lw=1)
    ax.plot([-lim, lim], [-lim, lim], ":", color=LGREY, lw=1.2)
    ax.scatter(x, yv, s=44, c=NAVY, zorder=3, edgecolor="white", linewidth=.8)
    for a, xx, yy in zip(D["areas"], x, yv):
        if abs(xx) > .3 or abs(yy) > .45:
            ax.annotate(a, (xx, yy), xytext=(5, 4), textcoords="offset points",
                        fontsize=7.8)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel(f"{A1}:  per-area log OR")
    ax.set_ylabel(f"{A2}:  per-area log OR")
    rP, pP = D["pearson_lor"]
    ax.set_title(f"...but the two animals do not agree on WHICH regions\n"
                 f"Pearson r = {rP:+.2f}, p = {pP:.2f}, n = {D['n']} areas\n"
                 f"same direction in {D['agree']}/{D['n']} "
                 f"(binomial vs chance, p = {D['p_binom']:.2f})",
                 fontsize=9.2, pad=7)
    ax.text(.03, .03, "dotted line = the two animals agree exactly",
            transform=ax.transAxes, va="bottom", fontsize=7.8, color=GREY)
    panel_tag(ax, "B", dx=-0.25, dy=1.18)

    # --- C  what survives ----------------------------------------------------
    ax = fig.add_subplot(gs[2])
    ax.axis("off")

    sub = ax.inset_axes([.16, .76, .74, .22])
    for k, key in enumerate((A1, A2)):
        p, se = C[key]["pool"], C[key]["se_pool"]
        sub.errorbar(np.exp(p), k, xerr=[[np.exp(p) - np.exp(p - 1.96 * se)],
                                         [np.exp(p + 1.96 * se) - np.exp(p)]],
                     fmt="o", color=(TEAL if k == 0 else AMBER), ms=8, capsize=4,
                     lw=2)
        sub.text(np.exp(p), k + .34, f"OR {np.exp(p):.2f}", ha="center", fontsize=8.4)
    sub.axvline(1, color="k", ls="--", lw=1.2)
    sub.set_yticks([0, 1]); sub.set_yticklabels([A1, A2], fontsize=8.4)
    sub.set_ylim(-.6, 1.8)
    sub.set_xlim(.95, 1.30)
    sub.set_xlabel("pooled odds ratio,\nmorphine vs melittin", fontsize=8.2)
    sub.set_title("What IS consistent: the overall shift", fontsize=9.2,
                  fontweight="bold")
    sub.tick_params(labelsize=8)

    E1, E2 = ST["E"][A1], ST["E"][A2]
    txt = (
        "READING A AND B TOGETHER\n"
        "\n"
        "A   Inside one animal the morphine effect is far\n"
        "      from uniform. I² of 82 % and 92 % means most\n"
        "      of the spread between areas is real, not\n"
        "      sampling noise.\n"
        "\n"
        "B   The two animals rank the areas completely\n"
        f"      differently: r = {D['pearson_lor'][0]:+.2f}, "
        f"{D['agree']}/{D['n']} directions agree,\n"
        f"      CIs overlap in {D['ci_overlap']}/{D['n']}. At n = 2 the region\n"
        "      pattern does not reproduce, so no single area\n"
        "      can be named on a slide.\n"
        "\n"
        "What survives both is the whole-cortex shift:\n"
        f"OR {np.exp(C[A1]['pool']):.2f} and {np.exp(C[A2]['pool']):.2f}, "
        "both above 1, same direction.\n"
        "\n"
        "ALSO TESTED, ALSO NOT REPRODUCIBLE\n"
        "Does morphine act most where melittin acts most?\n"
        f"   {A1}  r = {E1['r']:+.2f}  beyond the noise null\n"
        f"   {A2}  r = {E2['r']:+.2f}  inside it — one animal"
    )
    ax.text(-.02, .57, txt, transform=ax.transAxes, va="top", ha="left",
            fontsize=8.1, color="#222222", linespacing=1.42)
    panel_tag(ax, "C", dx=-0.05, dy=1.10)

    fig.suptitle("Slide 7  |  Is morphine region-specific?  Within an animal yes, "
                 "across animals not reproducible at n = 2",
                 fontsize=13.5, fontweight="bold", y=.965)
    fig.savefig(os.path.join(OUT, "S7_morphine_regions.png"), dpi=170)
    plt.close(fig)
    print("S7 done")


if __name__ == "__main__":
    slide6(); slide7()
    print("->", OUT)
