"""Merged recordings + Allen-area statistics for the two assays.

Face poke: same cells, peri-poke window vs pre-poke baseline (signed-rank p<0.05).
Melittin: injection cells vs baseline-session mean (independent EXTRACT IDs).
Do not treat the two percentages as the same measurement.

Inference for "which areas change more":
  recording-level Wilcoxon signed-rank of (area rate - that recording's overall rate),
  BH-FDR across areas. Cell-level Fisher is exploratory (cells are nested).
"""
from __future__ import annotations

import collections
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.stats import fisher_exact, mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures_v2")
os.makedirs(OUT, exist_ok=True)

CTX = {
    "MOs", "MOp", "SSp-m", "SSp-n", "SSp-bfd", "SSp-ul", "SSp-ll", "SSp-tr", "SSp-un",
    "SSs", "VISp", "VISa", "VISam", "VISal", "VISpm", "VISrl", "VISl", "VISpor",
    "RSPd", "RSPv", "RSPagl", "AUDp", "AUDd", "AUDv", "ACAd", "TEa",
}
MIN_POOLED = 200
MIN_PER_REC = 40
MIN_RECS = 8
RED, BLUE, GREY = "#C0392B", "#2C6FBB", "#555555"
COND_COL = {
    "Face poke": "#7B4FA8",
    "Melittin": "#0D5C63",
    "Melittin+morphine": "#D08C1E",
}


def poke_animal(movie: str) -> str:
    if movie.startswith("0004749770_m1_") or "d250813" in movie and "_m1_" in movie[:30]:
        return "m1"
    if movie.startswith("m200_"):
        return "m200"
    if movie.startswith("m201_"):
        return "m201"
    if movie.startswith("m3_"):
        return "m3"
    return movie.split("_")[0]


def load_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def as_bool(v) -> bool:
    return str(v).strip().lower() in {"1", "true", "yes"}


def load_atlas():
    A = sio.loadmat(os.path.join(HERE, "allen_dorsal_map.mat"))
    label = A["label_img"]
    bound = A["boundary"].astype(bool)
    bregma = A["bregma_rc"].ravel() + 1
    acr = [str(x[0]) for x in A["acronyms"].ravel()]
    ids = A["ids"].ravel()
    id2acr = {int(i): a for i, a in zip(ids, acr)}
    ctx_mask = np.isin(label, [i for i, a in id2acr.items() if a in CTX])
    return label, bound, bregma, ctx_mask


def attach_meta(poke, mel):
    for r in poke:
        r["_cond"] = "Face poke"
        r["animal"] = poke_animal(r["movie"])
        r["is_exc"] = as_bool(r["is_exc"])
        r["is_inh"] = as_bool(r["is_inh"])
        r["amp"] = float(r["amp"])
    for r in mel:
        r["_cond"] = "Melittin+morphine" if r["cohort"] == "melittin_morphine" else "Melittin"
        r["animal"] = r["mouse"]
        r["is_exc"] = as_bool(r["is_exc"])
        r["is_inh"] = as_bool(r["is_inh"])
        r["amp"] = float(r["amp_pct"]) if r.get("amp_pct") not in (None, "") else float(r["amp"])
    return poke, mel


def summarize(rows, name):
    n = len(rows)
    n_exc = sum(r["is_exc"] for r in rows)
    n_inh = sum(r["is_inh"] for r in rows)
    animals = sorted({r["animal"] for r in rows})
    recs = sorted({r["label"] for r in rows})
    print(f"{name:22s} {n:7d} cells  {n_exc:6d} up ({100*n_exc/n:.1f}%)  "
          f"{n_inh:6d} down ({100*n_inh/n:.1f}%)  "
          f"{len(recs)} recs  {len(animals)} animals {animals}")


def region_counts(rows):
    tot = collections.Counter()
    exc = collections.Counter()
    inh = collections.Counter()
    amp = collections.defaultdict(list)
    for r in rows:
        a = r["area_acronym"]
        if a not in CTX:
            continue
        tot[a] += 1
        if r["is_exc"]:
            exc[a] += 1
        if r["is_inh"]:
            inh[a] += 1
        amp[a].append(r["amp"])
    return tot, exc, inh, amp


def recording_deltas(rows, kind="exc"):
    """Per recording, area_rate - recording_rate, for areas with enough cells."""
    by = collections.defaultdict(list)
    recs = collections.defaultdict(list)
    for r in rows:
        if r["area_acronym"] in CTX:
            recs[r["label"]].append(r)
    for lab, rec_rows in recs.items():
        if kind == "exc":
            rec_rate = np.mean([r["is_exc"] for r in rec_rows])
        else:
            rec_rate = np.mean([r["is_inh"] for r in rec_rows])
        by_area = collections.defaultdict(list)
        for r in rec_rows:
            by_area[r["area_acronym"]].append(r)
        for area, ar in by_area.items():
            if len(ar) < MIN_PER_REC:
                continue
            if kind == "exc":
                rate = np.mean([x["is_exc"] for x in ar])
            else:
                rate = np.mean([x["is_inh"] for x in ar])
            by[area].append(rate - rec_rate)
    return by


def fisher_area_vs_rest(tot, hit, area):
    a_yes = hit.get(area, 0)
    a_no = tot[area] - a_yes
    rest_yes = sum(hit.values()) - a_yes
    rest_no = (sum(tot.values()) - tot[area]) - rest_yes
    if min(a_yes + a_no, rest_yes + rest_no) == 0:
        return np.nan
    table = [[a_yes, a_no], [rest_yes, rest_no]]
    try:
        _, p = fisher_exact(table, alternative="two-sided")
    except ValueError:
        p = np.nan
    return p


def fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    if ok.sum() == 0:
        return q
    _, q_ok, _, _ = multipletests(p[ok], method="fdr_bh")
    q[ok] = q_ok
    return q


def test_regions(rows, cond_name):
    tot, exc, inh, amp = region_counts(rows)
    n_all = sum(tot.values())
    p_exc = sum(exc.values()) / max(n_all, 1)
    p_inh = sum(inh.values()) / max(n_all, 1)
    d_exc = recording_deltas(rows, "exc")
    d_inh = recording_deltas(rows, "inh")
    areas = sorted(a for a, n in tot.items() if n >= MIN_POOLED)
    recs = []
    for a in areas:
        n = tot[a]
        e = exc.get(a, 0)
        i = inh.get(a, 0)
        dx = d_exc.get(a, [])
        di = d_inh.get(a, [])
        p_w_e = np.nan
        p_w_i = np.nan
        if len(dx) >= MIN_RECS:
            try:
                p_w_e = wilcoxon(dx).pvalue
            except ValueError:
                p_w_e = np.nan
        if len(di) >= MIN_RECS:
            try:
                p_w_i = wilcoxon(di).pvalue
            except ValueError:
                p_w_i = np.nan
        recs.append({
            "condition": cond_name,
            "area": a,
            "n": n,
            "n_exc": e,
            "n_inh": i,
            "pct_exc": 100 * e / n,
            "pct_inh": 100 * i / n,
            "enr_exc": (e / n) / p_exc if p_exc > 0 else np.nan,
            "enr_inh": (i / n) / p_inh if p_inh > 0 else np.nan,
            "mean_amp": float(np.mean(amp[a])),
            "n_recs_exc": len(dx),
            "n_recs_inh": len(di),
            "delta_exc": float(np.mean(dx)) if dx else np.nan,
            "delta_inh": float(np.mean(di)) if di else np.nan,
            "p_wilcox_exc": p_w_e,
            "p_wilcox_inh": p_w_i,
            "p_fisher_exc": fisher_area_vs_rest(tot, exc, a),
            "p_fisher_inh": fisher_area_vs_rest(tot, inh, a),
        })
    qe = fdr([r["p_wilcox_exc"] for r in recs])
    qi = fdr([r["p_wilcox_inh"] for r in recs])
    qfe = fdr([r["p_fisher_exc"] for r in recs])
    qfi = fdr([r["p_fisher_inh"] for r in recs])
    for r, a, b, c, d in zip(recs, qe, qi, qfe, qfi):
        r["q_wilcox_exc"] = a
        r["q_wilcox_inh"] = b
        r["q_fisher_exc"] = c
        r["q_fisher_inh"] = d
        r["sig_exc"] = bool(np.isfinite(a) and a < 0.05 and r["enr_exc"] > 1)
        r["sig_inh"] = bool(np.isfinite(b) and b < 0.05 and r["enr_inh"] > 1)
    return recs, p_exc, p_inh


def mouse_consistency(rows, area, kind="exc"):
    by = collections.defaultdict(list)
    for r in rows:
        if r["area_acronym"] == area:
            by[r["animal"]].append(r["is_exc"] if kind == "exc" else r["is_inh"])
    overall = np.mean([r["is_exc"] if kind == "exc" else r["is_inh"] for r in rows])
    n_above = 0
    n_mice = 0
    for vals in by.values():
        if len(vals) < 30:
            continue
        n_mice += 1
        if np.mean(vals) > overall:
            n_above += 1
    return n_above, n_mice


def write_table(path, rows, extra_fields=None):
    if not rows:
        return
    fields = list(rows[0].keys())
    if extra_fields:
        for f in extra_fields:
            if f not in fields:
                fields.append(f)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print("wrote", path)


def stars(q):
    if not np.isfinite(q):
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ""


def response_map(rows, key, label, sigma_px=6):
    h, w = label.shape
    tot = np.zeros((h, w))
    sel = np.zeros((h, w))
    for r in rows:
        y = int(round(float(r["atlas_row"])))
        x = int(round(float(r["atlas_col"])))
        if 0 <= y < h and 0 <= x < w:
            tot[y, x] += 1
            if r[key]:
                sel[y, x] += 1
    from scipy.ndimage import gaussian_filter
    tot_s = gaussian_filter(tot, sigma_px)
    sel_s = gaussian_filter(sel, sigma_px)
    out = np.full((h, w), np.nan)
    if (tot_s > 0).any():
        ok = tot_s > np.percentile(tot_s[tot_s > 0], 25)
        out[ok] = 100 * sel_s[ok] / tot_s[ok]
    return out


def draw_atlas(ax, img, cmap, bound, bregma, ctx_mask, vmin=None, vmax=None):
    ax.imshow(np.where(ctx_mask, 0.93, np.nan), cmap="gray", vmin=0, vmax=1)
    im = ax.imshow(np.where(ctx_mask, img, np.nan), cmap=cmap, vmin=vmin, vmax=vmax)
    ax.imshow(np.ma.masked_where(~bound, bound), cmap="gray_r", alpha=0.35, interpolation="nearest")
    ax.plot(bregma[1], bregma[0], "k+", ms=8, mew=1.4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(60, 420)
    ax.set_ylim(400, 60)
    return im


def fig_ccf_maps(all_rows, label, bound, bregma, ctx_mask):
    conds = ["Face poke", "Melittin", "Melittin+morphine"]
    fig, axes = plt.subplots(2, 3, figsize=(15.2, 10.6))
    for j, cond in enumerate(conds):
        sub = [r for r in all_rows if r["_cond"] == cond]
        for i, (key, cmap, ttl) in enumerate([
            ("is_exc", "Reds", "Excited / up vs baseline"),
            ("is_inh", "Blues", "Inhibited / down vs baseline"),
        ]):
            ax = axes[i, j]
            m = response_map(sub, key, label)
            finite = m[np.isfinite(m)]
            vmax = float(np.percentile(finite, 98)) if finite.size else 1.0
            im = draw_atlas(ax, m, cmap, bound, bregma, ctx_mask, vmin=0, vmax=max(vmax, 5))
            n_an = len({r["animal"] for r in sub})
            ax.set_title(f"{cond}\n{ttl}\n{len(sub):,} cells, {n_an} animals",
                         fontsize=10, color=COND_COL[cond], fontweight="bold")
            cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
            cb.set_label("% of local cells", fontsize=8)
            cb.ax.tick_params(labelsize=7)
    fig.suptitle(
        "Merged recordings on Allen CCF   ·   anterior up, + = bregma\n"
        "Face poke = peri-poke vs pre-poke (same cells).  "
        "Melittin = injection cell vs baseline-session mean (IDs not matched).",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    path = os.path.join(OUT, "fig_merged_ccf_maps.png")
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print("wrote", path)


def fig_poke_bars(stats):
    order = sorted(stats, key=lambda r: r["enr_exc"])
    y = np.arange(len(order))
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 7.2), sharey=True)
    for ax, kind, col, title in (
        (axes[0], "exc", RED, "Excited vs pre-poke baseline"),
        (axes[1], "inh", BLUE, "Inhibited vs pre-poke baseline"),
    ):
        vals = [r[f"enr_{kind}"] for r in order]
        q = [r[f"q_wilcox_{kind}"] for r in order]
        ax.barh(y, [v - 1 for v in vals], left=1, color=col, edgecolor="white", height=0.72)
        for i, (v, qi) in enumerate(zip(vals, q)):
            ax.text(v + 0.03, i, f"{v:.2f}{stars(qi)}", va="center", fontsize=8, color=GREY)
        ax.axvline(1, color="k", lw=1)
        ax.set_xlabel("enrichment  (area %  /  all-cortex %)")
        ax.set_title(title, fontweight="bold", color=col)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([r["area"] for r in order], fontsize=9)
    fig.suptitle(
        "Face poke  ·  all 18 recordings merged  ·  which Allen areas change more\n"
        "1.0 = same rate as the rest of cortex.  * FDR q<0.05 from recording-level Wilcoxon.",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    path = os.path.join(OUT, "fig_poke_region_enrichment.png")
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print("wrote", path)


def fig_mel_bars(mel_s, mor_s):
    by_m = {r["area"]: r for r in mel_s}
    by_o = {r["area"]: r for r in mor_s}
    areas = sorted(
        set(by_m) | set(by_o),
        key=lambda a: (by_m[a]["pct_exc"] if a in by_m else by_o[a]["pct_exc"]),
    )
    y = np.arange(len(areas))
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 7.4), sharey=True)
    for ax, kind, title in (
        (axes[0], "exc", "% up vs baseline session"),
        (axes[1], "inh", "% down vs baseline session"),
    ):
        xm = [by_m[a][f"pct_{kind}"] if a in by_m else np.nan for a in areas]
        xo = [by_o[a][f"pct_{kind}"] if a in by_o else np.nan for a in areas]
        ax.barh(y - 0.18, xm, height=0.34, color=COND_COL["Melittin"], label="Melittin")
        ax.barh(y + 0.18, xo, height=0.34, color=COND_COL["Melittin+morphine"], label="Melittin+morphine")
        ax.set_xlabel(title)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, loc="lower right")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(areas, fontsize=9)
    fig.suptitle(
        "Melittin vs melittin+morphine  ·  injection cells vs that day's baseline\n"
        "Same measurement in both arms. Cells are not tracked from baseline into injection.",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    path = os.path.join(OUT, "fig_mel_vs_morphine_regions.png")
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print("wrote", path)


def fig_sig_table(poke_s, mel_s, mor_s):
    fig, ax = plt.subplots(figsize=(16.5, 8.4))
    ax.axis("off")
    rows = [["Assay", "Direction", "Allen areas with more change", "n recs"]]

    def line(stats, cond, kind, label):
        hits = [r for r in stats if r[f"sig_{kind}"]]
        hits = sorted(hits, key=lambda r: -r[f"enr_{kind}"])
        if hits:
            txt = ",  ".join(f"{r['area']} ({r[f'enr_{kind}']:.2f}×)" for r in hits[:8])
        else:
            # Reporting a bare "none after FDR" throws away the ranking, which is the part
            # that is actually informative when n is small. List the most enriched areas
            # instead, with their q values, and label them as not surviving FDR.
            top = sorted((r for r in stats
                          if np.isfinite(r[f"enr_{kind}"]) and r[f"enr_{kind}"] > 1),
                         key=lambda r: -r[f"enr_{kind}"])[:5]
            tested = any(np.isfinite(r[f"q_wilcox_{kind}"]) for r in stats)
            if not top:
                txt = "no area enriched above the cortex-wide rate"
            elif tested:
                lst = ",  ".join(
                    f"{r['area']} ({r[f'enr_{kind}']:.2f}×, q={r[f'q_wilcox_{kind}']:.2f})"
                    for r in top)
                txt = "none reach q<0.05 — most enriched:\n" + lst
            else:
                # With 4 recordings the signed-rank test was never run (MIN_RECS = 8), and
                # could not have helped anyway: at n = 4 its smallest possible two-sided p
                # is 0.125, so no effect size can reach q < 0.05. Saying "none after FDR"
                # would imply a test happened and failed.
                lst = ",  ".join(f"{r['area']} ({r[f'enr_{kind}']:.2f}×)" for r in top)
                txt = ("too few recordings to test (signed-rank floor p=0.125 at n=4)\n"
                       "ranked by enrichment: " + lst)
        n = int(np.nanmedian([r[f"n_recs_{kind}"] for r in stats]))
        rows.append([cond, label, txt, str(n)])

    line(poke_s, "Face poke", "exc", "more excited")
    line(poke_s, "Face poke", "inh", "more inhibited")
    line(mel_s, "Melittin", "exc", "more up")
    line(mel_s, "Melittin", "inh", "more down")
    line(mor_s, "Melittin+morphine", "exc", "more up")
    line(mor_s, "Melittin+morphine", "inh", "more down")

    table = ax.table(cellText=rows, loc="center", cellLoc="left",
                     colWidths=[0.16, 0.14, 0.60, 0.08])
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1.0, 2.6)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#DDDDDD")
        if r == 0:
            cell.set_facecolor("#1F3A5F")
            cell.set_text_props(color="white", weight="bold")
        elif r in (1, 2):
            cell.set_facecolor("#F4EEF8")
        elif r in (3, 4):
            cell.set_facecolor("#E8F2F2")
        else:
            cell.set_facecolor("#FFF6E8")
    ax.set_title(
        "Which Allen areas show more change than the rest of cortex\n"
        "Enrichment in parentheses. Wilcoxon on recordings, BH-FDR within each assay x direction.\n"
        "Where nothing survives FDR the most enriched areas are still listed, with their q values.",
        fontsize=13, fontweight="bold", pad=12,
    )
    fig.tight_layout()
    path = os.path.join(OUT, "fig_which_regions_table.png")
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print("wrote", path)


def compare_mel_morph(mel_rows, mor_rows, areas):
    out = []
    for a in areas:
        xm, xo = [], []
        for lab, rec in _by_rec(mel_rows).items():
            ar = [r for r in rec if r["area_acronym"] == a]
            if len(ar) >= MIN_PER_REC:
                xm.append(100 * np.mean([r["is_exc"] for r in ar]))
        for lab, rec in _by_rec(mor_rows).items():
            ar = [r for r in rec if r["area_acronym"] == a]
            if len(ar) >= MIN_PER_REC:
                xo.append(100 * np.mean([r["is_exc"] for r in ar]))
        p = np.nan
        if len(xm) >= 3 and len(xo) >= 2:
            try:
                p = mannwhitneyu(xm, xo, alternative="two-sided").pvalue
            except ValueError:
                p = np.nan
        out.append({
            "area": a,
            "mel_pct_exc": float(np.mean(xm)) if xm else np.nan,
            "mor_pct_exc": float(np.mean(xo)) if xo else np.nan,
            "delta": (float(np.mean(xo)) - float(np.mean(xm))) if xm and xo else np.nan,
            "n_recs_mel": len(xm),
            "n_recs_mor": len(xo),
            "p_mwu": p,
        })
    q = fdr([r["p_mwu"] for r in out])
    for r, qi in zip(out, q):
        r["q_mwu"] = qi
    out.sort(key=lambda r: -(r["delta"] if np.isfinite(r["delta"]) else -999))
    return out


def _by_rec(rows):
    d = collections.defaultdict(list)
    for r in rows:
        d[r["label"]].append(r)
    return d


def fig_mel_delta(comp):
    fig, ax = plt.subplots(figsize=(9.2, 6.8))
    y = np.arange(len(comp))
    delta = [r["delta"] if np.isfinite(r["delta"]) else 0 for r in comp]
    cols = [RED if d > 0 else BLUE for d in delta]
    ax.barh(y, delta, color=cols, edgecolor="white")
    for i, r in enumerate(comp):
        mark = stars(r["q_mwu"])
        if mark:
            ax.text(r["delta"] + (0.35 if r["delta"] > 0 else -0.35), i, mark,
                    va="center", ha="left" if r["delta"] > 0 else "right", fontsize=12)
    ax.set_yticks(y)
    ax.set_yticklabels([r["area"] for r in comp], fontsize=9)
    ax.axvline(0, color="k", lw=1)
    ax.set_xlabel("% up  (melittin+morphine)  minus  % up  (melittin)\n"
                  "recording-level means; * FDR q<0.05, Mann–Whitney 8 vs 4 recordings")
    ax.set_title("Does morphine change the regional melittin response?",
                 fontsize=13, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    path = os.path.join(OUT, "fig_morphine_delta_by_region.png")
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print("wrote", path)


def main():
    poke = load_csv(os.path.join(HERE, "n18", "cell_ccf_all.csv"))
    mel = load_csv(os.path.join(HERE, "melittin_ccf", "cell_ccf_all.csv"))
    poke, mel = attach_meta(poke, mel)
    all_rows = poke + mel
    summarize(poke, "Face poke")
    summarize([r for r in mel if r["_cond"] == "Melittin"], "Melittin")
    summarize([r for r in mel if r["_cond"] == "Melittin+morphine"], "Melittin+morphine")

    poke_s, pe, pi = test_regions(poke, "Face poke")
    mel_rows = [r for r in mel if r["_cond"] == "Melittin"]
    mor_rows = [r for r in mel if r["_cond"] == "Melittin+morphine"]
    mel_s, _, _ = test_regions(mel_rows, "Melittin")
    mor_s, _, _ = test_regions(mor_rows, "Melittin+morphine")

    for r in poke_s:
        a, n = mouse_consistency(poke, r["area"], "exc")
        r["mice_exc_above"] = f"{a}/{n}"
        a, n = mouse_consistency(poke, r["area"], "inh")
        r["mice_inh_above"] = f"{a}/{n}"

    write_table(os.path.join(OUT, "poke_region_stats.csv"), poke_s)
    write_table(os.path.join(OUT, "melittin_region_stats.csv"), mel_s)
    write_table(os.path.join(OUT, "morphine_region_stats.csv"), mor_s)

    areas = sorted({r["area"] for r in mel_s} | {r["area"] for r in mor_s})
    comp = compare_mel_morph(mel_rows, mor_rows, areas)
    write_table(os.path.join(OUT, "melittin_vs_morphine_by_region.csv"), comp)

    label, bound, bregma, ctx_mask = load_atlas()
    fig_ccf_maps(all_rows, label, bound, bregma, ctx_mask)
    fig_poke_bars(poke_s)
    fig_mel_bars(mel_s, mor_s)
    fig_sig_table(poke_s, mel_s, mor_s)
    fig_mel_delta(comp)

    print("\nFACE POKE  overall {:.1f}% exc / {:.1f}% inh".format(100 * pe, 100 * pi))
    print("  more excited (FDR):",
          ", ".join(r["area"] for r in sorted(poke_s, key=lambda x: -x["enr_exc"]) if r["sig_exc"]) or "none")
    print("  more inhibited (FDR):",
          ", ".join(r["area"] for r in sorted(poke_s, key=lambda x: -x["enr_inh"]) if r["sig_inh"]) or "none")
    print("MELITTIN more up (FDR):",
          ", ".join(r["area"] for r in sorted(mel_s, key=lambda x: -x["enr_exc"]) if r["sig_exc"]) or "none")
    print("MORPHINE more up (FDR):",
          ", ".join(r["area"] for r in sorted(mor_s, key=lambda x: -x["enr_exc"]) if r["sig_exc"]) or "none")
    n_sig = sum(1 for r in comp if np.isfinite(r["q_mwu"]) and r["q_mwu"] < 0.05)
    print(f"morphine vs melittin areas FDR-significant: {n_sig}/{len(comp)}")
    print("->", OUT)


if __name__ == "__main__":
    main()
