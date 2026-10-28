"""Comparison figures: face poke vs melittin vs melittin+morphine (same axes, same colours).

Reads the face-poke pooled summary CSV and the melittin per-recording results.mat files,
and writes the figure set used by build_combined_paradigm_pptx.py.
"""
import csv
import collections
import glob
import statistics as st
from pathlib import Path

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

ROOT = Path(r"C:\Users\hsollim\behavior_task\crystal_skull_analysis")
OUT = ROOT / "combined_paradigm"
OUT.mkdir(exist_ok=True)

RED, BLUE, GREY = "#C0392B", "#2C6FBB", "#555555"
COL = {"Face poke": "#7B4FA8", "Melittin": "#0D5C63", "Melittin+morphine": "#D08C1E"}
rng = np.random.default_rng(0)

# ---------------------------------------------------------------- face poke
def fp_mouse(movie):
    return "m1" if movie.startswith("0004749770_m1_") else movie.split("_")[0]


fp = []
src = ROOT / "stim_triggered_figures/combined_0813_0818/pooled_stim_summary.csv"
for r in csv.DictReader(open(src)):
    n = int(r["n_neurons"])
    fp.append(dict(par="Face poke", mouse=f"{fp_mouse(r['movie'])}-{r['day']}", n=n,
                   exc=100 * float(r["frac_exc"]), inh=100 * float(r["frac_inh"]),
                   n_exc=round(float(r["frac_exc"]) * n), n_inh=round(float(r["frac_inh"]) * n),
                   p=float(r["p_pop"]), d=np.nan, up_amp=np.nan, dn_amp=np.nan))


# ---------------------------------------------------------------- melittin
def mstr(R, k):
    return "".join(chr(c) for c in np.array(R[k]).flatten())


mel = []
for path in sorted(glob.glob(str(ROOT / "melittin_lab_meeting/per_recording/*/results.mat"))):
    with h5py.File(path, "r") as f:
        R = f["R"]
        g = lambda k: float(np.array(R[k]).flatten()[0])
        coh = mstr(R, "cohort")
        par = "Melittin+morphine" if coh == "melittin_morphine" else "Melittin"
        bl, inj = g("mean_dff_bl"), g("mean_dff_inj")
        amp = np.array(R["amp"]).flatten()
        up = np.array(R["is_up"]).flatten().astype(bool)
        dn = np.array(R["is_down"]).flatten().astype(bool)
        mel.append(dict(par=par, mouse=mstr(R, "mouse"), n=int(g("n_neurons_inj")),
                        exc=100 * g("frac_up"), inh=100 * g("frac_down"),
                        n_exc=int(g("n_up")), n_inh=int(g("n_down")),
                        p=g("p_pop"), d=100 * (inj - bl) / bl,
                        up_amp=100 * float(np.mean(amp[up])) / bl,
                        dn_amp=100 * float(np.mean(amp[dn])) / bl))

rows = fp + mel
PARS = ["Face poke", "Melittin", "Melittin+morphine"]
GROUPS = ["Melittin", "Melittin+morphine"]


def sub(par):
    return [r for r in rows if r["par"] == par]


def sd(v):
    return st.stdev(v) if len(v) > 1 else 0.0


def by_mouse(par, key):
    agg = collections.defaultdict(list)
    for r in sub(par):
        agg[r["mouse"]].append(r[key])
    return {m: st.mean(v) for m, v in agg.items()}


# ---------------------------------------------------------------- statistics
_a = [r["d"] for r in sub("Melittin")]
_b = [r["d"] for r in sub("Melittin+morphine")]
P_REC = mannwhitneyu(_a, _b, alternative="two-sided", method="exact").pvalue
_am = list(by_mouse("Melittin", "d").values())
_bm = list(by_mouse("Melittin+morphine", "d").values())
P_MOUSE = mannwhitneyu(_am, _bm, alternative="two-sided", method="exact").pvalue
_ua = [r["up_amp"] for r in sub("Melittin")]
_ub = [r["up_amp"] for r in sub("Melittin+morphine")]
P_UPAMP = mannwhitneyu(_ua, _ub, alternative="two-sided", method="exact").pvalue
_fa = [r["exc"] for r in sub("Melittin")]
_fb = [r["exc"] for r in sub("Melittin+morphine")]
P_FUP_REC = mannwhitneyu(_fa, _fb, alternative="two-sided", method="exact").pvalue
_fam = list(by_mouse("Melittin", "exc").values())
_fbm = list(by_mouse("Melittin+morphine", "exc").values())
P_FUP_MOUSE = mannwhitneyu(_fam, _fbm, alternative="two-sided", method="exact").pvalue

# ============================================ FIG 1 : fractions, 3 paradigms
fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))
ax = axes[0]
x = np.arange(len(PARS))
w = 0.35
for i, (key, col, lab) in enumerate([("exc", RED, "Excited / up"), ("inh", BLUE, "Inhibited / down")]):
    means = [st.mean([r[key] for r in sub(p)]) for p in PARS]
    errs = [sd([r[key] for r in sub(p)]) for p in PARS]
    ax.bar(x + (i - 0.5) * w, means, w, yerr=errs, capsize=4, color=col, label=lab, edgecolor="white")
    for j, p in enumerate(PARS):
        v = [r[key] for r in sub(p)]
        ax.scatter(np.full(len(v), x[j] + (i - 0.5) * w) + rng.uniform(-.06, .06, len(v)),
                   v, s=14, c="k", alpha=.45, zorder=3, linewidths=0)
ax.set_xticks(x)
ax.set_xticklabels([f"{p}\n(n={len(sub(p))} rec)" for p in PARS])
ax.set_ylabel("% of cells (mean +/- SD across recordings)")
ax.set_title("Responsive-cell fractions are paradigm-specific", fontweight="bold")
ax.legend(frameon=False)
ax.spines[["top", "right"]].set_visible(False)

ax = axes[1]
for p in PARS:
    v = sub(p)
    ax.scatter([r["exc"] for r in v], [r["inh"] for r in v], s=70, color=COL[p], label=p,
               edgecolor="white", linewidth=.8, zorder=3)
ax.set_xlabel("% excited / up")
ax.set_ylabel("% inhibited / down")
ax.set_title("Each recording", fontweight="bold")
ax.legend(frameon=False, loc="upper left")
ax.grid(alpha=.25)
ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("Fig 1  |  Excited vs inhibited cells across the three conditions",
             fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .94])
fig.savefig(OUT / "fig1_fractions_three_paradigms.png", dpi=170)
plt.close(fig)

# ============================================ FIG 2 : per mouse
fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), sharey=True)
for ax, p in zip(axes, PARS):
    agg = collections.OrderedDict()
    for r in sub(p):
        a = agg.setdefault(r["mouse"], dict(n=0, e=0, i=0))
        a["n"] += r["n"]
        a["e"] += r["n_exc"]
        a["i"] += r["n_inh"]
    ks = list(agg)
    xx = np.arange(len(ks))
    w = 0.36
    ax.bar(xx - w / 2, [100 * agg[k]["e"] / agg[k]["n"] for k in ks], w, color=RED, edgecolor="white")
    ax.bar(xx + w / 2, [100 * agg[k]["i"] / agg[k]["n"] for k in ks], w, color=BLUE, edgecolor="white")
    ax.set_xticks(xx)
    ax.set_xticklabels([f"{k}\n{agg[k]['n']:,} cells" for k in ks], fontsize=8, rotation=20, ha="right")
    ax.set_title(f"{p}  ({len(ks)} mice)", color=COL[p], fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("% of cells (pooled within mouse)")
axes[0].legend(handles=[plt.Rectangle((0, 0), 1, 1, color=RED), plt.Rectangle((0, 0), 1, 1, color=BLUE)],
               labels=["Excited / up", "Inhibited / down"], frameon=False, fontsize=9)
fig.suptitle("Fig 2  |  Cells excited vs inhibited, per mouse", fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .92])
fig.savefig(OUT / "fig2_per_mouse.png", dpi=170)
plt.close(fig)

# ============================================ FIG 3 : drug effect
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
ax = axes[0]
for j, p in enumerate(GROUPS):
    v = [r["d"] for r in sub(p)]
    ax.bar(j, st.mean(v), 0.55, yerr=sd(v), capsize=5, color=COL[p], edgecolor="white")
    ax.scatter(np.full(len(v), j) + rng.uniform(-.09, .09, len(v)), v, s=42, c="k",
               alpha=.6, zorder=3, linewidths=0)
ax.axhline(0, color="k", lw=1)
ax.set_xticks(range(2))
ax.set_xticklabels([f"{p}\n(n={len(sub(p))} rec)" for p in GROUPS])
ax.set_ylabel("Change in population mean dF/F\n(injection vs own baseline, %)")
ax.set_title(f"Per recording - exact Mann-Whitney p = {P_REC:.2f} (n.s.)", fontweight="bold", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)

ax = axes[1]
for j, p in enumerate(GROUPS):
    for k, v in by_mouse(p, "d").items():
        ax.scatter(j, v, s=110, color=COL[p], edgecolor="k", zorder=3)
        ax.annotate(k, (j, v), textcoords="offset points", xytext=(10, -3), fontsize=8, color=GREY)
ax.axhline(0, color="k", lw=1)
ax.set_xlim(-.5, 1.7)
ax.set_xticks(range(2))
ax.set_xticklabels(GROUPS)
ax.set_ylabel("Change in population mean dF/F (%)")
ax.set_title(f"Per mouse - n = 4 vs 2, exact Mann-Whitney p = {P_MOUSE:.2f} (n.s.)",
             fontweight="bold", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)
fig.suptitle("Fig 3  |  Drug effect - melittin only vs melittin + morphine",
             fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .90])
fig.savefig(OUT / "fig3_drug_effect.png", dpi=170)
plt.close(fig)

# ============================================ FIG 4 : which way does morphine move it?
fig, axes = plt.subplots(1, 2, figsize=(13.2, 4.9))
ax = axes[0]
MK = {}
for j, p in enumerate(GROUPS):
    v = [r["exc"] for r in sub(p)]
    ax.bar(j, st.mean(v), .5, yerr=sd(v), capsize=5, color=COL[p], alpha=.9, edgecolor="white")
    # one marker shape per mouse, so the mirror pairing is visible
    for r in sub(p):
        mk = MK.setdefault(r["mouse"], ["o", "s", "^", "D", "v", "P"][len(MK) % 6])
        ax.scatter(j + rng.uniform(-.11, .11), r["exc"], s=52, marker=mk, c="k",
                   alpha=.75, zorder=3, linewidths=0)
for j, p in enumerate(GROUPS):
    v = list(by_mouse(p, "exc").values())
    ax.scatter(np.full(len(v), j + 0.30), v, s=95, marker="_", c=RED, linewidths=2.5, zorder=4)
ax.set_xticks(range(2))
ax.set_xticklabels([f"Melittin\n{len(sub('Melittin'))} rec / 4 mice",
                    f"Melittin + morphine\n{len(sub('Melittin+morphine'))} rec / 2 mice"])
ax.set_ylabel("% of cells above baseline (UP)")
ax.set_title("More cells sit above baseline with morphine\n"
             + f"{st.mean(_fa):.1f}% -> {st.mean(_fb):.1f}%   "
             + f"p = {P_FUP_REC:.3f} by recording,  p = {P_FUP_MOUSE:.2f} by mouse",
             fontweight="bold", fontsize=11)
ax.spines[["top", "right"]].set_visible(False)
fig.text(0.045, 0.015, "marker shape = mouse; red dash = mouse mean.  The two points per mouse are the two "
         "mirrors of one session, so they are not independent samples.",
         fontsize=8.5, color=GREY, va="bottom")

ax = axes[1]
ax.axis("off")
hdr = f"{'measure':<21}{'melittin':>10}{'+morphine':>11}{'  p rec':>9}{'  p mouse':>10}"
lines = [
    'Does morphine move the melittin response,',
    'and in which direction?',
    "",
    hdr,
    "-" * len(hdr),
    f"{'fraction UP':<21}{st.mean(_fa):9.1f}%{st.mean(_fb):10.1f}%{P_FUP_REC:9.3f}{P_FUP_MOUSE:10.2f}",
    f"{'fraction DOWN':<21}{st.mean([r['inh'] for r in sub('Melittin')]):9.1f}%"
    f"{st.mean([r['inh'] for r in sub('Melittin+morphine')]):10.1f}%"
    f"{mannwhitneyu([r['inh'] for r in sub('Melittin')], [r['inh'] for r in sub('Melittin+morphine')], alternative='two-sided', method='exact').pvalue:9.3f}"
    f"{mannwhitneyu(list(by_mouse('Melittin','inh').values()), list(by_mouse('Melittin+morphine','inh').values()), alternative='two-sided', method='exact').pvalue:10.2f}",
    f"{'UP amplitude':<21}{st.mean(_ua):9.0f}%{st.mean(_ub):10.0f}%{P_UPAMP:9.3f}"
    f"{mannwhitneyu(list(by_mouse('Melittin','up_amp').values()), list(by_mouse('Melittin+morphine','up_amp').values()), alternative='two-sided', method='exact').pvalue:10.2f}",
    f"{'net population':<21}{st.mean(_a):+9.1f}%{st.mean(_b):+10.1f}%{P_REC:9.3f}{P_MOUSE:10.2f}",
    "",
    "Direction: every measure that moves, moves toward",
    "MORE activity with morphine. Melittin alone suppresses",
    "the cortex; with morphine that suppression is gone.",
    "",
    "But only fraction UP reaches p < 0.05, and only when",
    "the two mirrors of a session are counted as separate",
    "samples. At the animal level (n = 4 vs 2) nothing is",
    "significant.",
]
ax.text(0.0, 0.99, "\n".join(lines), va="top", ha="left", fontsize=10.2, color=GREY,
        family="monospace", linespacing=1.42, transform=ax.transAxes)
fig.suptitle("Fig 4  |  Which way does morphine move the melittin response?",
             fontsize=15, fontweight="bold")
fig.tight_layout(rect=[0, 0.05, 1, .90])
fig.savefig(OUT / "fig4_morphine_direction.png", dpi=170)
plt.close(fig)

# ---------------------------------------------------------------- tables
with open(OUT / "combined_summary.csv", "w", newline="", encoding="utf-8") as fh:
    wr = csv.writer(fh)
    wr.writerow(["paradigm", "mouse", "n_cells", "pct_exc", "pct_inh", "n_exc", "n_inh",
                 "p_pop", "rel_dff_change_pct", "up_amp_pctBL", "dn_amp_pctBL"])
    for r in rows:
        f2 = lambda v: "" if (isinstance(v, float) and np.isnan(v)) else f"{v:.2f}"
        wr.writerow([r["par"], r["mouse"], r["n"], f"{r['exc']:.2f}", f"{r['inh']:.2f}",
                     r["n_exc"], r["n_inh"], f"{r['p']:.3e}", f2(r["d"]),
                     f2(r["up_amp"]), f2(r["dn_amp"])])

with open(OUT / "stats.txt", "w", encoding="utf-8") as fh:
    for p in PARS:
        v = sub(p)
        fh.write(f"{p}: n_rec={len(v)} mice={len(set(r['mouse'] for r in v))} "
                 f"cells={sum(r['n'] for r in v)} "
                 f"exc={st.mean([r['exc'] for r in v]):.1f}+/-{sd([r['exc'] for r in v]):.1f}% "
                 f"inh={st.mean([r['inh'] for r in v]):.1f}+/-{sd([r['inh'] for r in v]):.1f}%\n")
    fh.write(f"net dF/F change: mel {st.mean(_a):+.1f}+/-{sd(_a):.1f}%  "
             f"mel+mor {st.mean(_b):+.1f}+/-{sd(_b):.1f}%  p_rec={P_REC:.3f} p_mouse={P_MOUSE:.3f}\n")
    fh.write(f"UP-cell amplitude: mel {st.mean(_ua):.0f}+/-{sd(_ua):.0f}%BL  "
             f"mel+mor {st.mean(_ub):.0f}+/-{sd(_ub):.0f}%BL  p={P_UPAMP:.3f}\n")

print("figures ->", OUT)
for p in PARS:
    v = sub(p)
    print(f"  {p:20s} rec={len(v):2d} mice={len(set(r['mouse'] for r in v))} "
          f"cells={sum(r['n'] for r in v):,} exc={st.mean([r['exc'] for r in v]):5.1f}% "
          f"inh={st.mean([r['inh'] for r in v]):5.1f}%")
print(f"  p_rec={P_REC:.3f}  p_mouse={P_MOUSE:.3f}  p_upamp={P_UPAMP:.3f}")
