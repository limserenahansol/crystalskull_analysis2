"""Melittin vs melittin+morphine as a WITHIN-ANIMAL comparison.

Hansol confirmed that 1014m01b and 1014m02w are the same animals that appear in the
melittin-only cohort as CS1014-1b and CS1014-2w. Every previous version of this
comparison treated the two arms as independent groups of mice (n = 4 vs 2, unpaired
Mann-Whitney). For the two animals that appear in both arms that was the wrong design -
a paired comparison removes between-animal variance, which is the dominant noise source
here.

What this does NOT fix: there are still only two paired animals. A two-sided sign test on
n = 2 has a floor of p = 0.5, and a paired Wilcoxon on n = 2 cannot go below p = 0.5
either. So the pairing improves the estimate, not the evidence. The per-animal deltas are
reported because they are interpretable; no p value is quoted for the paired test because
none would be meaningful.

Output -> figures_v2/fig17_paired_morphine.png  and  paired_morphine.csv
"""
from __future__ import annotations

import collections
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures_v2")

# confirmed by Hansol 2026-09-29
PAIRS = {"CS1014-1b": "1014m01b", "CS1014-2w": "1014m02w"}
RED, BLUE, TEAL, AMBER, GREY = "#C0392B", "#2C6FBB", "#0D5C63", "#D08C1E", "#555555"
CTX = {"MOs", "MOp", "SSp-m", "SSp-n", "SSp-bfd", "SSp-ul", "SSp-ll", "SSp-tr", "SSp-un",
       "SSs", "VISp", "VISa", "VISam", "VISal", "VISpm", "VISrl", "VISl", "VISpor",
       "RSPd", "RSPv", "RSPagl", "AUDp", "AUDd", "AUDv", "ACAd", "TEa"}

rows = list(csv.DictReader(open(os.path.join(HERE, "melittin_ccf", "cell_ccf_all.csv"))))
for r in rows:
    r["exc"] = r["is_exc"] == "1"
    r["inh"] = r["is_inh"] == "1"
    r["amp_pct"] = float(r["amp_pct"])

by_mouse = collections.defaultdict(list)
for r in rows:
    by_mouse[r["mouse"]].append(r)

print("cells per mouse:")
for m, v in sorted(by_mouse.items()):
    arm = "melittin+morphine" if v[0]["cohort"] == "melittin_morphine" else "melittin only"
    print(f"  {m:12s} {len(v):7,}  {arm}")

# ---------------------------------------------------------------- whole-cortex, paired
print("\n" + "=" * 74)
print("WITHIN-ANIMAL: melittin alone -> melittin + morphine")
print("=" * 74)
recs = []
for mel_id, mor_id in PAIRS.items():
    a, b = by_mouse.get(mel_id), by_mouse.get(mor_id)
    if not a or not b:
        print(f"  {mel_id} / {mor_id}: missing"); continue
    d = {
        "animal": f"{mel_id} = {mor_id}",
        "n_cells_mel": len(a), "n_cells_mor": len(b),
        "pct_up_mel": 100 * np.mean([x["exc"] for x in a]),
        "pct_up_mor": 100 * np.mean([x["exc"] for x in b]),
        "pct_down_mel": 100 * np.mean([x["inh"] for x in a]),
        "pct_down_mor": 100 * np.mean([x["inh"] for x in b]),
        "amp_mel": float(np.mean([x["amp_pct"] for x in a])),
        "amp_mor": float(np.mean([x["amp_pct"] for x in b])),
    }
    d["d_up"] = d["pct_up_mor"] - d["pct_up_mel"]
    d["d_down"] = d["pct_down_mor"] - d["pct_down_mel"]
    d["d_amp"] = d["amp_mor"] - d["amp_mel"]
    recs.append(d)

print(f"\n{'animal':<26}{'% up':>18}{'% down':>18}{'mean amp (% base)':>22}")
print(f"{'':<26}{'mel -> +mor':>18}{'mel -> +mor':>18}{'mel -> +mor':>22}")
print("-" * 84)
for d in recs:
    print(f"{d['animal']:<26}"
          f"{d['pct_up_mel']:6.1f} ->{d['pct_up_mor']:6.1f} "
          f"({d['d_up']:+.1f})"
          f"{d['pct_down_mel']:6.1f} ->{d['pct_down_mor']:6.1f} "
          f"({d['d_down']:+.1f})"
          f"{d['amp_mel']:+8.1f} ->{d['amp_mor']:+7.1f} ({d['d_amp']:+.1f})")

print(f"\nboth animals move the same way on % up:   "
      f"{'yes' if all(d['d_up'] > 0 for d in recs) else 'no'}")
print(f"both animals move the same way on amp:    "
      f"{'yes' if len({np.sign(d['d_amp']) for d in recs}) == 1 else 'no'}")
print("\nn = 2 paired animals. A paired Wilcoxon or sign test on n = 2 has a two-sided")
print("floor of p = 0.5, so no p value is quoted - the pairing sharpens the estimate,")
print("it does not supply evidence. Two more paired animals would.")

# ---------------------------------------------------------------- by region, paired
print("\n" + "=" * 74)
print("WITHIN-ANIMAL BY REGION  (delta = morphine minus melittin, percentage points)")
print("=" * 74)
reg = collections.defaultdict(dict)
for mel_id, mor_id in PAIRS.items():
    for tag, mid in (("mel", mel_id), ("mor", mor_id)):
        per = collections.defaultdict(list)
        for r in by_mouse.get(mid, []):
            if r["area_acronym"] in CTX:
                per[r["area_acronym"]].append(r)
        for a, v in per.items():
            if len(v) >= 100:
                reg[a][f"{tag}_{mel_id}"] = 100 * np.mean([x["exc"] for x in v])

area_rows = []
for a, d in reg.items():
    ds = [d[f"mor_{m}"] - d[f"mel_{m}"] for m in PAIRS if f"mor_{m}" in d and f"mel_{m}" in d]
    if len(ds) == 2:
        area_rows.append((a, ds[0], ds[1], float(np.mean(ds)),
                          "yes" if ds[0] * ds[1] > 0 else "no"))
area_rows.sort(key=lambda t: -t[3])
print(f"{'area':<10}{'d 1b':>9}{'d 2w':>9}{'mean':>9}  same direction in both?")
for a, d1, d2, m, agree in area_rows:
    print(f"{a:<10}{d1:+9.1f}{d2:+9.1f}{m:+9.1f}  {agree}")
n_ag = sum(1 for r in area_rows if r[4] == "yes")
print(f"\n{n_ag} of {len(area_rows)} areas move the same way in both animals "
      f"({100*n_ag/max(len(area_rows),1):.0f}%)")

with open(os.path.join(OUT, "paired_morphine.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["area", "delta_up_CS1014-1b", "delta_up_CS1014-2w", "mean_delta",
                "same_direction"])
    w.writerows(area_rows)

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(15, 7.2))
gs = fig.add_gridspec(1, 3, width_ratios=[.8, .8, 1.8], wspace=.36)

for k, (key_mel, key_mor, ttl, ylab) in enumerate(
        [("pct_up_mel", "pct_up_mor", "% of cells up vs baseline", "% up"),
         ("amp_mel", "amp_mor", "mean amplitude", "% of own baseline")]):
    ax = fig.add_subplot(gs[k])
    for d in recs:
        ax.plot([0, 1], [d[key_mel], d[key_mor]], "-o", color=GREY, ms=9, lw=2,
                markerfacecolor="w", markeredgewidth=2)
        ax.annotate(d["animal"].split(" = ")[0], (1, d[key_mor]),
                    xytext=(8, 0), textcoords="offset points", fontsize=9, va="center")
    ax.set_xlim(-.35, 1.75); ax.set_xticks([0, 1])
    ax.set_xticklabels(["melittin\nalone", "melittin +\nmorphine"])
    ax.set_ylabel(ylab)
    ax.set_title(ttl, fontsize=11, fontweight="bold")
    ax.axhline(0, color="k", ls=":", lw=1)
    ax.spines[["top", "right"]].set_visible(False)

ax = fig.add_subplot(gs[2])
y = np.arange(len(area_rows))
ax.barh(y - .19, [r[1] for r in area_rows], .36, color=TEAL, label="CS1014-1b = 1014m01b")
ax.barh(y + .19, [r[2] for r in area_rows], .36, color=AMBER, label="CS1014-2w = 1014m02w")
ax.set_yticks(y)
ax.set_yticklabels([r[0] for r in area_rows], fontsize=8.5)
ax.invert_yaxis()
ax.axvline(0, color="k", lw=1)
ax.set_xlabel("change in % up,  morphine minus melittin  (percentage points)")
ax.set_title(f"By Allen area, each animal separately\n"
             f"{n_ag}/{len(area_rows)} areas agree in direction",
             fontsize=11, fontweight="bold")
ax.legend(frameon=False, fontsize=9, loc="lower right")
ax.spines[["top", "right"]].set_visible(False)

fig.suptitle("Fig 17  |  Melittin vs melittin + morphine, WITHIN animal\n"
             "the same two mice under both arms  -  n = 2 pairs, so no p value is quoted",
             fontsize=13.5, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .84])
fig.savefig(os.path.join(OUT, "fig17_paired_morphine.png"), dpi=160)
plt.close(fig)
print(f"\nwrote {OUT}/fig17_paired_morphine.png")
