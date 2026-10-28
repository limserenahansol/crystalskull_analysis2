"""Exact cell counts per Allen region, per condition, per mouse.

Counts, not percentages. Three levels of grouping, because "how many inhibited neurons in
SSp" can mean the leaf area (SSp-bfd) or the parent (all SSp-* summed):

  leaf    SSp-bfd, SSp-m, VISp, MOs ...          as labelled per cell
  parent  SSp, VIS, MO, RSP, AUD, SSs, TEa ...   leaves rolled up
  per-mouse   the same, broken down by animal

Outputs -> counts/
  counts_by_area_leaf.csv
  counts_by_area_parent.csv
  counts_by_area_and_mouse.csv
"""
import collections
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "counts")
os.makedirs(OUT, exist_ok=True)


def parent_of(acr):
    """Roll a leaf acronym up to its functional parent."""
    if acr.startswith("SSp-"):
        return "SSp"
    if acr.startswith("VIS") and acr != "VISC":
        return "VIS"
    if acr in ("MOp", "MOs"):
        return "MO"
    if acr.startswith("RSP"):
        return "RSP"
    if acr.startswith("AUD"):
        return "AUD"
    return acr


def load(path, cond_fn, mouse_fn):
    out = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            out.append(dict(area=r["area_acronym"], cond=cond_fn(r), mouse=mouse_fn(r),
                            exc=r["is_exc"] == "1", inh=r["is_inh"] == "1"))
    return out


def poke_mouse(r):
    mv = r["movie"]
    return "m1" if mv.startswith("0004749770_m1_") else mv.split("_")[0]


rows = load(os.path.join(HERE, "n18", "cell_ccf_all.csv"),
            lambda r: "Face poke", poke_mouse)
rows += load(os.path.join(HERE, "melittin_ccf", "cell_ccf_all.csv"),
             lambda r: "Melittin+morphine" if r["cohort"] == "melittin_morphine" else "Melittin",
             lambda r: r["mouse"])

CONDS = ["Face poke", "Melittin", "Melittin+morphine"]
SKIP = {"(outside)"}


def tally(key_fn):
    t = collections.defaultdict(lambda: dict(n=0, e=0, i=0))
    for r in rows:
        if r["area"] in SKIP:
            continue
        k = key_fn(r)
        t[k]["n"] += 1
        t[k]["e"] += r["exc"]
        t[k]["i"] += r["inh"]
    return t


def write(path, tal, keycols):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(list(keycols) + ["n_cells", "n_excited", "n_inhibited", "n_nonsig",
                                    "pct_excited", "pct_inhibited"])
        for k in sorted(tal, key=lambda x: (x[0], -tal[x]["n"])):
            d = tal[k]
            w.writerow(list(k) + [d["n"], d["e"], d["i"], d["n"] - d["e"] - d["i"],
                                  f"{100*d['e']/d['n']:.2f}", f"{100*d['i']/d['n']:.2f}"])


leaf = tally(lambda r: (r["cond"], r["area"]))
par = tally(lambda r: (r["cond"], parent_of(r["area"])))
pm = tally(lambda r: (r["cond"], parent_of(r["area"]), r["mouse"]))

write(os.path.join(OUT, "counts_by_area_leaf.csv"), leaf, ["condition", "area"])
write(os.path.join(OUT, "counts_by_area_parent.csv"), par, ["condition", "area_group"])
write(os.path.join(OUT, "counts_by_area_and_mouse.csv"), pm,
      ["condition", "area_group", "mouse"])


def show(tal, title, keylen=2, minn=100):
    print(f"\n{'='*78}\n{title}\n{'='*78}")
    for cond in CONDS:
        sub = {k: v for k, v in tal.items() if k[0] == cond and v["n"] >= minn}
        if not sub:
            continue
        tot = sum(v["n"] for v in sub.values())
        te = sum(v["e"] for v in sub.values())
        ti = sum(v["i"] for v in sub.values())
        print(f"\n{cond}   ({tot:,} cells, {te:,} excited, {ti:,} inhibited)")
        print(f"  {'area':<10} {'cells':>8} {'excited':>10} {'inhibited':>12} {'no change':>11}")
        print("  " + "-" * 54)
        for k in sorted(sub, key=lambda x: -sub[x]["n"]):
            d = sub[k]
            print(f"  {k[1]:<10} {d['n']:>8,} {d['e']:>6,} ({100*d['e']/d['n']:4.1f}%) "
                  f"{d['i']:>6,} ({100*d['i']/d['n']:4.1f}%) {d['n']-d['e']-d['i']:>10,}")
        print(f"  {'TOTAL':<10} {tot:>8,} {te:>6,} ({100*te/tot:4.1f}%) "
              f"{ti:>6,} ({100*ti/tot:4.1f}%) {tot-te-ti:>10,}")


show(par, "COUNTS BY AREA GROUP  (SSp = all SSp-* summed, VIS = all visual, MO = MOp+MOs)")
show(leaf, "COUNTS BY LEAF AREA", minn=200)

print(f"\n\nCSVs written to {OUT}")
for f in sorted(os.listdir(OUT)):
    print("  ", f)
