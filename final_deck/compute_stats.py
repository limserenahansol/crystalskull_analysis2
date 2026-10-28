"""Recompute every number that goes on the 5 PI slides, from the cell tables.

Nothing here is copied from an earlier script. If a number on a slide does not appear
in this printout, it does not go on the slide.
"""
from __future__ import annotations
import collections, csv, json, math, os, re
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
CCF  = os.path.join(os.path.dirname(HERE), "ccf_registration")

def load(path, movie_key):
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        r["exc"] = r["is_exc"] == "1"
        r["inh"] = r["is_inh"] == "1"
        r["AP"]  = float(r["AP_mm"]); r["ML"] = float(r["ML_mm"])
        r["rec"] = r[movie_key]
    return rows

fp  = load(os.path.join(CCF, "n18", "cell_ccf_all.csv"), "movie")
mel = load(os.path.join(CCF, "melittin_ccf", "cell_ccf_all.csv"), "mouse")

# animal identity for face poke comes from the movie prefix, NOT the 'mouse'
# column (that column holds the FOV label m1f / m2c and collapses 4 animals to 3).
def fp_animal(movie):
    return re.split(r"_d\d{6}_", movie)[0]

out = {}

# ---------------------------------------------------------------- registration
tm = list(csv.DictReader(open(os.path.join(CCF, "transforms", "transform_matrices.csv"))))
scale = np.array([math.sqrt(abs(float(r["a11"])*float(r["a22"]) - float(r["a12"])*float(r["a21"])))*25
                  for r in tm])
dice = np.array([float(r["dice"]) for r in tm]) if "dice" in tm[0] else None
out["reg"] = dict(n=len(tm), um_med=float(np.median(scale)), um_min=float(scale.min()),
                  um_max=float(scale.max()), spread=float(scale.max()/scale.min()))
if dice is not None:
    out["reg"].update(dice_med=float(np.median(dice)), dice_min=float(dice.min()))
print("REGISTRATION", json.dumps(out["reg"], indent=2))

# ---------------------------------------------------------------- face poke
recs = sorted({r["rec"] for r in fp})
anis = sorted({fp_animal(r["rec"]) for r in fp})
out["fp"] = dict(n_cells=len(fp), n_recs=len(recs), n_animals=len(anis),
                 animals=anis,
                 pct_exc=100*np.mean([r["exc"] for r in fp]),
                 pct_inh=100*np.mean([r["inh"] for r in fp]))
print(f"\nFACE POKE  {len(fp):,} cells  {len(recs)} recordings  {len(anis)} animals {anis}")
print(f"  excited  {out['fp']['pct_exc']:.1f} %   inhibited {out['fp']['pct_inh']:.1f} %")

# A-P gradient: per recording, inhibited fraction anterior vs posterior of that
# recording's own median AP. Paired across recordings -> Wilcoxon signed-rank.
ant, pos = [], []
for rc in recs:
    v = [r for r in fp if r["rec"] == rc]
    m = np.median([r["AP"] for r in v])
    a = [r["inh"] for r in v if r["AP"] >= m]     # AP_mm positive = anterior
    p = [r["inh"] for r in v if r["AP"] <  m]
    ant.append(100*np.mean(a)); pos.append(100*np.mean(p))
ant, pos = np.array(ant), np.array(pos)
W, pW = stats.wilcoxon(ant, pos)
out["fp"].update(ap_ant=ant.tolist(), ap_pos=pos.tolist(),
                 ap_same=int((ant > pos).sum()), ap_n=len(ant),
                 ap_W=float(W), ap_p=float(pW),
                 ap_delta=float(np.mean(ant-pos)))
print(f"  A-P gradient (inhibited %): anterior {ant.mean():.1f} vs posterior {pos.mean():.1f}"
      f"  delta {np.mean(ant-pos):+.1f} pp")
print(f"    same direction in {int((ant>pos).sum())}/{len(ant)} recordings, "
      f"Wilcoxon signed-rank W={W:.0f} p={pW:.2g}")

# ---------------------------------------------------------------- melittin
mo = [r for r in mel if r["cohort"] == "melittin_morphine"]
mn = [r for r in mel if r["cohort"] == "melittin_only"]
for tag, v in (("melittin only", mn), ("melittin + morphine", mo)):
    print(f"\n{tag:22s} {len(v):,} cells  {len(set(x['mouse'] for x in v))} mice "
          f"{sorted(set(x['mouse'] for x in v))}")
    print(f"  up {100*np.mean([x['exc'] for x in v]):.1f} %   "
          f"down {100*np.mean([x['inh'] for x in v]):.1f} %")
out["mel"] = dict(
    only=dict(n=len(mn), mice=sorted({x["mouse"] for x in mn}),
              up=100*np.mean([x["exc"] for x in mn]), down=100*np.mean([x["inh"] for x in mn])),
    mor =dict(n=len(mo), mice=sorted({x["mouse"] for x in mo}),
              up=100*np.mean([x["exc"] for x in mo]), down=100*np.mean([x["inh"] for x in mo])))

# paired, within animal
PAIRS = {"CS1014-1b": "1014m01b", "CS1014-2w": "1014m02w"}
bym = collections.defaultdict(list)
for r in mel: bym[r["mouse"]].append(r)
pair_rows = []
for a, b in PAIRS.items():
    A, B = bym[a], bym[b]
    pair_rows.append(dict(animal=a, n_mel=len(A), n_mor=len(B),
        up_mel=100*np.mean([x["exc"] for x in A]), up_mor=100*np.mean([x["exc"] for x in B]),
        dn_mel=100*np.mean([x["inh"] for x in A]), dn_mor=100*np.mean([x["inh"] for x in B]),
        amp_mel=float(np.mean([float(x["amp_pct"]) for x in A])),
        amp_mor=float(np.mean([float(x["amp_pct"]) for x in B]))))
for d in pair_rows:
    d["d_up"] = d["up_mor"]-d["up_mel"]; d["d_amp"] = d["amp_mor"]-d["amp_mel"]
out["paired"] = pair_rows
print("\nPAIRED (same animal, both arms)")
for d in pair_rows:
    print(f"  {d['animal']:10s} up {d['up_mel']:5.1f} -> {d['up_mor']:5.1f} ({d['d_up']:+.1f} pp)   "
          f"amp {d['amp_mel']:+5.1f} -> {d['amp_mor']:+5.1f} ({d['d_amp']:+.1f})")
print(f"  mean delta up = {np.mean([d['d_up'] for d in pair_rows]):+.1f} pp, n = 2 pairs, no p value")

# ---------------------------------------------------------------- enrichment
CTX_MIN = 300
def enrich(rows, field):
    base = np.mean([r[field] for r in rows])
    per = collections.defaultdict(list)
    for r in rows:
        if r["area_acronym"] not in ("root", "", "nan"): per[r["area_acronym"]].append(r[field])
    o = []
    for a, v in per.items():
        if len(v) >= CTX_MIN:
            k = np.mean(v)/base
            # binomial test of that area against the rest of the cortex
            x = int(np.sum(v)); n = len(v)
            p = stats.binomtest(x, n, base).pvalue
            o.append((a, k, n, p))
    # Benjamini-Hochberg
    o.sort(key=lambda t: t[3])
    m = len(o)
    q, prev = [], 1.0
    for i, t in enumerate(reversed(o)):
        prev = min(prev, t[3]*m/(m-i)); q.append(prev)
    q = q[::-1]
    return [(a, k, n, p, qq) for (a, k, n, p), qq in zip(o, q)]

for tag, rows, field in (("FACE POKE inhibited", fp, "inh"),
                         ("MELITTIN up",         mn, "exc")):
    e = enrich(rows, field)
    e.sort(key=lambda t: -t[1])
    print(f"\n{tag}  enrichment vs that paradigm's own cortex-wide rate "
          f"(areas with >= {CTX_MIN} cells, n = {len(e)} areas, binomial + BH-FDR)")
    for a, k, n, p, qq in e[:6]:
        print(f"   {a:9s} {k:5.2f}x  n={n:6d}  q={qq:.2g}")
    print("   ...")
    for a, k, n, p, qq in e[-4:]:
        print(f"   {a:9s} {k:5.2f}x  n={n:6d}  q={qq:.2g}")
    out.setdefault("enrich", {})[tag] = [dict(area=a, k=k, n=n, q=qq) for a, k, n, p, qq in e]

json.dump(out, open(os.path.join(HERE, "stats.json"), "w"), indent=1)
print(f"\nwrote {os.path.join(HERE, 'stats.json')}")
