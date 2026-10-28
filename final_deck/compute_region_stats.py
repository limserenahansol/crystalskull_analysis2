"""Region-level analysis of the melittin and melittin+morphine arms.

Five questions, in order:
  A  which Allen areas does melittin drive, up and down
  B  does melittin show the posterior-weighted inhibition face poke shows
  C  within one animal, is the morphine effect uniform across cortex or region-specific
  D  do the two paired animals agree on WHICH regions morphine changes
  E  does morphine act most where melittin acts most, or is that regression to the mean

Unit of analysis matters here and is easy to get wrong. mirror1 and mirror2 are two
separate recordings of the SAME imaging session, each with its own EXTRACT run and its
own registration - both cover the full bilateral FOV, so they are technical replicates,
not independent samples and not spatial halves. Tests that need an independent n use
sessions (6) or animals (4 melittin-only, 2 morphine); anything quoted over 12 labels
is labelled as such.

Writes region_stats.json for the figure script.

    python compute_region_stats.py
"""
from __future__ import annotations

import collections
import csv
import json
import os

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
CCF = os.path.join(os.path.dirname(HERE), "ccf_registration")

# isocortex. Used wherever a claim is about cortex; SCzo, root and the outside bin
# are excluded there, but kept in the cortex-wide base rate so the melittin numbers
# stay comparable with the face-poke slide, which was computed the same way.
ISO = {"MOs", "MOp", "SSp-m", "SSp-n", "SSp-bfd", "SSp-ul", "SSp-ll", "SSp-tr",
       "SSp-un", "SSs", "VISp", "VISa", "VISam", "VISal", "VISpm", "VISrl",
       "VISl", "VISpor", "RSPd", "RSPv", "RSPagl", "AUDp", "AUDd", "AUDv",
       "AUDpo", "ACAd", "ACAv", "TEa", "PL", "ILA", "ORBl", "ORBm", "ORBvl"}

PAIRS = [("CS1014-1b", "1014m01b"), ("CS1014-2w", "1014m02w")]
OUT = {}


def load(path, reckey):
    rows = list(csv.DictReader(open(path)))
    for r in rows:
        r["exc"] = r["is_exc"] == "1"
        r["inh"] = r["is_inh"] == "1"
        r["AP"] = float(r["AP_mm"])
        r["rec"] = r[reckey]
    return rows


MEL = load(os.path.join(CCF, "melittin_ccf", "cell_ccf_all.csv"), "label")
FP = load(os.path.join(CCF, "n18", "cell_ccf_all.csv"), "movie")
MEL_ONLY = [r for r in MEL if r["cohort"] == "melittin_only"]
MEL_MOR = [r for r in MEL if r["cohort"] == "melittin_morphine"]


def bh(pvals):
    """Benjamini-Hochberg q values, returned in the input order."""
    p = np.asarray(pvals, float)
    o = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    prev = 1.0
    for rank, i in enumerate(o[::-1]):
        prev = min(prev, p[i] * m / (m - rank))
        q[i] = prev
    return q


def enrich(rows, field, minn=300):
    """Fold enrichment per area vs that paradigm's own cortex-wide rate.

    Same recipe as the face-poke slide: base rate over every cell in the arm, areas
    filtered to those with at least `minn` cells, two-sided binomial per area, BH
    across the areas tested.
    """
    base = float(np.mean([r[field] for r in rows]))
    per = collections.defaultdict(list)
    for r in rows:
        if r["area_acronym"] not in ("root", "", "nan"):
            per[r["area_acronym"]].append(r[field])
    areas = [(a, v) for a, v in per.items() if len(v) >= minn]
    ps = [stats.binomtest(int(np.sum(v)), len(v), base).pvalue for _, v in areas]
    qs = bh(ps)
    out = [dict(area=a, k=float(np.mean(v)) / base, rate=100 * float(np.mean(v)),
                n=len(v), q=float(qq))
           for (a, v), qq in zip(areas, qs)]
    out.sort(key=lambda d: -d["k"])
    return out, base


# ==================================================================== A
print("=" * 78)
print("A   WHICH AREAS DOES MELITTIN DRIVE   (melittin-only arm)")
print("=" * 78)
print(f"{len(MEL_ONLY):,} cells, {len(set(r['mouse'] for r in MEL_ONLY))} mice, "
      f"{len(set(r['rec'] for r in MEL_ONLY))} recordings "
      f"({len(set(r['rec'] for r in MEL_ONLY))//2} sessions x 2 mirrors)")

up, base_up = enrich(MEL_ONLY, "exc")
dn, base_dn = enrich(MEL_ONLY, "inh")
print(f"\ncortex-wide: {100*base_up:.1f} % driven UP, {100*base_dn:.1f} % driven DOWN")
print(f"\nUP enrichment, {len(up)} areas with >= 300 cells")
for d in up:
    print(f"   {d['area']:9s} {d['k']:5.2f}x  rate {d['rate']:5.1f} %  "
          f"n={d['n']:6,}  q={d['q']:.1e}")
print(f"\nDOWN enrichment, {len(dn)} areas")
for d in dn[:6]:
    print(f"   {d['area']:9s} {d['k']:5.2f}x  rate {d['rate']:5.1f} %  "
          f"n={d['n']:6,}  q={d['q']:.1e}")
print("   ...")
for d in dn[-4:]:
    print(f"   {d['area']:9s} {d['k']:5.2f}x  rate {d['rate']:5.1f} %  "
          f"n={d['n']:6,}  q={d['q']:.1e}")

# are UP and DOWN maps mirror images, or independent?
du = {d["area"]: d["k"] for d in up}
dd = {d["area"]: d["k"] for d in dn}
sh = sorted(set(du) & set(dd))
rho_ud, p_ud = stats.spearmanr([du[a] for a in sh], [dd[a] for a in sh])
print(f"\nUP vs DOWN enrichment across the {len(sh)} shared areas: "
      f"Spearman rho = {rho_ud:+.2f}, p = {p_ud:.2g}")
print("   (strongly negative would mean the two maps are mirror images of each other,")
print("    i.e. one number, not two)")

print("\nper mouse, to check no single animal carries the pooled number:")
permouse = {}
for m in sorted({r["mouse"] for r in MEL_ONLY}):
    v = [r for r in MEL_ONLY if r["mouse"] == m]
    permouse[m] = dict(n=len(v), up=100 * float(np.mean([x["exc"] for x in v])),
                       dn=100 * float(np.mean([x["inh"] for x in v])))
    print(f"   {m:11s} n={len(v):6,}  up {permouse[m]['up']:5.1f} %  "
          f"down {permouse[m]['dn']:5.1f} %")

OUT["A"] = dict(up=up, down=dn, base_up=base_up, base_down=base_dn,
                rho_updown=float(rho_ud), p_updown=float(p_ud), per_mouse=permouse,
                n_cells=len(MEL_ONLY))

# ==================================================================== B
print("\n" + "=" * 78)
print("B   DOES MELITTIN SHOW THE POSTERIOR INHIBITION GRADIENT FACE POKE SHOWS")
print("=" * 78)


def gradient(rows, field, by="rec"):
    """Split each recording at its OWN median AP, so the test does not depend on
    where that FOV happened to sit. Returns (anterior %, posterior %) per unit."""
    ant, pos, keys = [], [], []
    for k in sorted({r[by] for r in rows}):
        v = [r for r in rows if r[by] == k]
        m = np.median([r["AP"] for r in v])
        a = [r[field] for r in v if r["AP"] >= m]
        p = [r[field] for r in v if r["AP"] < m]
        if not a or not p:
            continue
        ant.append(100 * float(np.mean(a)))
        pos.append(100 * float(np.mean(p)))
        keys.append(k)
    return np.array(ant), np.array(pos), keys


def report(tag, ant, pos, keys):
    if len(ant) < 2:
        print(f"   {tag}: n = {len(ant)}, no test possible")
        return dict(n=len(ant), ant=ant.tolist(), pos=pos.tolist(), keys=keys)
    W, p = stats.wilcoxon(ant, pos)
    same = int((pos > ant).sum())
    print(f"   {tag:38s} anterior {ant.mean():5.1f} %  posterior {pos.mean():5.1f} %  "
          f"delta {pos.mean()-ant.mean():+5.1f} pp   posterior higher in "
          f"{same}/{len(ant)}   W={W:.0f} p={p:.2g}")
    return dict(n=len(ant), ant=ant.tolist(), pos=pos.tolist(), keys=keys,
                W=float(W), p=float(p), same=same,
                delta=float(pos.mean() - ant.mean()))


# session = the mouse, since each mouse contributes one session x 2 mirrors per arm
for r in MEL:
    r["session"] = r["mouse"]

gradB = {}
print("\nface poke, for reference (the claim on slide 2):")
a, p_, k = gradient(FP, "inh")
gradB["facepoke_inh"] = report("face poke, inhibited, 18 recordings", a, p_, k)

print("\nmelittin only:")
a, p_, k = gradient(MEL_ONLY, "inh")
gradB["mel_inh_12"] = report("DOWN, 8 recordings (mirrors separate)", a, p_, k)
a, p_, k = gradient(MEL_ONLY, "inh", by="session")
gradB["mel_inh_4"] = report("DOWN, 4 sessions (mirrors pooled)", a, p_, k)
a, p_, k = gradient(MEL_ONLY, "exc")
gradB["mel_exc_12"] = report("UP, 8 recordings (mirrors separate)", a, p_, k)
a, p_, k = gradient(MEL_ONLY, "exc", by="session")
gradB["mel_exc_4"] = report("UP, 4 sessions (mirrors pooled)", a, p_, k)

print("\nmelittin + morphine:")
a, p_, k = gradient(MEL_MOR, "inh")
gradB["mor_inh_4"] = report("DOWN, 4 recordings (mirrors separate)", a, p_, k)
a, p_, k = gradient(MEL_MOR, "exc")
gradB["mor_exc_4"] = report("UP, 4 recordings (mirrors separate)", a, p_, k)

# continuous check, pooled, as a sanity test on the median split
for tag, rows in (("melittin only", MEL_ONLY), ("melittin+morphine", MEL_MOR),
                  ("face poke", FP)):
    ap = np.array([r["AP"] for r in rows])
    dn_ = np.array([r["inh"] for r in rows], float)
    rho, pp = stats.spearmanr(ap, dn_)
    gradB[f"cont_{tag}"] = dict(rho=float(rho), p=float(pp), n=len(rows))
    print(f"   continuous, pooled: {tag:20s} Spearman(AP, is_down) "
          f"rho = {rho:+.3f}, p = {pp:.1e}, n = {len(rows):,}")

OUT["B"] = gradB

# ==================================================================== C
print("\n" + "=" * 78)
print("C   WITHIN ONE ANIMAL, IS THE MORPHINE EFFECT UNIFORM OR REGION-SPECIFIC")
print("=" * 78)
print("Cochran's Q on the per-area log odds ratio. Significant Q = the effect differs")
print("by region. The sample here is CELLS within one animal, which is valid for")
print("asking whether regions differ inside that animal and says nothing about")
print("whether the pattern generalises - that is question D.")

bym = collections.defaultdict(list)
for r in MEL:
    bym[r["mouse"]].append(r)

het = {}
for mel_id, mor_id in PAIRS:
    rowsA, rowsB = bym[mel_id], bym[mor_id]
    per = {}
    for tag, rr in (("mel", rowsA), ("mor", rowsB)):
        d = collections.defaultdict(lambda: [0, 0])
        for r in rr:
            if r["area_acronym"] in ISO:
                d[r["area_acronym"]][0] += int(r["exc"])
                d[r["area_acronym"]][1] += 1
        per[tag] = d
    areas = sorted(a for a in per["mel"]
                   if per["mel"][a][1] >= 100 and per["mor"].get(a, [0, 0])[1] >= 100)
    rec = []
    for a in areas:
        c, nmel = per["mel"][a]
        k, nmor = per["mor"][a]
        # Haldane-Anscombe 0.5 correction guards the log when a cell count is 0
        aa, bb = k + .5, nmor - k + .5
        cc, ddv = c + .5, nmel - c + .5
        lor = np.log((aa * ddv) / (bb * cc))
        se = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / ddv)
        rec.append(dict(area=a, lor=float(lor), se=float(se),
                        n_mel=nmel, n_mor=nmor,
                        rate_mel=100 * c / nmel, rate_mor=100 * k / nmor,
                        d_pp=100 * (k / nmor - c / nmel),
                        p=float(2 * stats.norm.sf(abs(lor / se)))))
    qs = bh([r["p"] for r in rec])
    for r, qq in zip(rec, qs):
        r["q"] = float(qq)
    w = np.array([1 / r["se"] ** 2 for r in rec])
    lo = np.array([r["lor"] for r in rec])
    pool = float((w * lo).sum() / w.sum())
    se_pool = float(np.sqrt(1 / w.sum()))
    Q = float((w * (lo - pool) ** 2).sum())
    df = len(rec) - 1
    pQ = float(stats.chi2.sf(Q, df))
    I2 = float(max(0.0, (Q - df) / Q) * 100) if Q > 0 else 0.0
    rec.sort(key=lambda r: -r["lor"])
    print(f"\n   {mel_id} = {mor_id}   {len(rec)} isocortical areas "
          f"(>= 100 cells in both arms)")
    print(f"   pooled log OR = {pool:+.3f}  (OR {np.exp(pool):.3f}, 95% CI "
          f"{np.exp(pool-1.96*se_pool):.3f}-{np.exp(pool+1.96*se_pool):.3f})")
    print(f"   Cochran's Q = {Q:.1f}, df = {df}, p = {pQ:.2e}, I2 = {I2:.0f} %")
    print(f"     strongest increase  {rec[0]['area']:8s} logOR {rec[0]['lor']:+.2f}  "
          f"{rec[0]['rate_mel']:.1f} -> {rec[0]['rate_mor']:.1f} %  q={rec[0]['q']:.1e}")
    print(f"     strongest decrease  {rec[-1]['area']:8s} logOR {rec[-1]['lor']:+.2f}  "
          f"{rec[-1]['rate_mel']:.1f} -> {rec[-1]['rate_mor']:.1f} %  "
          f"q={rec[-1]['q']:.1e}")
    nsig = sum(1 for r in rec if r["q"] < .05)
    print(f"     {nsig}/{len(rec)} areas individually significant after BH")
    het[mel_id] = dict(areas=rec, pool=pool, se_pool=se_pool, Q=Q, df=df, pQ=pQ,
                       I2=I2, n_sig=nsig)

OUT["C"] = het

# ==================================================================== D
print("\n" + "=" * 78)
print("D   DO THE TWO ANIMALS AGREE ON WHICH REGIONS MORPHINE CHANGES")
print("=" * 78)

m1 = {r["area"]: r for r in het["CS1014-1b"]["areas"]}
m2 = {r["area"]: r for r in het["CS1014-2w"]["areas"]}
common = sorted(set(m1) & set(m2))
x_lor = np.array([m1[a]["lor"] for a in common])
y_lor = np.array([m2[a]["lor"] for a in common])
x_pp = np.array([m1[a]["d_pp"] for a in common])
y_pp = np.array([m2[a]["d_pp"] for a in common])

rP, pP = stats.pearsonr(x_lor, y_lor)
rS, pS = stats.spearmanr(x_lor, y_lor)
rP2, pP2 = stats.pearsonr(x_pp, y_pp)
rS2, pS2 = stats.spearmanr(x_pp, y_pp)
agree = int(np.sum(np.sign(x_lor) == np.sign(y_lor)))
pbin = float(stats.binomtest(agree, len(common), 0.5).pvalue)
ov = 0
for a in common:
    l1, s1 = m1[a]["lor"], m1[a]["se"]
    l2, s2 = m2[a]["lor"], m2[a]["se"]
    if not (l1 + 1.96 * s1 < l2 - 1.96 * s2 or l2 + 1.96 * s2 < l1 - 1.96 * s1):
        ov += 1

print(f"\n   {len(common)} areas present in both animals")
print(f"   log OR   Pearson r = {rP:+.2f} (p = {pP:.2f})   "
      f"Spearman rho = {rS:+.2f} (p = {pS:.2f})")
print(f"   delta pp Pearson r = {rP2:+.2f} (p = {pP2:.2f})   "
      f"Spearman rho = {rS2:+.2f} (p = {pS2:.2f})")
print(f"   same direction in {agree}/{len(common)} areas, "
      f"binomial vs 0.5 p = {pbin:.2f}")
print(f"   95% CIs overlap in {ov}/{len(common)} areas")

print(f"\n   {'area':9s} {'1b logOR':>10s} {'rank':>5s}   {'2w logOR':>10s} {'rank':>5s}"
      f"   same dir")
r1 = {a: i + 1 for i, a in enumerate(sorted(common, key=lambda a: -m1[a]["lor"]))}
r2 = {a: i + 1 for i, a in enumerate(sorted(common, key=lambda a: -m2[a]["lor"]))}
for a in sorted(common, key=lambda a: -(m1[a]["lor"] + m2[a]["lor"]) / 2):
    print(f"   {a:9s} {m1[a]['lor']:+10.2f} {r1[a]:5d}   {m2[a]['lor']:+10.2f} "
          f"{r2[a]:5d}   {'yes' if np.sign(m1[a]['lor'])==np.sign(m2[a]['lor']) else 'no'}")
rho_rank, p_rank = stats.spearmanr([r1[a] for a in common], [r2[a] for a in common])
print(f"\n   rank agreement: Spearman rho = {rho_rank:+.2f}, p = {p_rank:.2f}")

OUT["D"] = dict(areas=common,
                lor1=x_lor.tolist(), lor2=y_lor.tolist(),
                pp1=x_pp.tolist(), pp2=y_pp.tolist(),
                se1=[m1[a]["se"] for a in common], se2=[m2[a]["se"] for a in common],
                pearson_lor=[float(rP), float(pP)],
                spearman_lor=[float(rS), float(pS)],
                pearson_pp=[float(rP2), float(pP2)],
                spearman_pp=[float(rS2), float(pS2)],
                agree=agree, n=len(common), p_binom=pbin, ci_overlap=ov)

# ==================================================================== E
print("\n" + "=" * 78)
print("E   DOES MORPHINE ACT MOST WHERE MELITTIN ACTS MOST")
print("=" * 78)
print("x = that area's melittin-alone UP rate, y = morphine minus melittin (pp).")
print("A negative correlation is EXPECTED even under no real effect, because the noise")
print("in x enters y with a minus sign. The permutation null below quantifies that.")

rng = np.random.default_rng(0)
Eout = {}
for mel_id, _ in PAIRS:
    rec = het[mel_id]["areas"]
    x = np.array([r["rate_mel"] for r in rec])
    y = np.array([r["d_pp"] for r in rec])
    rP, pP = stats.pearsonr(x, y)
    rS, pS = stats.spearmanr(x, y)

    # null: the true up-rate is the same in both arms (the pooled rate), so any
    # correlation is the shared-noise artifact alone
    nmel = np.array([r["n_mel"] for r in rec])
    nmor = np.array([r["n_mor"] for r in rec])
    kmel = np.round(np.array([r["rate_mel"] for r in rec]) / 100 * nmel).astype(int)
    kmor = np.round(np.array([r["rate_mor"] for r in rec]) / 100 * nmor).astype(int)
    pool_p = (kmel + kmor) / (nmel + nmor)
    null = np.empty(2000)
    for i in range(2000):
        km = rng.binomial(nmel, pool_p) / nmel
        ko = rng.binomial(nmor, pool_p) / nmor
        null[i] = stats.pearsonr(100 * km, 100 * (ko - km))[0]
    pct = float((null <= rP).mean())
    print(f"\n   {mel_id}   {len(rec)} areas")
    print(f"     observed Pearson r = {rP:+.2f} (p = {pP:.3f}), "
          f"Spearman rho = {rS:+.2f} (p = {pS:.3f})")
    print(f"     null from shared noise alone: mean r = {null.mean():+.2f}, "
          f"95% range {np.percentile(null,2.5):+.2f} to {np.percentile(null,97.5):+.2f}")
    print(f"     observed sits at the {100*pct:.0f}th percentile of that null "
          f"-> {'NOT distinguishable from the artifact' if pct > .05 else 'beyond the artifact'}")
    Eout[mel_id] = dict(x=x.tolist(), y=y.tolist(),
                        areas=[r["area"] for r in rec],
                        r=float(rP), p=float(pP), rho=float(rS), rho_p=float(pS),
                        null_mean=float(null.mean()),
                        null_lo=float(np.percentile(null, 2.5)),
                        null_hi=float(np.percentile(null, 97.5)),
                        pct=pct)

OUT["E"] = Eout

json.dump(OUT, open(os.path.join(HERE, "region_stats.json"), "w"), indent=1)
print(f"\nwrote {os.path.join(HERE, 'region_stats.json')}")
