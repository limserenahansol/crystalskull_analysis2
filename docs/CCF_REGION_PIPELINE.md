# Allen CCF registration and region-level analysis

How a crystal-skull FOV becomes a table of cells with Allen area labels, and what you
can and cannot conclude from that table.

Two paradigms run through this pipeline:

| arm | cells | animals | recordings | response test |
|---|---|---|---|---|
| face poke | 43,156 | 4 | 18 | within-cell, trial-paired |
| melittin only | 81,834 | 4 | 8 (4 sessions × 2 mirrors) | across-session, cells not tracked |
| melittin + morphine | 62,533 | 2 | 4 (2 sessions × 2 mirrors) | across-session, cells not tracked |

The two response tests are **not** the same and their percentages are not comparable.
See [Response classification](#response-classification) before quoting any of them.

---

## 1. Registration — no manual landmarks

`ccf_registration/make_allen_dorsal_map.py` → `allen_dorsal_map.mat`

Projects the BrainGlobe `allen_mouse_25um_v1.2` annotation volume (orientation `asr`,
528 × 320 × 456) top-down. The first labelled voxel along each dorsal ray gives the
surface area; cortical layers collapse to their parent area, so `SSp-bfd2/3` and
`SSp-bfd5` both become `SSp-bfd`. Bregma lands at atlas row 216, column 228.

`ccf_registration/ccf_register_by_mask.m` → `transforms/*.mat`

A similarity transform — rotation, one isotropic scale, translation — fitted from the
FOV **silhouette only**. No cell coordinates enter the fit.

1. rotation = auto-detected sagittal midline (`ccf_autodetect_midline.m`) plus the
   fixed 90° between the FOV frame (A–P horizontal) and the atlas frame (A–P vertical)
2. scale + translation = maximise Dice overlap between the binarised tissue mask and
   the atlas cortex mask
3. A–P sense = run both, keep the better Dice

Each `transforms/<label>.mat` carries `dice`, `um_per_px`, `midline_angle_deg`,
`ap_flipped`, `flipLR`, `fov_size`, and a `_qc.png` beside it.

**Why not bregma and lambda.** The skull is removed in this preparation, so neither is
visible. An initial clicked-landmark pass gave a 1.69× spread in fitted scale and an
implied FOV width of 24.6 mm, which is impossible. Mask fitting gives 1.05×.

**Result over all 30 recordings.** Dice median 0.931 (range 0.912–0.942); scale median
5.48 µm/px (5.38–5.65, spread 1.05×); midline angles −4.0° to +7.5°. Cerebellum
mislabels among posterior cells dropped from 25 % to 0 %, and three known stereotaxic
coordinates each land in the right area (3.8 P / 2.5 L → `VISp`; 1.5 P / 3.2 L →
`SSp-bfd`; 1.0 A / 1.5 L → `MOp`).

**A fixed bug worth knowing about.** The first midline detector succeeded on 12/30 and
reported confidence 1.00 everywhere, because the score was clipped at 1 and the
connected-component step split one hemisphere. The projection-based detector with an
unclipped score succeeds on 30/30.

## 2. Labelling cells

```
ccf_registration/ccf_label_facepoke_n18.m   ->  n18/cell_ccf_all.csv
ccf_registration/ccf_label_melittin.m       ->  melittin_ccf/cell_ccf_all.csv
```

Applies the transform to each EXTRACT footprint centroid and writes one row per cell:
`cell, atlas_row, atlas_col, AP_mm, ML_mm, hemisphere, area_acronym, area_name, label,
mouse, [cohort|movie], amp, pval, [amp_pct], is_exc, is_inh`.

`AP_mm` is positive anterior. `label` is the recording; `mouse` is the animal in the
melittin tables but the **FOV label** (`m1f` / `m2c`) in the face-poke table, where the
animal has to come from the movie-name prefix instead — using the column collapses the
4 animals to 3.

`ccf_label_melittin.m` carries its own copy of `local_frr_spatial` because the melittin
`frr` layout differs from what `extract_load_ext_spatial` expects.

### Mirror views are not independent recordings

`mirror1` and `mirror2` are two separate recordings of the **same imaging session**,
each with its own EXTRACT run and its own registration, both covering the full
bilateral FOV (ML −5 to +5 mm). They are technical replicates. Any test that needs an
independent n uses sessions or animals, never the 12 labels.

### Response classification

**Face poke** — within-cell and trial-paired (`extract_stim_peristim_core.m`). Per cell
and per poke, amplitude = mean dF/F over 0 to +2 s minus mean dF/F over −2 to 0 s;
the per-cell value is the mean of that difference across trials, tested with a
two-sided t **across trials**, uncorrected p < 0.05. Peri-stimulus window −2 to +3 s;
only pokes whose full window lies inside the recording are used.
→ 8.6 % excited, 27.6 % inhibited.

**Melittin** — baseline and injection are separate EXTRACT runs, so cells are **not
tracked across the injection**. Amplitude = that injection cell's mean dF/F minus the
grand mean over all baseline cells, tested with a t **across time bins**, uncorrected
p < 0.05.
→ 28.3 % up, 63.0 % down.

The melittin test uses time bins as its sample, so its degrees of freedom run into the
thousands and nearly every cell crosses p < 0.05 in one direction or the other —
91.3 % of cells are called responsive. That is why the up-enrichment and
down-enrichment maps correlate at Spearman ρ = −0.98: they are one map with the sign
flipped.

So **28 % / 63 % and 8.6 % / 27.6 % are not the same measurement.** What stays valid is
any comparison where the identical test is applied to both sides — between areas within
one arm, between the two drug arms, or between paradigms after each is divided by its
own cortex-wide rate.

The fix, if someone wants these arms directly comparable: register the baseline and
injection EXTRACT runs to each other so cells can be tracked, then run the melittin arm
as a within-cell paired test. That is the single biggest methodological upgrade
available on this dataset.

## 3. Region-level statistics

`final_deck/compute_stats.py` and `final_deck/compute_region_stats.py` recompute every
number that appears on a slide, straight from the cell tables, into `stats.json` and
`region_stats.json`. Nothing downstream hard-codes a number.

**Fold enrichment.** An area's response rate divided by that arm's own cortex-wide
rate. Areas with ≥ 300 cells only. Two-sided binomial test of each area against the
cortex-wide rate, Benjamini–Hochberg across the areas tested (19 for face poke, 27 for
melittin). Dividing by each arm's own base rate is what makes the paradigms comparable
despite base rates of 8.6 % and 28.3 %.

**A–P gradient.** Each recording is split at its **own** median `AP_mm`, so the result
does not depend on where that FOV happened to sit, then Wilcoxon signed-rank on the
paired anterior/posterior fractions.

Face poke: inhibition is higher posteriorly in 18/18 recordings, +7.2 pp, W = 0,
p = 7.6e-06. Melittin: +2.2 pp, 3/4 mice, p = 0.38 — and note a two-sided signed-rank
test on n = 4 **floors at p = 0.125**, so it could not have reached significance even
if all four animals had agreed. Not demonstrable ≠ absent.

**Within-animal region heterogeneity** (morphine). Per isocortical area with ≥ 100
cells in both arms, the 2 × 2 table (up / not-up) × (melittin / melittin+morphine)
gives a log odds ratio and its SE, with a Haldane–Anscombe 0.5 correction so a zero
count cannot break the log. Inverse-variance pooling gives the animal's overall log OR.
Cochran's Q = Σ wᵢ(logORᵢ − logOR_pool)², wᵢ = 1/SEᵢ², df = K−1, χ²; I² = (Q−df)/Q.

The sample in that test is **cells inside one animal**. It is the right unit for asking
whether regions differ within that animal and the wrong unit for asking whether
morphine is region-specific in general. Treating cells as independent animals is
pseudoreplication and would turn Q = 129 into a false headline.

**Cross-animal reproducibility** (morphine). The per-area log ORs of the two paired
animals, correlated across the 22 shared areas, with direction agreement tested against
chance by a two-sided binomial and per-area CI overlap counted.

## 4. What the morphine arm actually shows

Within the two animals recorded under both arms (`CS1014-1b` = `1014m01b`,
`CS1014-2w` = `1014m02w`):

| | CS1014-1b | CS1014-2w |
|---|---|---|
| driven up % | 28.5 → 31.3 (+2.8) | 28.9 → 31.6 (+2.6) |
| driven down % | 63.9 → 62.5 (−1.4) | 60.2 → 58.7 (−1.5) |
| no response % | 7.6 → 6.2 | 10.9 → 9.8 |
| mean amplitude (% of baseline) | −5.7 → −3.3 (+2.4) | +0.5 → +2.7 (+2.3) |
| pooled odds ratio | 1.15 [1.10–1.19] | 1.16 [1.12–1.20] |

All four measures move the same way: morphine makes the melittin response slightly
**larger**, not smaller.

Within each animal the effect is strongly region-dependent (Q = 129 and 278, I² = 82 %
and 92 %). Across the two animals it is not: r = +0.10 (p = 0.65), 11/22 directions
agree, CIs overlap in 14/22. **No individual area can be named from n = 2.**

### Two limits that are not statistical

**There is no morphine-alone arm and no morphine + face-poke arm.** Morphine exists
only as melittin + morphine.

**The two arms are seven weeks apart** — melittin alone on 20260218/19, melittin +
morphine on 20260408 — and there is no melittin + vehicle session on the second day. A
+2.7 pp shift is equally consistent with morphine, with window maturation over seven
weeks, or with a second melittin exposure. This design cannot separate them. A same-day
vehicle arm, ideally a within-session crossover, would.

## 5. Figures and decks

```
final_deck/make_final_figures.py    S1..S4      registration, face poke, morphine, dissociation
final_deck/make_s5.py               S5          behaviour video + grimace feasibility
final_deck/make_region_figures.py   S6, S7      melittin by region, morphine region-specificity
final_deck/make_slide4b.py          S4b         what morphine changed, + a one-slide pptx
final_deck/make_final_deck.py                   the 5-slide deck + supplements S1-S5
final_deck/make_region_deck.py                  the 2 region slides + supplement S6, and a merged deck
final_deck/preview_deck.py                      layout proof: redraws each slide and flags overflow
```

Decks are built with **python-pptx only**. Never use PowerPoint COM — it closes
whatever the user has open.

There is no LibreOffice on the analysis machine, so a `.pptx` cannot be rasterised
directly. `preview_deck.py` walks the shapes, redraws each slide in matplotlib at the
real geometry and font sizes, and outlines in red any text box that runs past its
slide. It is not pixel-identical to PowerPoint; it exists to answer one question — does
anything overflow.

Run order for a full rebuild:

```bash
cd final_deck
python compute_stats.py          # -> stats.json, reg_qc.csv
python compute_region_stats.py   # -> region_stats.json
python make_final_figures.py     # -> figures/S1..S4
python make_s5.py                # -> figures/S5
python make_region_figures.py    # -> figures/S6, S7
python make_slide4b.py           # -> figures/S4b + one-slide pptx
python make_final_deck.py        # -> the 5-slide deck
python make_region_deck.py       # -> the add-on and merged decks
python preview_deck.py           # -> preview/, and an overflow report
```

## 6. Behaviour video

`behavior_sync/ccf_behavior_index.py` → `behavior_calcium_manifest.csv`

Indexes the 19 face-video sessions on the share against the calcium recordings:
5 `.mp4` from 20260408, 14 loose `.jpg` sequences from 20260218/19, 15 of 19 with a
matching calcium folder.

**Frame index × 1/35 s is not wall clock.** `1014m01b 20260408 Injection` holds 1,504 s
of frames spread over 3,761 s of wall clock, so that recording is segmented. Always use
the per-frame timestamps.

Automated grimace scoring does not transfer to this footage. With the published GrimACE
weights (Zenodo 15119195) in the authors' own code, the face cropper accepts 6/12
sampled frames of GrimACE's own example video and **0/12** of ours. Our frames are
bimodal — dark IR background, bright mouse, SD 1.5× theirs — against a scorer that
hard-codes an input mean of 0.5175 and SD of 0.14, so a scalar gain cannot fix it: it
moves the mean and the SD together. This is a rig mismatch, not a software problem.
The fix is a light background, diffuse IR with no specular hot spot, and wider framing
so both ears stay in frame; then either GrimACE or PainFace should run with no
retraining. See `grimace/README.md`.

## 7. What is not in this repository

Regenerable bulk and unpublished material are excluded by `.gitignore`:

- `behavior_sync/painface_upload/*.mp4` — re-encoded clips, ~700 MB, rebuild with
  `grimace/prepare_painface_video.py`
- `behavior_sync/frames/` — extracted preview frames, rebuild with
  `behavior_sync/ccf_extract_face_frames.py`
- `ccf_registration/cache/*.mat` — full-res FOV images, rebuild with
  `ccf_cache_fov_images`
- `grimace/models/`, `grimace/grimace_output/`, `grimace/demo_output/`
- `final_deck/preview*/` — layout proofs, rebuild with `preview_deck.py`

The per-cell tables (`n18/cell_ccf_all.csv`, `melittin_ccf/cell_ccf_all.csv`) and the
slide decks hold unpublished data and are kept out of this public repository
deliberately. Every script here regenerates its own outputs from them once they are in
place.
