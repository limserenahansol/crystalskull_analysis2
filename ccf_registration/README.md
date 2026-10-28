# Allen CCF registration for the crystal-skull top-view FOV

Registers a top-down 1p crystal-skull FOV to the Allen CCFv3 dorsal cortical surface, so
every EXTRACT cell gets an Allen area label and a bregma-referenced AP/ML coordinate.

Built from the `allen_mouse_25um` volume already cached on this machine at
`~/.brainglobe/allen_mouse_25um_v1.2` — nothing is downloaded.

## Files

| File | What it does |
|---|---|
| `make_allen_dorsal_map.py` | One-time. Projects the CCF annotation volume top-down and writes `allen_dorsal_map.mat` + a preview PNG. Already run. |
| `ccf_register_fov.m` | Interactive registration. Returns and saves the transform. |
| `ccf_apply_to_cells.m` | Cell centroids → atlas coords, AP/ML in mm, Allen area names. |
| `ccf_area_centroid.m` | Atlas coordinates of named areas. Call with no arguments to list all 58 surface areas. |
| `test_ccf_roundtrip.m` | Synthetic check. Currently passes: transform recovered to <1 µm, 100% area agreement. |
| `ccf_register_facepoke_n18.m` | Interactive batch: all 18 face-poke FOVs → CCF. Resume-safe. |
| `ccf_label_facepoke_n18.m` | After transforms exist: Allen area for every EXTRACT cell. |
| `ccf_n18_atlas_gallery.m` | QC: warp registered FOVs onto the dorsal map. |
| `ccf_cache_fov_images.m` | Cache EXTRACT summaries for poke + melittin FOVs. |
| `ccf_autodetect_midline.m` | Find the sagittal sinus automatically. |
| `ccf_register_assisted.m` | Fewer-click registration once the midline is cached. |

## The atlas map

`allen_dorsal_map.mat` is 528 × 456 px at 25 µm (13.2 × 11.4 mm), rows = AP with anterior
first, cols = ML with left first. For each (AP, ML) column it takes the first labelled voxel
walking down from the top — literally what you see looking down on the brain — and collapses
cortical layers to their parent area, since from above you only ever hit layer 1. Bregma sits
at row 216, col 228.

Sanity checks that passed when it was built: VISp at 3.8 mm posterior / 2.5 mm lateral → `VISp`,
SSp-bfd at 1.5P / 3.2L → `SSp-bfd`, MOp at 1.0A / 1.5L → `MOp`.

Look at `allen_dorsal_map_preview.png` before you register anything, so you know what you are
aiming at.

## Face-poke n18 FOVs (13 Aug + 18 Aug)

Do **not** warp the 18-panel montage as one image. Each recording is its own FOV.
The batch loads the EXTRACT `summary_image` (same pixels as `example_fov.png`, P left / A right).

```matlab
addpath('C:\Users\hsollim\behavior_task\crystal_skull_analysis\ccf_registration')
addpath('C:\Users\hsollim\behavior_task\crystal_skull_analysis')

% Interactive. One cpselect window per recording. Resume-safe.
ccf_register_facepoke_n18                    % all 18
ccf_register_facepoke_n18("Only", "m1f")     % one mouse
ccf_register_facepoke_n18("Mode", "bregma")  % 2 clicks instead of 4-6

ccf_n18_atlas_gallery
ccf_label_facepoke_n18
```

Transforms and QC land in `ccf_registration/n18/<label>/`. The atlas in cpselect is
rotated to **P left, A right** so it matches the lab maps (left hemisphere toward the top).
Cyan cross = bregma, yellow = lambda. The crystal-skull window usually covers MOs / MOp /
SSp / RSP / VIS — not the olfactory bulbs — so click the midline sinus, not MOB.

A free check after labelling: face-poke excitation should pile up in mouth / nose / barrel
cortex (`SSp-n`, `SSp-m`, `SSp-bfd`). If the peak is in MOs or VISp, redo that FOV.

## Superimposed figures (all FOVs on CCF)

After transforms exist:

```matlab
ccf_plot_superimposed
```

Writes `ccf_registration/figures/`. Start with `fig_ccf_poke_superimposed.png` (coverage + mean poke response). Crystal-skull windows cover most of dorsal cortex, so outlines sit near the atlas edge; the useful maps are coverage and the red/blue response overlay.

## Three ways to register, worst to best

### 1. `mode="bregma"` — two clicks

```matlab
T = ccf_register_fov(summaryImage);              % click bregma, then lambda
```

Fits a similarity transform (rotation, uniform scale, translation). Fastest, and what most
widefield papers do. Accuracy is limited by how well you can see the sutures: **typically
200–400 µm**, which is fine for "SSp vs MOs" and not fine for "SSp-bfd vs SSp-m".

Two points exactly determine a similarity transform, so the reported residual is 0 by
construction and tells you nothing. Do not read it as accuracy.

### 2. `mode="landmarks"` — cpselect, 4–6 points

```matlab
T = ccf_register_fov(summaryImage, "mode", "landmarks", "transform", "affine");
```

Opens MATLAB's control-point tool with your FOV and the atlas side by side. **For a top view
the vasculature beats the sutures**: the superior sagittal sinus gives you the midline and
therefore rotation and ML = 0, and the two junctions where it meets the transverse sinuses
pin the posterior end. Those are the highest-contrast features in a crystal-skull image and
they are anatomically tight.

Use `affine` (needs ≥3 points) rather than `projective` unless the FOV genuinely has
perspective distortion — with an overfitted transform the borders will look beautiful and be
wrong.

With more points than the transform needs, the printed RMS residual is a real error estimate.

### 3. `mode="points"` — functional registration, the accurate one

Skull and sinus landmarks are proxies. The cortex itself moves relative to them between
animals by a few hundred microns. If you can spare one short mapping session, use sensory
maps as the anchors:

| Stimulus | Anchor area |
|---|---|
| full-field flash | `VISp` |
| whisker deflection | `SSp-bfd` |
| forelimb / hindlimb tap | `SSp-ul` / `SSp-ll` |
| tone | `AUDp` |

```matlab
P = ccf_area_centroid(["VISp","SSp-bfd","SSp-ul"], "hemisphere","R");
T = ccf_register_fov(img, "mode","points", "transform","affine", ...
                     "fovPoints", myPeaks, "ccfPoints", [P.x, P.y]);
```

This typically halves the error versus skull landmarks, and it is the only version that is
defensible if a reviewer asks how you know a cell is in SSp-bfd rather than SSp-m.

**A free check on the data you already have:** the face-poke evoked excitation should peak in
SSp mouth/nose/barrel. Register with bregma/lambda, then run `ccf_apply_to_cells` on the
excited cells and see where the peak lands. If it comes out in MOs or VISp, the registration
is wrong — and if it lands in SSp, you have validated the registration for free on all 18
recordings.

## Labelling cells

```matlab
tbl = ccf_apply_to_cells(T, cellXY);
% cell | atlas_row | atlas_col | AP_mm | ML_mm | hemisphere | area_acronym | area_name
summary(categorical(tbl.area_acronym))
```

AP is positive anterior, ML positive right, both in mm from bregma.

## Two things to watch

**mirror2 is flipped.** The analysis scripts flip mirror2 left-right (`flipLR = contains(rec.mirror,
"mirror2")`). Register each mirror separately and pass `"flipLR", true` for mirror2, or make
sure your cell coordinates are already in the flipped frame. Getting this wrong mirrors the
hemispheres and every area label is wrong in a way that looks plausible.

**Check the implied scale.** `ccf_register_fov` prints µm per FOV pixel. If that does not match
what you expect from the objective and sensor, the fit is wrong regardless of how good the
overlay looks.

## Honest limits

The CCF dorsal projection is a flat top-down view of a curved surface. Near the lateral edges,
where cortex curves away steeply, a 2D affine transform cannot be right — errors grow toward
the periphery of the FOV. For medial and dorsal areas (MOs, MOp, SSp, RSP, VISp) this is not a
practical problem. For the lateral edge (AUDp, TEa, VISpor) treat the labels as approximate.
A proper fix is a cortical flat-map, which is a much bigger job and probably not worth it here.
