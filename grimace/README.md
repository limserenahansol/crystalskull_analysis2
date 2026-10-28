# Automated grimace scoring for crystal-skull recordings

The published algorithms stay in their own software. This folder **ingests their scores** and
aligns them to EXTRACT / poke times.

| Paper | Algorithm | Where it runs | Scale |
|---|---|---|---|
| McCoy et al., *Pain* 2024 | **PainFace** | cloud, [painface.net](http://painface.net) | 4 FAUs, 0–8 (black C57BL/6) |
| McCoy et al., *Pain* 2025 | PainFace, white coat | same | 5 FAUs, 0–10 |
| Sturman et al., *Lab Animal* 2026 | **GrimACE** | GrimACE box + app ([GitHub](https://github.com/ETHZ-INS/GrimACE_manuscript)) | 5 FAUs, 0–10 + top-view pose |

GrimACE was validated **after craniotomy with head implants**, which is the relevant case for
crystal skull. PainFace was trained on intact faces in a bright chamber; a skull-off / implant
FOV will score poorly unless the face is still fully visible and well lit.

Existing 13 Aug / 18 Aug face-poke movies **had no behavior camera**. Poke time is still the
cardboard blocker. Melittin / morphine face videos live on the share (do not copy the zips):

| Date | Assay | Behavior | Matched EXTRACT |
|---|---|---|---|
| 20260408 | melittin + morphine | `\\stg-tnr50a1.stanford.edu\crystal_skull\20260408\behavior` | `omer_hazon\morphine_20260408\1014m01b` / `1014m02w` |
| 20260218 | melittin only | `...\20260218\behavior` | `omer_hazon\meltinonly\CS1014-1b` / `CS1014-2w` (`1014m01b`/`1014m02w` same mice). `mouse1f` has no EXTRACT in the current CCF set |
| 20260219 | melittin only | `...\20260219\behavior` | `fcs_02041b`/`fcs_02042w` → `CS0204-1b`/`CS0204-2w`. `mouse2c` unmatched |

```bash
python grimace/build_behavior_manifest.py
python grimace/score_grimace_videos.py
```

Then in MATLAB: `grimace_run_matched`. Match table: `grimace/behavior_calcium_match.csv`.
Scores use official `rwitz/GrimACE` weights (face crop + 5-AU ViT). Those weights were trained
on GrimACE-box front IR video — inspect `grimace/grimace_output/figures/` before treating
MGS as validated.

20260218 and 20260219 are JPEG dumps. Listing the whole folder hangs, so scoring reads only
the first JPEG names `dir /b` returns (then stops). 20260408 is MP4 and is sampled across the file.
Figures: `grimace/grimace_output/figures/fig01`–`fig07`.

## 1. Score the video (do not skip this)

**PainFace** — 1080p color MP4, 30 fps, ~1500 lux, mouse facing the camera. Upload at
painface.net, use the black or white build that matches coat color, export the per-frame CSV.

**GrimACE** — official app is built for their box. We also run the published
`rwitz/GrimACE` weights offline (`score_grimace_videos.py`) on sampled mesoscope
frames. Treat those MGS values as provisional until QC crops look like a front face.

Prepare an MP4 if needed:

```bash
python grimace/prepare_painface_video.py path\to\face.avi
```

## 2. Align to imaging

```matlab
addpath('C:\Users\hsollim\behavior_task\crystal_skull_analysis\grimace')
addpath('C:\Users\hsollim\behavior_task\crystal_skull_analysis')

G = grimace_load("painface_export.csv", "Source", "painface");
M = grimace_metrics(G);                 % mean MGS, % high-grimace (5-8)
S = grimace_align_pokes(G, pokeTimesS); % peri-poke grimace, same windows as neural STA
grimace_plot(G, M, S, "grimace_out");
```

`pokeTimesS` is poke onset in **seconds from video t = 0**. If the face camera started with
the mesoscope, convert EXTRACT frames with `onset / fs + videoOffsetS`.

Session-level (melittin BL vs INJ, no pokes):

```matlab
Mb = grimace_metrics(grimace_load("bl.csv", "Source", "painface"));
Mi = grimace_metrics(grimace_load("inj.csv", "Source", "painface"));
```

## 3. Metrics (McCoy 2024)

On **MaxFAU** frames only (all 4 or 5 action units visible):

- mean MGS
- % time in low (0–2) / medium (3–4) / high (5–max) bins
- peri-poke: grimace in `resp_win` vs `base_win` (default 0.4–1.2 s vs −2–0 s)

Need ≥ 75 MaxFAU frames per 10 min (PainFace guideline) or the mean is not reported.

## CSV columns this loader accepts

PainFace (any subset; names are case-insensitive):

`frame`, `time_s` or `time`, `mgs` or `score` or `grimace`, `n_fau` or `fau_count`,
optional `orbital`, `nose`, `ear`/`ears`, `whisker`/`whiskers`

GrimACE: `frame`, optional `time_s`, per-AU columns `orbital`, `nose`, `cheek`, `ear`,
`whisker` (0/1/2 or `not rateable`). Total MGS = sum of rateable AUs.

Run `grimace_demo` once to write an example CSV and figures under `grimace/demo_output/`.
