# Crystal Skull Analysis (`crystalskull_analysis2`)

Population-level analysis of **EXTRACT** calcium imaging outputs for comparing **Baseline vs Injection** conditions in crystal skull / mesoscope whole dorsal cortex recordings.

**GitHub:** [github.com/limserenahansol/crystalskull_analysis2](https://github.com/limserenahansol/crystalskull_analysis2)

## Overview

This repository provides a MATLAB pipeline that:

- Loads EXTRACT outputs (HDF5/v7.3 `.mat` files)
- Computes dF/F (8th-percentile rolling baseline) and z-scores
- Derives event rates (z > 2 upward crossings per second)
- Computes **per–1 s bin fraction active** (“population synchrony” / recruitment; see **docs/METHODS.md**)
- Produces **7 publication-ready figures** and **`summary_stats.txt`** (bootstrap CIs, Wilcoxon, KS, Cohen’s *d*)
- Supports **batch mode**: scan a parent folder for multiple experiment IDs (e.g. `CS1014-1b`, `CS0204-1b`) and write **`out_root/<experiment_id>/`** for each

## Documentation (methods & legends)

| Document | Purpose |
|----------|---------|
| [**docs/METHODS.md**](docs/METHODS.md) | Definitions: dF/F, z-score, event rate, fraction active; bootstrap vs Wilcoxon; figure-numbering note |
| [**docs/DATASET_STRUCTURE.md**](docs/DATASET_STRUCTURE.md) | EXTRACT file layout and HDF5 fields |
| [**docs/FIGURE_LEGENDS.md**](docs/FIGURE_LEGENDS.md) | Copy-paste figure legend text |
| [**PUSH_INSTRUCTIONS.md**](PUSH_INSTRUCTIONS.md) | Push to `crystalskull_analysis2` |

### Why report **bootstrap 95% CI** and **Wilcoxon rank-sum**?

- **Bootstrap CI** — Uncertainty for a **summary** (e.g. mean fraction active across **1 s bins**, or mean event rate across **neurons**). Non-overlapping CIs are a clear **effect + precision** visual.
- **Wilcoxon** (`ranksum`) — Standard **two-sample test** on the **observations** (bins or neurons) without assuming normality; gives a **p-value** reviewers often expect.

They are **complementary** (estimation vs distribution-level test). Details and caveats (e.g. temporal autocorrelation of bins): **docs/METHODS.md**.

## Requirements

- **MATLAB** (R2016b or later recommended)
- **Statistics and Machine Learning Toolbox** (`ranksum`, `bootstrp`, `pca`, `kstest2`, `skewness`)

## Quick Start

1. **Clone the repository**
   ```bash
   git clone https://github.com/limserenahansol/crystalskull_analysis2.git
   cd crystalskull_analysis2
   ```

2. **Edit the CONFIG section** at the top of `extract_analysis.m`:
   - **`batch_mode`** — `true` = auto-discover experiments under **`root_hazon`**; `false` = only **`exp_list_manual{1}`**
   - **`root_hazon`** — parent folder containing one subfolder per experiment (each with `Baseline/` and `Injection/`)
   - **`bl_subfolder`** / **`inj_subfolder`** — mirror path under each condition (e.g. `mirror2path`)
   - **`out_root`** — defaults to **`extract_figures/`** next to the script; outputs go to **`out_root/<experiment_id>/`**
   - **`exp_list_filter`** — optional cell array to limit which folder names run (default `{}` = all valid)

3. **Run** in MATLAB:
   ```matlab
   run('extract_analysis.m')
   ```

Parameters: `z_thresh`, `time_bin`, `n_boot`, `n_corr_sample`, etc. Incomplete experiments (missing EXTRACT files) are **skipped**; the batch continues.

**Memory (large mesoscope movies):** If MATLAB reports **Out of memory** during dF/F, lower **`dff_neuron_block`** (e.g. `512` or `256`). **`use_single_trace = true`** (default) stores traces and dF/F in **single** precision (~half the RAM of double). Set **`use_single_trace = false`** only if you need full double precision.

## Dataset Structure

Your data must follow the EXTRACT output layout. See **[docs/DATASET_STRUCTURE.md](docs/DATASET_STRUCTURE.md)** for a detailed description.

**Typical layout (each experiment folder under `root_hazon`):**

```
<experiment_id>/   e.g. CS1014-1b
├── Baseline/
│   └── mirror1path/   (or mirror2path, m1fmousemirror2, etc.)
│       ├── M_moco_frr.mat      # Full-rate traces (T) + spatial weights (S)
│       ├── M_moco_ds_ext.mat   # Temporal weights, summary image
│       ├── M_summary.mat       # Framerate
│       └── cell_map.png        # Cell map image
└── Injection/
    └── mirror1path/
        └── (same files)
```

## Cross-day cohort comparison (Δ = injection − baseline)

When **baseline differs across days**, compare **within-mouse change** (injection minus baseline on the same day) instead of raw cross-day baselines.

- **`extract_session_metrics.m`** — metrics for one mirror folder (Baseline or Injection), same preprocessing as `extract_analysis2_core` (dF/F, z-score, event rate, 1 s population bins, sampled pairwise correlations, PCA variance).
- **`extract_cohort_delta_pipeline.m`** — **whole traces** (default: no time crop on BL or INJ; set `use_last_sec_inj` in `params` if you want). Edit CONFIG: two cohorts (root + mouse IDs), `mirror_subfolder`. Writes **`cohort_delta_per_mouse.csv`** (includes **`pair_idx`** when paired), **`delta_group_comparison.txt`**, **`fig_cohort_delta_comparison.png`** under **`out_dir`**.
- **`extract_cohort_delta_pipeline_truncated.m`** — same cohort logic with **`match_inj_last_to_baseline_dur`**: **injection = last N s** where **N = baseline duration** for that mouse (after any baseline time crops). Default baseline = full file. Separate **`out_dir`**.
- **`paired_cohort_compare`** (in both cohort scripts, default **true**) — treats **`cohort_a_mice{i}`** and **`cohort_b_mice{i}`** as the **same biological mouse** on different study days. Group tests use **paired *t*-test** (`ttest`) and **Wilcoxon signed-rank** (`signrank`) on those Δ vectors; the figure draws **gray lines** between pair members. Set **`paired_cohort_compare = false`** for unrelated mice per cohort (unpaired **ranksum** only).
- **`extract_cohort_delta_from_session_pair.m`** — one mouse’s Δ from BL+INJ paths (shared by cohort pipelines).
- **`extract_cohort_delta_pair_stats.m`** — paired vs unpaired *p*-values.

Per-mouse **fig1–fig7** (same style as `extract_analysis2_core`) come from **`extract_analysis2_truncated_batch`** / **`extract_analysis2_core`** on each session folder; the cohort scripts only summarize **scalar Δ** across mice.

**`extract_cohort_neuron_delta_figures.m`** — pool **all injection neurons** across mice in each cohort and compare **Δ = injection − baseline session mean** (per neuron) with **KS** and **rank-sum** (large *N*). CONFIG **`time_window_mode`**: **`'matched'`** (injection last *N* s = baseline duration, like **`extract_analysis2_truncated_batch`**) or **`'whole'`** (full BL + full INJ, like whole-trace batch). Outputs go under **`fullfile(out_dir_parent, time_window_mode)`** unless **`out_dir_override`** is set. Run twice with each mode for both versions. Outputs **`fig_cohort_delta_population.png`** (Δ population timecourses, pooled synchrony bins), **`fig_cohort_delta_neuron_cdf.png`** (fig3-style ECDFs for Δ dF/F, event rate, std, skew), **`fig_cohort_delta_heatmaps.png`** (ΔZ: injection *z* minus baseline population-mean *z* over time), **`fig_cohort_delta_spatial.png`** (cohort-mean Δ event rate and Δ mean dF/F maps), and **`cohort_neuron_delta_stats.txt`**. Neurons cluster within mice (non-independent); treat *p*-values as exploratory unless you add mixed models or mouse-level bootstrap.

For **injection length = baseline length** (per mouse), set **`match_inj_last_to_baseline_dur = true`** in `params` (used by `extract_analysis2_core`, `extract_analysis2_truncated*`, and `extract_cohort_delta_from_session_pair`). After any baseline time crop, **N = baseline duration in seconds**; the injection trace uses the **last N seconds** (at injection frame rate). If the baseline file is 900 s, injection uses the last 900 s; if 800 s, last 800 s.

To cap baseline to the **first** 600 s and still match injection to that cap, set e.g. **`use_first_sec_bl = 600`** and keep **`match_inj_last_to_baseline_dur = true`** (injection last 600 s).

To fix injection length manually instead, set **`match_inj_last_to_baseline_dur = false`** and set **`use_last_sec_inj`** as before.

Legacy symmetric `skip_first_sec` / `use_last_sec` still applies the **same** crop to both arms if you do not set the `*_bl` / `*_inj` fields. See **`extract_resolve_time_window_params.m`**.

**`extract_analysis2_truncated.m`** / **`extract_analysis2_truncated_batch.m`** default to **`match_inj_last_to_baseline_dur = true`** with a **full** baseline trace so fig2 axes match per mouse (longer baseline ⇒ longer matched injection tail).

## Stimulus-triggered (facial-poke) analysis

Separate from the Baseline-vs-Injection pipeline, **`extract_stim_triggered_analysis.m`** (+ driver **`extract_stim_triggered_batch.m`**) aligns EXTRACT calcium activity to **facial-poke stimulus times** detected from the mean-frame-intensity trace (Omer's `estimated_stim_timing.mat`).

**Stim file (`estimated_stim_timing.mat`)** — one entry per movie:

| Variable | Meaning |
|----------|---------|
| `fnames{i}` | movie file name (e.g. `m201_d250813_..._mirror2.h5`) |
| `peak_idx{i}` | detected positive intensity-peak frame indices |
| `peak_pre_idx{i}` | minimum frame index ~1 s **before** each peak (poke-onset estimate) |
| `peak_post_idx{i}` | minimum frame index ~1 s **after** each peak (poke-offset estimate) |

Frame indices are in the **original movie frame space** (33 Hz). The code assumes the EXTRACT full-rate traces (`/T` in `*moco_frr.mat`) have one row per movie frame at the same rate; if the trace frame count differs, pass `params.n_frames_movie` to rescale, otherwise out-of-range pokes are clipped with a warning.

**Onset anchor (`params.event_field`)** — which index is `t = 0`: `'peak'` (default), `'pre'`, or `'post'`. The mean `pre`/`post` offsets are drawn as dotted markers on the PSTHs regardless, so all of Omer's timing info is shown. Run the batch three times (changing `event_field`) if you want all anchors.

**Run:** edit the CONFIG block of `extract_stim_triggered_batch.m` (`stim_mat_file`, `ext_roots`, `out_root`, `params`) then `run('extract_stim_triggered_batch.m')`. Movies are auto-mapped to EXTRACT folders by matching the movie basename against `*moco_frr.mat` files; use the optional `manual_map` for exceptions.

**Per-movie outputs** (in `out_root/<movie>/`):

| File | Description |
|------|-------------|
| `fig1_stim_population.png` | Population PSTH: mean z-score and fraction-active (mean ± SEM across pokes), poke onset/offset markers, resp-vs-base *p* |
| `fig2_neuron_sta_heatmap.png` | Per-neuron stimulus-triggered average (sorted by response amplitude) + group-mean STA (excited/inhibited/all) |
| `fig3_responsive_cells.png` | Excited/inhibited/n.s. counts (paired pre-vs-post signed-rank or *t*-test), amplitude distribution, *p*-value CDF |
| `fig4_spatial_response.png` | Response amplitude (post − pre, z) projected onto the FOV |
| `fig5_trial_consistency.png` | Per-poke population raster + inter-poke-interval histogram |
| `stim_summary.txt` | Neurons, pokes used, ITI, window params, responsive counts, population *p* |

**Cohort outputs** (in `out_root/`): `cohort_stim_summary.csv` (per-movie scalars) and `cohort_stim_population.png` (per-movie + cohort-mean pooled PSTH).

**Response test:** for each neuron, mean z in the **response window** (`resp_win`, default `[0 2] s`) vs the **baseline window** (`base_win`, default `[-2 0] s`) across pokes, via Wilcoxon **signed-rank** (default) or paired *t*-test; excited = *p* < α & amp > 0, inhibited = *p* < α & amp < 0. Neurons are non-independent within a movie — treat single-movie *p*-values as exploratory.

### Crystal-skull EXTRACT (`*_ext.mat`) layout

For 13 Aug / 18 Aug movies, cell extraction is a single `*_moco_ds_pre_ext.mat` (no separate `moco_frr`). Use **`extract_stim_ext_analysis.m`** / **`extract_stim_ext_batch.m`**. Traces are `/temporal_weights` at 33 Hz. Spatial maps reconstruct EXTRACT `spatial_weights` from the ndSparse CSC under `/#refs#` (`extract_load_ext_spatial.m`). **13 Aug `peak_idx`** = positive brightness peak; **18 Aug `peak_idx`** = single-frame intensity dip. Do not copy 13 Aug mW values onto 18 Aug (`pwr NA`).

Pooled n=18 figures: **`write_n18_clean_figs.m`**, **`extract_stim_pool_days.m`**, **`extract_stim_ap_takehome.m`**. PPTs: `build_crystal_skull_combined_pptx.py`, `build_example_neurons_pptx.py`.

### Example neurons (traces + FOV maps)

**`extract_stim_example_neurons.m`** picks 4 cells per recording (strongest significant excited and inhibited in the anterior and posterior halves), writes EXTRACT traces with poke times, peri-poke STA, and contrast-stretched crystal-skull FOV maps (P left, A right; mirror2 flipped). Outputs: `stim_triggered_figures/combined_0813_0818/example_neurons/` including `CrystalSkull_face_poke_example_traces_maps_n18.pptx`.

## Outputs

| File | Description |
|------|-------------|
| `fig1_cellmaps_FOV.png` | Summary images, cell maps, spatial mean activity |
| `fig2_population_activity.png` | Population timecourses; **synchrony histogram + bootstrap CI + Wilcoxon *p***; mean event rate + CI + *p*; dF/F; neuron counts |
| `fig3_neuron_distributions.png` | CDFs of mean dF/F, event rate, std, skewness |
| `fig4_network_correlations.png` | Pairwise correlation matrices and distributions |
| `fig5_PCA_dimensionality.png` | Scree plot, neural trajectories |
| `fig6_activity_heatmaps.png` | Z-scored activity heatmaps (sorted by peak) |
| `fig7_spatial_event_rate.png` | Event rate projected onto FOV |
| `summary_stats.txt` | Population metrics, medians, statistical tests |

With **batch mode**, each row above lives under **`out_root/<experiment_id>/`**.

**Manuscript vs file names:** Synchrony and population mean event rate are in **`fig2_population_activity.png`**, not in `fig5_PCA_dimensionality.png` (PCA only). If your report uses “Figure 5E/F” for those panels, state that **panel letters refer to the composite manuscript figure**.

## Analysis Pipeline

1. **Load** full-rate traces `T` (frames × neurons) and sparse spatial weights `S` from `M_moco_frr.mat`
2. **Compute dF/F** using 8th-percentile rolling baseline (15% window)
3. **Z-score** each neuron
4. **Event rate** = count of z > 2 upward crossings per second
5. **Population vectors** = mean z-score and fraction active per 1-s bin
6. **Pairwise correlations** (500 neurons sampled)
7. **PCA** on binned population activity
8. **Spatial maps** = project mean dF/F and event rate onto FOV via `S`

## Statistical Tests

- **Wilcoxon rank-sum** — mean dF/F, event rate, pairwise correlations, **fraction active (per-bin series)**; compares Baseline vs Injection **distributions** of observations.
- **Bootstrap 95% CI** — mean fraction active across **bins**; mean event rate and mean dF/F across **neurons** (`n_boot` resamples, percentile method).
- **Kolmogorov–Smirnov** — selected per-neuron distribution comparisons (see `summary_stats.txt` and fig3).
- **Cohen's d** — effect sizes where reported in `summary_stats.txt`.

## License

See [LICENSE](LICENSE) if provided.

## Citation

If you use this pipeline, please cite the EXTRACT algorithm and this repository.
