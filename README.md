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
- **`extract_cohort_delta_pipeline.m`** — **whole traces** (default: no time crop on BL or INJ; set `use_last_sec_inj` in `params` if you want). Edit CONFIG: two cohorts (root + mouse IDs), `mirror_subfolder`. Writes **`cohort_delta_per_mouse.csv`**, **`delta_group_comparison.txt`**, **`fig_cohort_delta_comparison.png`** under **`out_dir`** (default under your `2p` folder; change in CONFIG).
- **`extract_cohort_delta_pipeline_truncated.m`** — same cohort logic, but **baseline full** and **injection last 600 s** only (matches `extract_analysis2_truncated`). Uses a **separate `out_dir`** so it does not overwrite the whole-trace run.

For **full baseline vs last 10 min injection** (matched duration, different days), set in `params`:

- `skip_first_sec_bl = 0`, `use_last_sec_bl = inf` (baseline uncropped)
- `skip_first_sec_inj = 0`, `use_last_sec_inj = 600` (injection: last 10 min only)

Legacy symmetric `skip_first_sec` / `use_last_sec` still applies the **same** crop to both arms if you do not set the `*_bl` / `*_inj` fields. See **`extract_resolve_time_window_params.m`**.

**`extract_analysis2_truncated.m`** uses that asymmetric pattern by default.

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
