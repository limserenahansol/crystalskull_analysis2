# Crystal Skull Analysis

Population-level analysis of **EXTRACT** calcium imaging outputs for comparing **Baseline vs Injection** conditions in crystal skull / mesoscope whole dorsal cortex recordings.

## Overview

This repository provides a MATLAB pipeline that:

- Loads EXTRACT outputs (HDF5/v7.3 `.mat` files)
- Computes dF/F (8th-percentile rolling baseline) and z-scores
- Derives event rates (z > 2 threshold crossings per second)
- Produces **7 publication-ready figures** and a **summary statistics** file

## Requirements

- **MATLAB** (R2016b or later recommended)
- No additional toolboxes required (uses base MATLAB + `Statistics and Machine Learning Toolbox` for `ranksum`, `bootstrp`, `pca`, `kstest2`, `skewness`)

## Quick Start

1. **Clone the repository**
   ```bash
   git clone https://github.com/limserenahansol/crystalskull_analysis.git
   cd crystalskull_analysis
   ```

2. **Edit the CONFIG section** at the top of `extract_analysis.m`:
   - `base_dir` — root path to your experiment (e.g. `\\server\crystal_skull\experiment_id`)
   - `bl_subfolder` — subfolder under `Baseline/` (e.g. `mirror1path`)
   - `inj_subfolder` — subfolder under `Injection/` (e.g. `mirror1path`)
   - `out_dir` — where figures and `summary_stats.txt` are saved

3. **Run the script** in MATLAB:
   ```matlab
   run('extract_analysis.m')
   ```

## Dataset Structure

Your data must follow the EXTRACT output layout. See **[docs/DATASET_STRUCTURE.md](docs/DATASET_STRUCTURE.md)** for a detailed description.

**Typical layout:**

```
base_dir/
├── Baseline/
│   └── mirror1path/   (or m1fmousemirror2, etc.)
│       ├── M_moco_frr.mat      # Full-rate traces (T) + spatial weights (S)
│       ├── M_moco_ds_ext.mat   # Temporal weights, summary image
│       ├── M_summary.mat       # Framerate
│       └── cell_map.png        # Cell map image
└── Injection/
    └── mirror1path/
        └── (same files)
```

## Outputs

| File | Description |
|------|-------------|
| `fig1_cellmaps_FOV.png` | Summary images, cell maps, spatial mean activity |
| `fig2_population_activity.png` | Population timecourses, event rate, dF/F, neuron counts |
| `fig3_neuron_distributions.png` | CDFs of mean dF/F, event rate, std, skewness |
| `fig4_network_correlations.png` | Pairwise correlation matrices and distributions |
| `fig5_PCA_dimensionality.png` | Scree plot, neural trajectories |
| `fig6_activity_heatmaps.png` | Z-scored activity heatmaps (sorted by peak) |
| `fig7_spatial_event_rate.png` | Event rate projected onto FOV |
| `summary_stats.txt` | Population metrics, medians, statistical tests |

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

- **Wilcoxon rank-sum** for mean dF/F, event rate, pairwise correlations, fraction active
- **Kolmogorov–Smirnov** for distribution comparisons
- **Bootstrap 95% CI** for population-level metrics
- **Cohen's d** for effect sizes

## License

See [LICENSE](LICENSE) if provided.

## Citation

If you use this pipeline, please cite the EXTRACT algorithm and this repository.
