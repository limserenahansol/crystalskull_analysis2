# Analysis methods (Baseline vs Injection)

This document matches **`extract_analysis.m`** and is intended for **methods sections**, **figure legends**, and **collaborator handoffs**.

## Preprocessing

1. **Traces** — Full-frame-rate matrix `T` (frames × neurons) from EXTRACT `M_moco_frr.mat`.
2. **dF/F** — For each neuron, fluorescence normalized by an **8th-percentile rolling baseline** over a window of **15% of recording length** (minimum 300 frames).
3. **Z-score** — Each neuron’s dF/F trace is z-scored (mean and std across time).

## Event rate (per neuron)

- **Threshold** — `z_thresh` (default **2**).
- **Events** — **Upward** threshold crossings: transitions from `z ≤ 2` to `z > 2` on consecutive frames.
- **Rate** — Total crossings divided by **recording duration (s)** → **events/s per neuron**.
- **Population bar (fig2, panel 4)** — **Mean** per-neuron rate; error bars = **95% bootstrap CI** for that mean (`n_boot` resamples of **neurons**).
- **Inference** — **Two-sided Wilcoxon rank-sum** (`ranksum`) on **per-neuron** event rates (Baseline vs Injection).

## “Population synchrony” / fraction active (per 1 s bin)

**Not** pairwise spike synchrony or correlation structure. It is **per-time-bin recruitment**:

- **Bin width** — `time_bin` seconds (default **1** s).
- **Active in a bin** — Neuron active if **`z > z_thresh` in at least one frame** in that bin.
- **Fraction active** — Per bin: fraction of neurons meeting the criterion → a **time series** in `[0, 1]`.
- **Histogram (fig2, panel 3)** — Overlaid **PDF** histograms of per-bin values; vertical lines = **mean across bins**; title brackets = **95% bootstrap CI** for that mean (`n_boot` resamples of **bins**).
- **Wilcoxon** — Two-sided `ranksum` on **all baseline bins vs all injection bins** (same statistic as in `summary_stats.txt`).

### Why both **bootstrap CI** and **Wilcoxon**?

| Report | Role |
|--------|------|
| **Bootstrap 95% CI** | Uncertainty for the **mean** fraction across bins; **separation** of CIs between conditions is an intuitive precision/effect summary. |
| **Wilcoxon rank-sum** | Standard **p-value** for whether the two **distributions of bin values** differ (without assuming normality). |

**Caveat:** 1 s bins are **autocorrelated**; bootstrap and Wilcoxon assume **exchangeable/independent** samples. For strict inference, consider **block bootstrap** or hierarchical models.

## Figure numbering vs manuscript

Repo outputs **`fig2_population_activity.png`** for synchrony + mean event rate. **`fig5_PCA_dimensionality.png`** is **PCA only**. If a manuscript uses “Figure 5E/F” for synchrony and event-rate panels, note that **panel letters refer to the composite figure**, not the `fig5_*.png` filename.

## Other outputs

See **`summary_stats.txt`**: KS tests, Cohen’s *d*, PCA variance, etc., as implemented in the script.
