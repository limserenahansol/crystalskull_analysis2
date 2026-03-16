# Dataset Structure

This document describes the expected layout and contents of EXTRACT output files used by `extract_analysis.m`.

## Directory Layout

```
base_dir/
├── Baseline/
│   └── <bl_subfolder>/     e.g. mirror1path, m1fmousemirror2
│       ├── M_moco_frr.mat
│       ├── M_moco_ds_ext.mat
│       ├── M_summary.mat
│       └── cell_map.png
└── Injection/
    └── <inj_subfolder>/
        ├── M_moco_frr.mat
        ├── M_moco_ds_ext.mat
        ├── M_summary.mat
        └── cell_map.png
```

- `base_dir`: Root path (e.g. `\\server\crystal_skull\experiment_id` or `CS1014-1b`)
- `Baseline/` and `Injection/`: Condition folders
- `<bl_subfolder>` / `<inj_subfolder>`: Typically mirror/FOV identifiers (e.g. `mirror1path`)

## File Descriptions

### M_moco_frr.mat (HDF5 / v7.3)

Full-rate motion-corrected traces and spatial footprints.

| HDF5 path | Description | Dimensions |
|-----------|-------------|------------|
| `/T` | Fluorescence traces | `frames × neurons` (transposed to `neurons × frames` in MATLAB) |
| `/S/ir` | Row indices for sparse spatial matrix | 1D |
| `/S/jc` | Column pointers (CSC format) | `(neurons+1) × 1` |
| `/S/data` | Non-zero values for sparse spatial matrix | 1D |

- `T`: Raw fluorescence; used to compute dF/F.
- `S`: Sparse spatial weights (CSC). Reconstructed as `sparse(ir+1, col_idx, data, n_pixels, n_neurons)` and used to project per-neuron metrics onto the FOV.

### M_moco_ds_ext.mat (HDF5 / v7.3)

EXTRACT temporal components and metadata.

| HDF5 path | Description | Dimensions |
|-----------|-------------|------------|
| `/temporal_weights` | Temporal components | `neurons × time` (downsampled) |
| `/info/summary_image` | Mean/summary FOV image | `frame_px × frame_px` |

- `temporal_weights`: Used for optional analyses; main pipeline uses `T` from `M_moco_frr.mat`.
- `info/summary_image`: Displayed in Figure 1.

### M_summary.mat (HDF5 / v7.3)

| HDF5 path | Description |
|-----------|-------------|
| `/framerate` | Acquisition frame rate (Hz) |

### cell_map.png

Image of detected cell ROIs (e.g. from EXTRACT or manual segmentation). Displayed in Figure 1.

## Frame Size

The script assumes a square FOV of `frame_px × frame_px` (default `2304`). Adjust `frame_px` in the CONFIG section if your data differ.

## Data Format Notes

- EXTRACT v7.3 saves `.mat` files as HDF5. Use `h5read()` instead of `load()`.
- Sparse `S` uses 0-based row indices in `ir`; the script adds 1 for MATLAB indexing.
- `T` is read as `frames × neurons` and transposed to `neurons × frames` for downstream use.
