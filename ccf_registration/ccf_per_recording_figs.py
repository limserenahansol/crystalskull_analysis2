"""Per-recording superimposed figures: raw FOV + Allen borders + responsive cells.

One PNG per recording in figures_v2/per_recording/, plus a contact sheet per paradigm.
Cells are drawn at their FOV pixel position, so this shows the raw data and the registration
together - if a panel looks wrong, that recording's transform is wrong.
"""
import csv
import collections
import os

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import affine_transform, binary_dilation

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures_v2", "per_recording")
os.makedirs(OUT, exist_ok=True)

A = sio.loadmat(os.path.join(HERE, "allen_dorsal_map.mat"))
BOUND = A["boundary"].astype(bool)
RED, BLUE = "#E03131", "#1C7ED6"


def stretch(img):
    v = img[np.isfinite(img)]
    lo, hi = np.percentile(v, [1, 99.5])
    return np.clip((img - lo) / max(hi - lo, 1e-9), 0, 1)


# The saved transform is a MATLAB affinetform2d object, which scipy cannot read, so the
# matrices are exported to transform_matrices.csv by a one-line MATLAB step.
TFORM = {}
with open(os.path.join(HERE, "transforms", "transform_matrices.csv"), newline="") as fh:
    for r in csv.DictReader(fh):
        TFORM[r["label"]] = np.array([
            [float(r["a11"]), float(r["a12"]), float(r["a13"])],
            [float(r["a21"]), float(r["a22"]), float(r["a23"])],
            [0.0, 0.0, 1.0]])


def warp_boundary_to_fov(M, shape):
    """Draw the atlas boundaries in FOV pixel space.

    scipy.ndimage.affine_transform samples output -> input. Here the output grid is the FOV
    and the input image is the atlas, so the mapping it needs is FOV -> atlas, which is the
    FORWARD transform M. Passing the inverse (the intuitive but wrong choice) pushes the
    borders almost entirely off the frame.
    """
    mat = np.array([[M[1, 1], M[1, 0]], [M[0, 1], M[0, 0]]])
    off = np.array([M[1, 2], M[0, 2]])
    bw = affine_transform(BOUND.astype(float), mat, offset=off,
                          output_shape=shape, order=1) > 0.15
    # the atlas is ~4x coarser than the FOV, so thicken slightly to stay visible
    return binary_dilation(bw, np.ones((3, 3)))


def load_cells(path, keycol=None):
    by = collections.defaultdict(list)
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            by[r["label"]].append(r)
    return by


poke_cells = load_cells(os.path.join(HERE, "n18", "cell_ccf_all.csv"))
mel_cells = load_cells(os.path.join(HERE, "melittin_ccf", "cell_ccf_all.csv"))

# face-poke CSV labels have no "poke_" prefix; transform/cache files do
jobs = []
for lab, rows in poke_cells.items():
    jobs.append(("poke_" + lab, lab, rows, "Face poke"))
for lab, rows in mel_cells.items():
    coh = "Melittin+morphine" if rows[0]["cohort"] == "melittin_morphine" else "Melittin"
    jobs.append((lab, lab, rows, coh))

print(f"{len(jobs)} recordings")

panels = collections.defaultdict(list)
for cache_lab, csv_lab, rows, cond in sorted(jobs):
    cpath = os.path.join(HERE, "cache", cache_lab + ".mat")
    if not os.path.exists(cpath):
        print("  skip (no cached image):", cache_lab)
        continue
    img = sio.loadmat(cpath)["img"].astype(float)
    if "mirror2" in cache_lab:
        img = np.fliplr(img)
    if cache_lab not in TFORM:
        print("  skip (no transform matrix):", cache_lab)
        continue
    M = TFORM[cache_lab]
    bw = warp_boundary_to_fov(M, img.shape)

    # cells: recompute FOV position from the atlas position via the inverse transform
    Minv = np.linalg.inv(M)
    ar = np.array([float(r["atlas_row"]) for r in rows])
    ac = np.array([float(r["atlas_col"]) for r in rows])
    ones = np.ones_like(ar)
    fx = Minv[0, 0] * ac + Minv[0, 1] * ar + Minv[0, 2] * ones
    fy = Minv[1, 0] * ac + Minv[1, 1] * ar + Minv[1, 2] * ones
    exc = np.array([r["is_exc"] == "1" for r in rows])
    inh = np.array([r["is_inh"] == "1" for r in rows])

    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    ax.imshow(stretch(img), cmap="gray")
    ax.imshow(np.ma.masked_where(~bw, bw), cmap="autumn_r", alpha=.55,
              interpolation="nearest")
    ax.scatter(fx[inh], fy[inh], s=2.2, c=BLUE, alpha=.42, linewidths=0)
    ax.scatter(fx[exc], fy[exc], s=2.6, c=RED, alpha=.62, linewidths=0)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{csv_lab}\n{cond}   {len(rows):,} cells   "
                 f"exc {100*exc.mean():.1f}%  inh {100*inh.mean():.1f}%",
                 fontsize=9)
    fig.tight_layout()
    dst = os.path.join(OUT, cache_lab + ".png")
    fig.savefig(dst, dpi=110)
    plt.close(fig)
    panels[cond].append((csv_lab, dst))

# contact sheets
for cond, items in panels.items():
    n = len(items)
    ncol = 6 if n > 8 else 4
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.7 * ncol, 2.8 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax in axes:
        ax.axis("off")
    for ax, (lab, path) in zip(axes, sorted(items)):
        ax.imshow(plt.imread(path))
        ax.set_title(lab, fontsize=5.5)
        ax.axis("off")
    fig.suptitle(f"{cond} - all recordings, raw FOV + Allen borders + responsive cells "
                 f"(red = excited/up, blue = inhibited/down)", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, .96])
    slug = cond.lower().replace(" ", "_").replace("+", "_")
    fig.savefig(os.path.join(HERE, "figures_v2", f"fig4_per_recording_{slug}.png"), dpi=140)
    plt.close(fig)
    print(f"  contact sheet: {cond} ({n} panels)")

print("->", OUT)
