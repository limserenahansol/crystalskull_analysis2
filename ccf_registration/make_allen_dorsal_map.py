"""Build a top-down (dorsal-surface) Allen CCF parcellation for widefield / crystal-skull FOVs.

Reads the BrainGlobe allen_mouse_25um atlas already on this machine and writes
allen_dorsal_map.mat, which ccf_register_fov.m consumes.

Projection: for every (AP, ML) column, walk down the DV axis from the top and take the first
labelled voxel. That is literally "what you would see looking down on the brain", which is the
right target for a top-view 1p FOV. Layer-specific labels are collapsed to their parent area
(SSp-bfd1 -> SSp-bfd) because from above you only ever hit layer 1.

Output arrays are in CCF space: rows = AP (anterior at row 1), cols = ML (left hemisphere at
col 1), 25 um per pixel.
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy.io import savemat

ATLAS = Path.home() / ".brainglobe" / "allen_mouse_25um_v1.2"
OUT = Path(__file__).resolve().parent
OUT.mkdir(exist_ok=True)

# Bregma in CCFv3, the value used by cortex-lab/allenCCF and most widefield pipelines.
BREGMA_AP_UM = 5400.0
# Bregma-to-lambda along the skull in an adult mouse. Used only to seed the GUI; the user
# clicks the real thing.
BREGMA_LAMBDA_UM = 4200.0


def load_structures():
    df = pd.read_csv(ATLAS / "structures.csv")
    by_id = {}
    for rec in df.to_dict("records"):
        by_id[int(rec["id"])] = dict(
            acronym=rec["acronym"], name=rec["name"],
            path=[int(x) for x in str(rec["structure_id_path"]).strip("/").split("/") if x])
    return by_id


def collapse_layers(by_id):
    """Map every structure id to the nearest ancestor that is not a cortical layer."""
    layer_re = re.compile(r",\s*layer\s", re.I)
    out = {}
    for sid, s in by_id.items():
        target = sid
        # walk up while the current node is a layer subdivision
        while True:
            node = by_id.get(target)
            if node is None or not layer_re.search(node["name"]):
                break
            path = node["path"]
            if len(path) < 2:
                break
            target = path[-2]
        out[sid] = target
    return out


def main():
    meta = json.loads((ATLAS / "metadata.json").read_text())
    res = float(meta["resolution"][0])
    assert meta["orientation"] == "asr", f"unexpected orientation {meta['orientation']}"

    ann = tifffile.imread(ATLAS / "annotation.tiff")          # (AP, DV, ML)
    ref = tifffile.imread(ATLAS / "reference.tiff")
    n_ap, n_dv, n_ml = ann.shape

    by_id = load_structures()
    collapse = collapse_layers(by_id)

    # ---- first labelled voxel from the top, per (AP, ML) column
    labelled = ann > 0
    has_any = labelled.any(axis=1)                             # (AP, ML)
    first_dv = labelled.argmax(axis=1)                         # (AP, ML)
    ap_idx, ml_idx = np.meshgrid(np.arange(n_ap), np.arange(n_ml), indexing="ij")
    surface = ann[ap_idx, first_dv, ml_idx].astype(np.int64)
    surface[~has_any] = 0
    depth_um = first_dv.astype(np.float32) * res               # how deep the surface sits
    depth_um[~has_any] = np.nan

    # also grab the anatomical grey-scale surface, useful as a visual backdrop
    surf_ref = ref[ap_idx, first_dv, ml_idx].astype(np.float32)
    surf_ref[~has_any] = 0

    # ---- collapse layers
    flat = np.zeros_like(surface)
    for sid in np.unique(surface):
        if sid == 0:
            continue
        flat[surface == sid] = collapse.get(int(sid), int(sid))

    used = [int(x) for x in np.unique(flat) if x != 0]
    acronyms = [by_id[i]["acronym"] if i in by_id else f"id{i}" for i in used]
    names = [by_id[i]["name"] if i in by_id else f"id{i}" for i in used]

    # ---- boundaries (label changes between 4-neighbours)
    b = np.zeros(flat.shape, bool)
    b[:-1, :] |= flat[:-1, :] != flat[1:, :]
    b[:, :-1] |= flat[:, :-1] != flat[:, 1:]
    boundary = b & (flat > 0)

    # ---- landmarks in pixel coordinates of this map (row = AP, col = ML)
    bregma_row = BREGMA_AP_UM / res
    midline_col = n_ml / 2.0
    lambda_row = (BREGMA_AP_UM + BREGMA_LAMBDA_UM) / res

    def at(row, col):
        r, c = int(round(row)), int(round(col))
        sid = int(flat[r, c]) if 0 <= r < n_ap and 0 <= c < n_ml else 0
        return by_id[sid]["acronym"] if sid in by_id else "(outside)"

    print(f"atlas {meta['name']} v{meta['version']}  {res:.0f} um  shape {ann.shape}")
    print(f"dorsal map {flat.shape}  (AP x ML)  = {n_ap*res/1000:.1f} x {n_ml*res/1000:.1f} mm")
    print(f"{len(used)} surface structures")
    print(f"bregma  row {bregma_row:.0f}, col {midline_col:.0f}  -> surface structure {at(bregma_row, midline_col)}")
    print(f"lambda  row {lambda_row:.0f}, col {midline_col:.0f}  -> surface structure {at(lambda_row, midline_col)}")
    for probe, (dap, dml) in {"V1 (3.8P, 2.5L)": (3800, 2500),
                              "SSp-bfd (1.5P, 3.2L)": (1500, 3200),
                              "MOp (1.0A, 1.5L)": (-1000, 1500)}.items():
        print(f"  check {probe:24s} -> {at(bregma_row + dap / res, midline_col + dml / res)}")

    savemat(str(OUT / "allen_dorsal_map.mat"), {
        "label_img": flat.astype(np.uint32),
        "boundary": boundary,
        "surface_ref": surf_ref,
        "surface_depth_um": depth_um,
        "ids": np.array(used, dtype=np.uint32),
        "acronyms": np.array(acronyms, dtype=object),
        "names": np.array(names, dtype=object),
        "res_um": res,
        "bregma_rc": np.array([bregma_row, midline_col]),
        "lambda_rc": np.array([lambda_row, midline_col]),
        "bregma_lambda_um": BREGMA_LAMBDA_UM,
        "note": ("rows=AP (anterior first), cols=ML (left first), 25um/px. "
                 "Row/col are 0-based here; ccf_register_fov.m adds 1 on load."),
    }, do_compression=True)
    print("wrote", OUT / "allen_dorsal_map.mat")

    # ---- preview
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 8.5))
    # shade by how deep the surface sits: reads as a smooth dome, unlike the reference
    # intensity at the surface voxel which moires badly
    shade = np.where(np.isfinite(depth_um), depth_um, np.nan)
    ax.imshow(shade, cmap="bone_r", interpolation="nearest")
    ax.imshow(np.ma.masked_where(~boundary, boundary), cmap="autumn", interpolation="nearest")
    ax.plot(midline_col, bregma_row, "c+", ms=15, mew=2.5)
    ax.plot(midline_col, lambda_row, "cx", ms=13, mew=2.5)
    ax.annotate("bregma", (midline_col, bregma_row), color="c", xytext=(10, -2),
                textcoords="offset points", fontsize=10, fontweight="bold")
    ax.annotate("lambda", (midline_col, lambda_row), color="c", xytext=(10, -2),
                textcoords="offset points", fontsize=10, fontweight="bold")
    # label each hemisphere separately - a bilateral centroid lands on the midline and
    # stacks every acronym on top of the next
    for sid in used:
        acr = by_id[sid]["acronym"] if sid in by_id else str(sid)
        for half, sl in (("L", slice(0, int(midline_col))), ("R", slice(int(midline_col), n_ml))):
            m = np.zeros_like(flat, bool)
            m[:, sl] = flat[:, sl] == sid
            if m.sum() < 700:
                continue
            rr, cc = np.nonzero(m)
            # medoid-ish: the in-region pixel closest to the centroid, so the text never
            # lands outside a horseshoe-shaped area
            d = (rr - rr.mean()) ** 2 + (cc - cc.mean()) ** 2
            k = int(np.argmin(d))
            ax.text(cc[k], rr[k], acr, fontsize=5.5, color="black", ha="center", va="center",
                    bbox=dict(facecolor="yellow", alpha=.65, pad=.4, edgecolor="none"))
    ax.set_title("Allen CCF dorsal surface (25 um)  -  rows = AP, cols = ML\n"
                 "shading = depth of the first labelled voxel", fontsize=11)
    ax.set_xlabel("ML pixel")
    ax.set_ylabel("AP pixel")
    fig.tight_layout()
    fig.savefig(OUT / "allen_dorsal_map_preview.png", dpi=180)
    print("wrote", OUT / "allen_dorsal_map_preview.png")


if __name__ == "__main__":
    main()
