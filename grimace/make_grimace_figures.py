"""Figures for GrimACE on crystal-skull face videos.

Reads per-session CSVs and QC crops. Does not copy the share videos.
Frame scores are not independent replicates — mouse means are the plotted unit.
"""
from __future__ import annotations

import csv
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "grimace_output"
QC = OUT / "qc_crops"
FIG = OUT / "figures"
AU = ["orbital", "nose", "cheek", "ear", "whisker"]
AU_LABEL = {
    "orbital": "Orbital",
    "nose": "Nose",
    "cheek": "Cheek",
    "ear": "Ear",
    "whisker": "Whisker",
}


def load_frames() -> pd.DataFrame:
    files = sorted(p for p in OUT.glob("*_grimace.csv") if p.name != "grimace_all_frames.csv")
    if not files:
        raise FileNotFoundError(f"no session CSVs in {OUT}")
    frames = [pd.read_csv(p) for p in files]
    df = pd.concat(frames, ignore_index=True)
    for c in AU + ["mgs", "n_fau", "face_conf", "frame_quality"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["condition"] = df.apply(_condition, axis=1)
    df["paired_epoch"] = df["epoch"].map(_paired_epoch)
    return df


def _paired_epoch(epoch: str) -> str:
    if str(epoch).startswith("Baseline"):
        return "Baseline"
    if str(epoch).startswith("Injection"):
        return "Injection"
    return str(epoch)


def _condition(r: pd.Series) -> str:
    assay = r["assay"]
    epoch = str(r["epoch"])
    if assay == "melittin_only" and epoch == "Baseline":
        return "Baseline\n(melittin day)"
    if assay == "melittin_only" and epoch == "Injection":
        return "Melittin"
    if assay == "melittin_only" and epoch == "Injection2":
        return "Melittin\ninject 2"
    if assay == "melittin_only" and epoch == "Injection3":
        return "Melittin\ninject 3"
    if assay == "melittin_morphine" and epoch.startswith("Baseline"):
        return "Baseline\n(morphine day)"
    if assay == "melittin_morphine" and epoch == "Injection":
        return "Melittin\n+ morphine"
    return f"{assay}\n{epoch}"


def _style() -> None:
    plt.rcParams.update({
        "font.size": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.dpi": 160,
        "savefig.bbox": "tight",
    })


COND_ORDER = [
    "Baseline\n(melittin day)",
    "Melittin",
    "Baseline\n(morphine day)",
    "Melittin\n+ morphine",
]
COND_COLOR = {
    "Baseline\n(melittin day)": "#4C78A8",
    "Melittin": "#F58518",
    "Melittin\ninject 2": "#E45756",
    "Melittin\ninject 3": "#B279A2",
    "Baseline\n(morphine day)": "#72B7B2",
    "Melittin\n+ morphine": "#54A24B",
}


def mouse_means(df: pd.DataFrame, conditions: list[str]) -> pd.DataFrame:
    sub = df[df["condition"].isin(conditions) & df["mgs"].notna() & (df["n_fau"] >= 3)].copy()
    g = (sub.groupby(["mouse_id", "condition", "assay"], as_index=False)
         .agg(mean_mgs=("mgs", "mean"),
              n_frames=("mgs", "size"),
              mean_n_fau=("n_fau", "mean"),
              **{f"au_{a}": (a, "mean") for a in AU}))
    return g


def fig_algorithm(models) -> None:
    cropper, scorer, _q = models
    raw_path = QC / "20260218_1014m01b_19mWpCM2_base" / "raw_frame.jpg"
    if not raw_path.is_file():
        raw_path = QC / "1014m01b_16mWpCM2_base_20260408_f200.jpg"
    raw = cv2.cvtColor(cv2.imread(str(raw_path)), cv2.COLOR_BGR2RGB)
    crop, conf = cropper.crop(cv2.cvtColor(raw, cv2.COLOR_RGB2BGR))
    crop_rgb = cv2.cvtColor(crop, cv2.COLOR_BGR2RGB)
    net = cv2.resize(crop_rgb, (224, 224), interpolation=cv2.INTER_LINEAR)
    pred = scorer.score(cv2.cvtColor(crop_rgb, cv2.COLOR_RGB2BGR))
    labels, vals = [], []
    for name, v in zip(AU, pred):
        if int(v) == 3:
            labels.append(AU_LABEL[name] + "\nNR")
            vals.append(0)
        else:
            labels.append(AU_LABEL[name])
            vals.append(int(v))
    mgs = int(np.sum([v for v in pred if int(v) != 3]))
    n_fau = int(np.sum([int(v) != 3 for v in pred]))

    fig = plt.figure(figsize=(12.2, 7.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.25, 0.85], width_ratios=[1.45, 1, 1, 1.15],
                          hspace=0.45, wspace=0.35)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])
    ax3 = fig.add_subplot(gs[0, 3])
    for ax, im, title in (
        (ax0, raw, "1. Raw IR frame"),
        (ax1, crop_rgb, f"2. Face crop\n(det {conf:.2f})"),
        (ax2, net, "3. Network input\n224 px"),
    ):
        ax.imshow(im, cmap="gray")
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
    bar_vals = [v if "NR" not in lb else 0 for lb, v in zip(labels, vals)]
    colors = ["#4C78A8" if "NR" not in lb else "#DDDDDD" for lb in labels]
    ax3.bar(range(5), bar_vals, color=colors)
    for i, lb in enumerate(labels):
        if "NR" in lb:
            ax3.text(i, 0.08, "NR", ha="center", va="bottom", fontsize=8, color="#666666")
    ax3.set_xticks(range(5), [AU_LABEL[a] for a in AU], rotation=30, ha="right")
    ax3.set_ylim(0, 2.4)
    ax3.set_ylabel("AU score (0, 1, or 2)")
    ax3.set_title(f"4. Five AU heads\nMGS = {mgs}  ({n_fau}/5 rateable)")
    ax3.axhline(0, color="k", lw=0.6)

    axn = fig.add_subplot(gs[1, :])
    axn.axis("off")
    axn.set_xlim(0, 10)
    axn.set_ylim(0, 2)
    steps = [
        (0.15, "Raw face\nIR, ~1024×768"),
        (2.15, "Square crop\nYOLO, else center"),
        (4.15, "ViT-B/16\n(x−0.5175)/0.14"),
        (6.15, "FC 2048 → 1024\n5 × 4-class heads"),
        (8.15, "Argmax\n0 / 1 / 2 / NR"),
    ]
    for x, text in steps:
        axn.add_patch(plt.Rectangle((x, 0.45), 1.7, 1.15, fill=False, lw=1.2, ec="#333333"))
        axn.text(x + 0.85, 1.02, text, ha="center", va="center", fontsize=10)
    for x in (1.9, 3.9, 5.9, 7.9):
        axn.annotate("", xy=(x + 0.2, 1.02), xytext=(x, 1.02),
                     arrowprops=dict(arrowstyle="->", color="#333333", lw=1.2))
    axn.set_title(
        "GrimACE (Sturman, Lab Animal 2026) applied offline  ·  weights rwitz/GrimACE  ·  "
        "MGS = sum of rateable action units (0–10). NR = not rateable.",
        loc="left", fontsize=11, pad=8,
    )
    fig.suptitle(
        "How GrimACE is applied to a crystal-skull face frame\n"
        "Example: 20260218  1014m01b baseline (melittin day, CS1014-1b)",
        fontsize=13, y=1.02,
    )
    fig.savefig(FIG / "fig01_algorithm_steps.png")
    plt.close(fig)
    print("wrote", FIG / "fig01_algorithm_steps.png")


def _qc_image(session: str, frame: int) -> Path | None:
    folder = None
    for p in QC.iterdir():
        if p.is_dir() and p.name.endswith(session) or p.name == f"{session}" or session in p.name:
            folder = p
            break
    # exact folder is {date}_{session}; caller passes that
    return None


def _saved_crop(row: pd.Series) -> Path | None:
    folder = QC / f"{int(row['behavior_date'])}_{row['behavior_session']}"
    hits = sorted(folder.glob(f"frame_*_f{int(row['frame'])}.jpg"))
    return hits[0] if hits else None


def fig_representatives(df: pd.DataFrame) -> None:
    """One scored crop per mouse × main condition, only frames we saved."""
    want = COND_ORDER
    mice = ["1014m01b", "1014m02w", "CS0204-1b", "CS0204-2w", "mouse1f", "mouse2c"]
    mice = [m for m in mice if m in set(df["mouse_id"])]
    fig, axes = plt.subplots(
        len(mice), len(want), figsize=(2.5 * len(want), 2.35 * len(mice)),
    )
    if len(mice) == 1:
        axes = np.array([axes])
    for j, cond in enumerate(want):
        axes[0, j].set_title(cond.replace("\n", " "), fontsize=10, pad=6)
    for i, mouse in enumerate(mice):
        for j, cond in enumerate(want):
            ax = axes[i, j]
            sub = df[(df["mouse_id"] == mouse) & (df["condition"] == cond)].copy()
            im = None
            title = "no session"
            if len(sub):
                sub["_path"] = sub.apply(_saved_crop, axis=1)
                saved = sub[sub["_path"].notna()].sort_values(
                    ["n_fau", "frame_quality"], ascending=False)
                if len(saved):
                    row = saved.iloc[0]
                    im = cv2.cvtColor(cv2.imread(str(row["_path"])), cv2.COLOR_BGR2RGB)
                    mgs = row["mgs"]
                    mgs_s = f"{mgs:.0f}" if pd.notna(mgs) else "—"
                    title = f"MGS {mgs_s} ({int(row['n_fau'])}/5)"
                else:
                    title = "crop not saved"
            ax.set_xticks([])
            ax.set_yticks([])
            if im is None:
                ax.set_facecolor("#F4F4F4")
                ax.text(0.5, 0.5, title, ha="center", va="center", fontsize=8, transform=ax.transAxes)
            else:
                ax.imshow(im)
                ax.set_xlabel(title, fontsize=8)
            if j == 0:
                ax.set_ylabel(mouse, fontsize=10)
    fig.suptitle(
        "Representative scored crops (highest number of rateable AUs among saved frames)\n"
        "Gray cells: that mouse has no recording in that condition",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(FIG / "fig02_representative_faces.png")
    plt.close(fig)
    print("wrote", FIG / "fig02_representative_faces.png")


def fig_mouse_mgs(df: pd.DataFrame) -> None:
    g = mouse_means(df, COND_ORDER)
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    mice = sorted(g["mouse_id"].unique())
    markers = ["o", "s", "D", "^", "v", "P"]
    for mi, mouse in enumerate(mice):
        sub = g[g["mouse_id"] == mouse]
        xs, ys = [], []
        for xi, cond in enumerate(COND_ORDER):
            hit = sub[sub["condition"] == cond]
            if len(hit):
                xs.append(xi)
                ys.append(float(hit["mean_mgs"].iloc[0]))
        if xs:
            ax.plot(xs, ys, "-", color="#888888", lw=0.8, zorder=1)
            ax.scatter(xs, ys, s=55, marker=markers[mi % len(markers)],
                       label=mouse, zorder=2)
    ax.set_xticks(range(len(COND_ORDER)), [c.replace("\n", " ") for c in COND_ORDER])
    ax.set_ylabel("Mean MGS (sum of rateable AUs)")
    ax.set_ylim(0, 4)
    ax.legend(title="Mouse", frameon=False, ncol=2, loc="upper right")
    ax.set_title(
        "Mouse-level GrimACE: baseline vs melittin vs melittin + morphine\n"
        "Axis stops at 4; the full MGS scale is 0–10. Every mouse mean in this sample is below 3."
    )
    fig.text(
        0.01, 0.01,
        "Point = one mouse (mean of frames with ≥3 rateable AUs). Frames are not independent. "
        "1014m01b/02w on the melittin day are CS1014; morphine day is 20260408. "
        "Whisker is often not rateable, which keeps the sum low.",
        fontsize=8, color="#333333", transform=fig.transFigure,
    )
    fig.subplots_adjust(bottom=0.18)
    fig.savefig(FIG / "fig03_mouse_mgs_conditions.png")
    plt.close(fig)
    print("wrote", FIG / "fig03_mouse_mgs_conditions.png")


def fig_au(df: pd.DataFrame) -> None:
    sub = df[df["condition"].isin(COND_ORDER)].copy()
    # Mouse mean of each AU, then mean ± range across mice.
    rows = []
    for cond in COND_ORDER:
        for mouse, sm in sub[sub["condition"] == cond].groupby("mouse_id"):
            rec = {"condition": cond, "mouse_id": mouse}
            for a in AU:
                rec[a] = sm[a].mean(skipna=True)
                rec[a + "_nr"] = sm[a].isna().mean()
            rows.append(rec)
    m = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))
    x = np.arange(len(AU))
    width = 0.18
    for i, cond in enumerate(COND_ORDER):
        hit = m[m["condition"] == cond]
        means = [hit[a].mean(skipna=True) for a in AU]
        axes[0].bar(x + (i - 1.5) * width, means, width,
                    color=COND_COLOR[cond], label=cond.replace("\n", " "))
    axes[0].set_xticks(x, [AU_LABEL[a] for a in AU])
    axes[0].set_ylabel("Mean AU score (0–2), mouse means")
    axes[0].set_ylim(0, 2.2)
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].set_title("Action-unit scores by condition")

    for i, cond in enumerate(COND_ORDER):
        hit = m[m["condition"] == cond]
        nr = [100 * hit[a + "_nr"].mean() for a in AU]
        axes[1].bar(x + (i - 1.5) * width, nr, width,
                    color=COND_COLOR[cond], label=cond.replace("\n", " "))
    axes[1].set_xticks(x, [AU_LABEL[a] for a in AU])
    axes[1].set_ylabel("Not-rateable frames (%)")
    axes[1].set_ylim(0, 100)
    axes[1].set_title("Which AUs the network refuses")
    fig.suptitle("GrimACE action units  ·  average of mouse means", fontsize=13)
    fig.tight_layout()
    fig.savefig(FIG / "fig04_action_units.png")
    plt.close(fig)
    print("wrote", FIG / "fig04_action_units.png")


def fig_paired_delta(df: pd.DataFrame) -> None:
    """INJ − BL within assay, one point per mouse. Injection epoch only (not inject 2/3)."""
    sub = df[df["epoch"].isin(["Baseline", "Injection", "Baseline2"]) & (df["n_fau"] >= 3) & df["mgs"].notna()]
    rows = []
    for (mouse, assay), sm in sub.groupby(["mouse_id", "assay"]):
        bl = sm[sm["epoch"].isin(["Baseline", "Baseline2"])]["mgs"]
        inj = sm[sm["epoch"] == "Injection"]["mgs"]
        if len(bl) and len(inj):
            rows.append({
                "mouse_id": mouse,
                "assay": "Melittin" if assay == "melittin_only" else "Melittin + morphine",
                "delta": float(inj.mean() - bl.mean()),
                "bl": float(bl.mean()),
                "inj": float(inj.mean()),
            })
    d = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    assays = ["Melittin", "Melittin + morphine"]
    for i, assay in enumerate(assays):
        hit = d[d["assay"] == assay]
        ax.scatter(np.full(len(hit), i) + np.linspace(-0.08, 0.08, len(hit) or 1)[:len(hit)],
                   hit["delta"], s=60, zorder=2)
        hit = hit.sort_values("delta")
        for k, (_, r) in enumerate(hit.iterrows()):
            ax.annotate(
                r["mouse_id"],
                (i, r["delta"]),
                textcoords="offset points",
                xytext=(12, (k - (len(hit) - 1) / 2) * 11),
                fontsize=8,
                va="center",
                arrowprops=dict(arrowstyle="-", color="#BBBBBB", lw=0.6),
            )
        if len(hit):
            ax.plot([i - 0.2, i + 0.2], [hit["delta"].mean()] * 2, color="k", lw=2, zorder=1)
    ax.axhline(0, color="#888888", lw=0.8)
    ax.set_xlim(-0.45, 1.55)
    ax.set_xticks([0, 1], assays)
    ax.set_ylim(-2, 2)
    ax.set_ylabel("Δ MGS (injection − baseline)")
    ax.set_title(
        "Within-mouse change after injection\n"
        "Full MGS scale is 0–10. Black bar = mean of mice. n = mice, not frames."
    )
    fig.savefig(FIG / "fig05_delta_inj_minus_bl.png")
    plt.close(fig)
    print("wrote", FIG / "fig05_delta_inj_minus_bl.png")
    d.to_csv(FIG / "mouse_delta_mgs.csv", index=False)


def fig_later_injections(df: pd.DataFrame) -> None:
    mouse = "1014m02w"
    order = [
        "Baseline\n(melittin day)",
        "Melittin",
        "Melittin\ninject 2",
        "Melittin\ninject 3",
        "Baseline\n(morphine day)",
        "Melittin\n+ morphine",
    ]
    sub = df[(df["mouse_id"] == mouse) & df["condition"].isin(order) & (df["n_fau"] >= 3)]
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    short = [
        "BL\nmelittin day", "Melittin", "Inject 2", "Inject 3",
        "BL\nmorphine day", "Melittin\n+ morphine",
    ]
    data, labels = [], []
    for cond, lab in zip(order, short):
        vals = sub.loc[sub["condition"] == cond, "mgs"].dropna().to_numpy()
        data.append(vals)
        labels.append(lab)
        if len(vals):
            ax.scatter(np.full(len(vals), len(data) - 1), vals, s=28, color=COND_COLOR.get(cond, "#333"),
                       alpha=0.85, zorder=2)
            ax.plot([len(data) - 1.15, len(data) - 0.85], [vals.mean(), vals.mean()], color="k", lw=2)
    ax.set_xticks(range(len(labels)), labels)
    ax.set_ylabel("Frame MGS (rateable-AU sum)")
    ax.set_ylim(0, 5)
    ax.set_title(
        "1014m02w across injections\n"
        "Full MGS scale is 0–10. Points are sampled frames; black ticks are frame means."
    )
    fig.savefig(FIG / "fig06_1014m02w_timepoints.png")
    plt.close(fig)
    print("wrote", FIG / "fig06_1014m02w_timepoints.png")


def fig_yield(df: pd.DataFrame) -> None:
    g = (df.groupby(["behavior_date", "mouse_id", "epoch", "condition"], as_index=False)
         .agg(n=("mgs", "size"),
              frac_face=("face_detected", "mean"),
              frac5=("n_fau", lambda s: np.mean(s >= 5)),
              frac3=("n_fau", lambda s: np.mean(s >= 3))))
    g.to_csv(FIG / "session_yield.csv", index=False)
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    labels = (g["behavior_date"].astype(str) + " " + g["mouse_id"].astype(str) + " " + g["epoch"].astype(str)).tolist()
    x = np.arange(len(g))
    ax.bar(x - 0.15, 100 * g["frac3"], 0.3, label="≥3 AUs rateable", color="#4C78A8")
    ax.bar(x + 0.15, 100 * g["frac5"], 0.3, label="5/5 AUs (MaxFAU)", color="#F58518")
    ax.set_xticks(x, labels, rotation=70, ha="right", fontsize=8)
    ax.set_ylabel("Sampled frames (%)")
    ax.set_ylim(0, 100)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    ax.set_title("How often a full GrimACE score is possible")
    fig.tight_layout()
    fig.savefig(FIG / "fig07_rateable_yield.png")
    plt.close(fig)
    print("wrote", FIG / "fig07_rateable_yield.png")


def main() -> None:
    _style()
    FIG.mkdir(parents=True, exist_ok=True)
    df = load_frames()
    df.to_csv(OUT / "grimace_all_frames.csv", index=False)
    print(f"{len(df)} frames, {df['behavior_session'].nunique()} sessions")
    from score_grimace_videos import (  # noqa: WPS433
        FaceCropper, FrameQualityModel, GrimaceScoreModel, MODEL_DIR, _device,
    )
    device = _device()
    cropper = FaceCropper(MODEL_DIR / "face_cropper.pt", device)
    scorer = GrimaceScoreModel(device)
    scorer.load(MODEL_DIR / "grimace_scorer.pt", device)
    qmodel = FrameQualityModel(MODEL_DIR / "frame_quality.pt", device)
    fig_algorithm((cropper, scorer, qmodel))
    fig_representatives(df)
    fig_mouse_mgs(df)
    fig_au(df)
    fig_paired_delta(df)
    fig_later_injections(df)
    fig_yield(df)


if __name__ == "__main__":
    main()
