"""Run the three GrimACE steps on saved raw frames and draw what each one returned."""
from __future__ import annotations

import sys
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from score_grimace_videos import (  # noqa: E402
    FaceCropper,
    FrameQualityModel,
    GrimaceScoreModel,
    MODEL_DIR,
    _device,
    _square_center,
)

QC = HERE / "grimace_output" / "qc_crops"
FIG = HERE / "grimace_output" / "figures"
AU = ["orbital", "nose", "cheek", "ear", "whisker"]
AU_KO = ["눈 둘레", "코", "볼", "귀", "수염"]

def _grab_mp4(video: Path, frame_i: int) -> np.ndarray | None:
    cap = cv2.VideoCapture(str(video))
    if not cap.isOpened():
        return None
    cap.set(cv2.CAP_PROP_POS_FRAMES, frame_i)
    ok, fr = cap.read()
    cap.release()
    return fr if ok else None


SHARE = Path(r"\\stg-tnr50a1.stanford.edu\crystal_skull\20260408\behavior")
MP4_FRAMES = [
    ("1014m01b f0 YOLO miss", SHARE / "1014m01b_16mWpCM2_base_20260408.mp4", 0),
    ("1014m01b f9450 YOLO hit", SHARE / "1014m01b_16mWpCM2_base_20260408.mp4", 9450),
    ("1014m01b f6825 YOLO hit", SHARE / "1014m01b_16mWpCM2_base_20260408.mp4", 6825),
]

RAWS = [
    ("1014m01b baseline", QC / "1014m01b_16mWpCM2_base_20260408_f200.jpg"),
    ("1014m02w inject", QC / "1014m02w_16mWpCM2_inject_20260408_f200.jpg"),
    ("1014m01b mel BL", QC / "20260218_1014m01b_19mWpCM2_base" / "raw_frame.jpg"),
    ("CS0204-2w inject", QC / "20260219_fcs_02042w_inject1_20260219" / "raw_frame.jpg"),
    ("mouse2c BL", QC / "20260219_mouse2c_base_20260219" / "raw_frame.jpg"),
    ("1014m02w mel inject3", QC / "20260218_1014m02w_19mWpCM2_inject3" / "raw_frame.jpg"),
]


def yolo_box(cropper: FaceCropper, frame: np.ndarray):
    gray = frame.mean(axis=2).astype(frame.dtype)
    gray = np.stack([gray] * 3, axis=2)
    res = cropper.model(gray, verbose=False, device=cropper.device)[0]
    if res.boxes is None or len(res.boxes) == 0:
        return None, 0.0
    conf = res.boxes.conf.cpu().numpy()
    i = int(np.argmax(conf))
    box = res.boxes.xyxy[i].cpu().numpy().astype(int)
    return box, float(conf[i])


def au_text(pred: np.ndarray) -> tuple[str, int, int]:
    bits = []
    mgs = 0
    n = 0
    for name, v in zip(AU, pred):
        if int(v) == 3:
            bits.append(f"{name} NR")
        else:
            bits.append(f"{name} {int(v)}")
            mgs += int(v)
            n += 1
    return "   ".join(bits), mgs, n


def main() -> None:
    device = _device()
    cropper = FaceCropper(MODEL_DIR / "face_cropper.pt", device, conf=0.10)
    scorer = GrimaceScoreModel(device)
    scorer.load(MODEL_DIR / "grimace_scorer.pt", device)
    qmodel = FrameQualityModel(MODEL_DIR / "frame_quality.pt", device)

    rows = []
    sources = [(label, path, None) for label, path in RAWS]
    for label, video, frame_i in MP4_FRAMES:
        sources.append((label, video, frame_i))
    for label, path, frame_i in sources:
        if frame_i is None:
            if not path.is_file():
                print("missing", path)
                continue
            bgr = cv2.imread(str(path))
        else:
            print("reading", path.name, "frame", frame_i, flush=True)
            bgr = _grab_mp4(path, frame_i)
        if bgr is None:
            print("unreadable", label)
            continue
        q = qmodel.score(bgr)
        box, conf = yolo_box(cropper, bgr)
        crop, used_conf = cropper.crop(bgr)
        pred = scorer.score(crop)
        rows.append({
            "label": label, "bgr": bgr, "q": q, "box": box, "conf": conf,
            "crop": crop, "pred": pred, "used_conf": used_conf,
        })
        text, mgs, n = au_text(pred)
        print(f"{label:22} q={q:.3f} yolo={conf:.3f} box={box} MGS={mgs} ({n}/5) {text}")

    FIG.mkdir(parents=True, exist_ok=True)
    from matplotlib import font_manager
    font_manager.fontManager.addfont(r"C:\Windows\Fonts\calibri.ttf")
    plt.rcParams.update({
        "font.size": 10,
        "font.family": "Calibri",
        "axes.unicode_minus": False,
        "figure.facecolor": "white",
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    })
    fig_pipeline(max(rows, key=lambda r: r["conf"]))
    fig_quality(sorted(rows, key=lambda r: r["q"])[:6])
    hits = sorted(rows, key=lambda r: r["conf"], reverse=True)
    fig_yolo((hits[:3] + hits[-3:]) if len(hits) > 6 else hits)
    fig_vit(rows)
    fig_verdict(rows)
    print("figures written")


def _rgb(bgr):
    return cv2.cvtColor(bgr, cv2.COLOR_BGR2RGB)


def fig_pipeline(r) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(16.5, 4.6))
    axes[0].imshow(_rgb(r["bgr"]))
    axes[0].set_title("1. Raw frame\nour IR face")
    vis = r["bgr"].copy()
    axes[1].imshow(_rgb(vis))
    axes[1].set_title(f"2. Quality network\nscore {r['q']:.2f}  (0-1)")
    boxed = r["bgr"].copy()
    if r["box"] is not None:
        x1, y1, x2, y2 = r["box"]
        cv2.rectangle(boxed, (x1, y1), (x2, y2), (0, 180, 255), 3)
        title = f"3. YOLO face box\nconfidence {r['conf']:.2f}"
    else:
        title = "3. YOLO\nno box"
    axes[2].imshow(_rgb(boxed))
    axes[2].set_title(title)
    axes[3].imshow(_rgb(r["crop"]))
    text, mgs, n = au_text(r["pred"])
    axes[3].set_title(f"4. ViT scores this crop\nMGS {mgs}  ({n}/5)")
    for a in axes:
        a.set_xticks([])
        a.set_yticks([])
    fig.suptitle(
        f"{r['label']}  |  Real model output. The box is the YOLO coordinate, not a hand drawing.\n{text}",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(FIG / "figS1_three_steps_one_frame.png")
    plt.close(fig)


def fig_quality(rows) -> None:
    order = sorted(rows, key=lambda r: r["q"])
    fig, axes = plt.subplots(1, len(order), figsize=(3.1 * len(order), 4.2))
    for ax, r in zip(axes, order):
        ax.imshow(_rgb(r["bgr"]))
        ax.set_title(f"{r['label']}\nquality {r['q']:.2f}", fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(
        "Quality-network scores on our frames (scale 0-1). Sharp front faces still score low.",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(FIG / "figS2_quality_network.png")
    plt.close(fig)


def fig_yolo(rows) -> None:
    n = len(rows)
    cols = 3
    rws = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rws, cols, figsize=(12.5, 4.0 * rws))
    axes = np.atleast_1d(axes).ravel()
    for ax, r in zip(axes, rows):
        boxed = r["bgr"].copy()
        if r["box"] is not None and r["conf"] >= 0.10:
            x1, y1, x2, y2 = r["box"]
            cv2.rectangle(boxed, (x1, y1), (x2, y2), (0, 180, 255), 3)
            status = "box used"
            color = "#0D5C63"
        elif r["box"] is not None:
            x1, y1, x2, y2 = r["box"]
            cv2.rectangle(boxed, (x1, y1), (x2, y2), (40, 40, 220), 3)
            status = "low confidence, center crop"
            color = "#9B1B30"
        else:
            status = "no box, center crop"
            color = "#9B1B30"
        ax.imshow(_rgb(boxed))
        ax.set_title(f"{r['label']}\nYOLO {r['conf']:.2f}  {status}", fontsize=10, color=color)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle(
        "YOLO boxes on our frames. Orange boxes passed the 0.10 cutoff and were used, but they sit on the paws or the top of the head.\n"
        "Clear front-face frames have no box.",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(FIG / "figS3_yolo_boxes.png")
    plt.close(fig)


def fig_vit(rows) -> None:
    # two clearest contrast: highest yolo conf and lowest, plus show scores
    show = sorted(rows, key=lambda r: r["conf"], reverse=True)[:4]
    fig, axes = plt.subplots(2, 4, figsize=(14.5, 7.2))
    for i, r in enumerate(show):
        axes[0, i].imshow(_rgb(r["crop"]))
        text, mgs, n = au_text(r["pred"])
        axes[0, i].set_title(f"{r['label']}\ncrop sent to the ViT", fontsize=11)
        axes[0, i].set_xticks([])
        axes[0, i].set_yticks([])
        vals = []
        colors = []
        for v in r["pred"]:
            if int(v) == 3:
                vals.append(0)
                colors.append("#DDDDDD")
            else:
                vals.append(int(v))
                colors.append("#4C78A8")
        axes[1, i].bar(range(5), vals, color=colors)
        axes[1, i].set_ylim(0, 2.4)
        axes[1, i].set_xticks(range(5), ["orbital", "nose", "cheek", "ear", "whisker"], rotation=30, ha="right", fontsize=10)
        axes[1, i].set_title(f"MGS {mgs} ({n}/5)", fontsize=10)
        for k, v in enumerate(r["pred"]):
            if int(v) == 3:
                axes[1, i].text(k, 0.08, "NR", ha="center", fontsize=9, color="#666666")
    fig.suptitle(
        "The vision transformer does not draw a box on each feature. It scores the whole crop. NR = not rateable.",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(FIG / "figS4_vit_scores.png")
    plt.close(fig)


def fig_verdict(rows) -> None:
    by_q = min(rows, key=lambda r: r["q"])
    by_box = max(rows, key=lambda r: r["conf"])
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.2))
    axes[0].imshow(_rgb(by_q["bgr"]))
    axes[0].set_title(f"Quality network: fails\nsharp face, score {by_q['q']:.2f} / 1", color="#9B1B30", fontsize=13)
    boxed = by_box["bgr"].copy()
    if by_box["box"] is not None:
        x1, y1, x2, y2 = by_box["box"]
        cv2.rectangle(boxed, (x1, y1), (x2, y2), (0, 180, 255), 3)
    axes[1].imshow(_rgb(boxed))
    axes[1].set_title(f"YOLO face crop: fails\nbest box ({by_box['conf']:.2f}) is on the paws", color="#9B1B30", fontsize=13)
    axes[2].imshow(_rgb(by_box["crop"]))
    text, mgs, n = au_text(by_box["pred"])
    axes[2].set_title(f"ViT score: not a face check\nMGS {mgs} ({n}/5) on the bad crop", color="#9B1B30", fontsize=13)
    axes[2].set_xlabel(text, fontsize=8)
    for a in axes:
        a.set_xticks([])
        a.set_yticks([])
    fig.suptitle("Check of the three steps on our frames", fontsize=16)
    fig.tight_layout()
    fig.savefig(FIG / "figS5_which_step_fails.png")
    plt.close(fig)


if __name__ == "__main__":
    main()
