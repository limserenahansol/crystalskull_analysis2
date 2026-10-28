"""Score crystal-skull face videos with official GrimACE weights (rwitz/GrimACE).

Samples frames in place (does not copy zip/mp4). Writes per-session CSVs that
grimace_load.m can ingest, plus QC crops.

These weights were trained on GrimACE-box front IR video. Mesoscope side/face
cameras are out of distribution — inspect QC crops before trusting MGS numbers.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MODEL_DIR = HERE / "models"
QC_DIR = HERE / "grimace_output" / "qc_crops"
SCORE_DIR = HERE / "grimace_output"
MATCH_CSV = HERE / "behavior_calcium_match.csv"
AU_NAMES = ("orbital", "nose", "cheek", "ear", "whisker")


def _device() -> str:
    try:
        import torch
        return "cuda" if torch.cuda.is_available() else "cpu"
    except ImportError:
        return "cpu"


def download_weights() -> None:
    MODEL_DIR.mkdir(parents=True, exist_ok=True)
    needed = ["frame_quality.pt", "face_cropper.pt", "grimace_scorer.pt"]
    missing = [n for n in needed if not (MODEL_DIR / n).is_file()]
    if not missing:
        return
    from huggingface_hub import hf_hub_download
    for name in missing:
        print(f"downloading {name} ...", flush=True)
        hf_hub_download(repo_id="rwitz/GrimACE", filename=name, local_dir=str(MODEL_DIR))


class GrimaceScoreModel:
    def __init__(self, device: str = "cpu"):
        import torch
        from torchvision.models import vit_b_16

        self.device = device
        self.base_model = vit_b_16(weights=None)
        self.base_model.heads = torch.nn.Identity()
        n_feat = self.base_model(torch.zeros(1, 3, 224, 224)).shape[1]
        self.hidden_layers = torch.nn.Sequential(
            torch.nn.Sequential(torch.nn.Dropout(0.0), torch.nn.Linear(n_feat, 2048), torch.nn.ReLU()),
            torch.nn.Sequential(torch.nn.Dropout(0.0), torch.nn.Linear(2048, 1024), torch.nn.ReLU()),
        )
        self.au_heads = torch.nn.ModuleList([
            torch.nn.Sequential(torch.nn.Linear(1024, 4)) for _ in range(5)
        ])
        self._mod = torch.nn.Module()
        self._mod.base_model = self.base_model
        self._mod.hidden_layers = self.hidden_layers
        self._mod.au_heads = self.au_heads

    def load(self, path: Path, device: str) -> None:
        import torch
        raw = torch.load(path, map_location=device, weights_only=False)
        if isinstance(raw, dict) and "state_dict" in raw:
            raw = raw["state_dict"]
        if not isinstance(raw, dict):
            raise TypeError(f"unexpected checkpoint type {type(raw)}")
        # Strip common prefixes from the official app dump.
        cleaned = {}
        for k, v in raw.items():
            k2 = k
            for pref in ("module.", "grimace_scorer_model.", "model."):
                if k2.startswith(pref):
                    k2 = k2[len(pref):]
            cleaned[k2] = v
        missing, unexpected = self._mod.load_state_dict(cleaned, strict=False)
        if missing:
            print("grimace_scorer missing keys:", missing[:8], flush=True)
        if unexpected:
            print("grimace_scorer unexpected keys:", unexpected[:8], flush=True)
        self._mod.to(device)
        self._mod.eval()
        self.device = device

    def preprocess(self, frame: np.ndarray):
        import cv2
        import torch
        frame = cv2.resize(frame, (224, 224), interpolation=cv2.INTER_LINEAR)
        if frame.ndim == 2 or frame.shape[2] == 1:
            frame = np.stack([frame] * 3, axis=2)
        x = torch.from_numpy(frame.astype(np.float32) / 255.0).permute(2, 0, 1)
        return (x - 0.5175) / 0.14

    def score(self, crop: np.ndarray) -> np.ndarray:
        import torch
        x = self.preprocess(crop).unsqueeze(0).to(self.device)
        with torch.no_grad():
            out = self._mod.base_model(x)
            out = self._mod.hidden_layers(out)
            out = torch.cat([h(out) for h in self._mod.au_heads], dim=1)
            out = out.softmax(dim=1).reshape(out.shape[0], 5, 4)
        return out.argmax(dim=2).squeeze(0).cpu().numpy()


class FrameQualityModel:
    def __init__(self, path: Path, device: str):
        import timm
        import torch

        self.device = device
        self.base = timm.create_model("mobilenetv3_large_100", in_chans=1, pretrained=False)
        n_feat = self.base(torch.zeros(1, 1, 224, 224)).shape[1]
        self.head = torch.nn.Sequential(torch.nn.Linear(n_feat, 1), torch.nn.Sigmoid())
        self.mod = torch.nn.Module()
        self.mod.base_model = self.base
        self.mod.score_head = torch.nn.Sequential(torch.nn.Identity(), self.head[0], self.head[1])
        raw = torch.load(path, map_location=device, weights_only=False)
        if isinstance(raw, dict) and "state_dict" in raw:
            raw = raw["state_dict"]
        if isinstance(raw, dict):
            cleaned = {k.replace("module.", ""): v for k, v in raw.items()}
            miss, unexp = self.mod.load_state_dict(cleaned, strict=False)
            if miss:
                print("frame_quality missing keys:", miss[:8], flush=True)
        self.mod.to(device).eval()

    def score(self, frame: np.ndarray) -> float:
        import cv2
        import torch
        h, w = frame.shape[:2]
        if h != w:
            d = abs(h - w)
            dl, dr = d // 2, d - d // 2
            frame = frame[dl:h - dr, :] if h > w else frame[:, dl:w - dr]
        frame = cv2.resize(frame, (224, 224), interpolation=cv2.INTER_LINEAR)
        if frame.ndim == 3:
            frame = frame.mean(axis=2, keepdims=True)
        elif frame.ndim == 2:
            frame = frame[:, :, None]
        x = torch.from_numpy(frame.astype(np.float32) / 255.0).permute(2, 0, 1)
        x = ((x - 0.6901) / 0.1750).unsqueeze(0).to(self.device)
        with torch.no_grad():
            return float(self.mod.score_head(self.mod.base_model(x)).squeeze().cpu())


class FaceCropper:
    def __init__(self, path: Path, device: str, conf: float = 0.10):
        from ultralytics import YOLO
        self.model = YOLO(str(path))
        self.device = device
        self.conf = conf

    def crop(self, frame: np.ndarray) -> tuple[np.ndarray | None, float]:
        # Mesoscope face cam is already a tight front view. If YOLO (trained on
        # a wide GrimACE-box frame) finds nothing, score the square center crop.
        gray = frame
        if frame.ndim == 3:
            gray = frame.mean(axis=2).astype(frame.dtype)
            gray = np.stack([gray] * 3, axis=2)
        res = self.model(gray, verbose=False, device=self.device)[0]
        if res.boxes is None or len(res.boxes) == 0:
            return _square_center(frame), 0.0
        conf = res.boxes.conf.cpu().numpy()
        i = int(np.argmax(conf))
        if conf[i] < self.conf:
            return _square_center(frame), float(conf[i])
        x1, y1, x2, y2 = res.boxes.xyxy[i].cpu().numpy().astype(int)
        side = max(x2 - x1, y2 - y1)
        cx, cy = (x1 + x2) // 2, (y1 + y2) // 2
        x1, y1 = cx - side // 2, cy - side // 2
        x2, y2 = x1 + side, y1 + side
        h, w = frame.shape[:2]
        pad = np.zeros((side, side, frame.shape[2] if frame.ndim == 3 else 1), dtype=frame.dtype)
        xa1, ya1 = max(0, x1), max(0, y1)
        xa2, ya2 = min(w, x2), min(h, y2)
        ox1, oy1 = xa1 - x1, ya1 - y1
        crop_src = frame[ya1:ya2, xa1:xa2]
        if crop_src.size == 0:
            return _square_center(frame), float(conf[i])
        if frame.ndim == 2:
            pad = np.zeros((side, side), dtype=frame.dtype)
            pad[oy1:oy1 + (ya2 - ya1), ox1:ox1 + (xa2 - xa1)] = crop_src
        else:
            pad[oy1:oy1 + (ya2 - ya1), ox1:ox1 + (xa2 - xa1)] = crop_src
        return pad, float(conf[i])


def _square_center(frame: np.ndarray) -> np.ndarray:
    h, w = frame.shape[:2]
    side = min(h, w)
    y0 = (h - side) // 2
    x0 = (w - side) // 2
    return frame[y0:y0 + side, x0:x0 + side]


def _to_gray3(frame: np.ndarray) -> np.ndarray:
    import cv2
    if frame.ndim == 2:
        return cv2.cvtColor(frame, cv2.COLOR_GRAY2BGR)
    if frame.shape[2] == 4:
        return cv2.cvtColor(frame, cv2.COLOR_BGRA2BGR)
    return frame


def sample_video(path: Path, n: int, every_s: float) -> list[tuple[int, float, np.ndarray]]:
    import cv2
    cap = cv2.VideoCapture(str(path))
    if not cap.isOpened():
        return []
    fps = cap.get(cv2.CAP_PROP_FPS) or 30.0
    nframes = int(cap.get(cv2.CAP_PROP_FRAME_COUNT) or 0)
    step = max(1, int(round(fps * every_s)))
    idxs = list(range(0, max(nframes, 1), step))
    if len(idxs) > n:
        sel = np.linspace(0, len(idxs) - 1, n).astype(int)
        idxs = [idxs[i] for i in sel]
    out = []
    for i in idxs:
        cap.set(cv2.CAP_PROP_POS_FRAMES, i)
        ok, fr = cap.read()
        if ok:
            out.append((i, i / fps, _to_gray3(fr)))
    cap.release()
    return out


def _first_jpeg_names(folder: Path, n: int, timeout_s: float = 25.0) -> list[str]:
    """Stop after n names. A full dir of the share JPEG dump hangs."""
    import subprocess
    # cmd cannot use a UNC path as cwd, so pass the full wildcard.
    proc = subprocess.Popen(
        ["cmd", "/c", "dir", "/b", "/a-d", str(folder / "*.jpg")],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
        bufsize=1,
    )
    names: list[str] = []
    try:
        assert proc.stdout is not None
        while len(names) < n:
            line = proc.stdout.readline()
            if not line:
                break
            name = line.strip()
            if name.lower().endswith(".jpg"):
                names.append(name)
    finally:
        proc.kill()
        proc.wait(timeout=timeout_s)
    return names


def sample_stills(folder: Path, n: int) -> list[tuple[int, float, np.ndarray]]:
    """Read only the first n JPEGs. Do not walk the rest of the dump."""
    import cv2
    names = _first_jpeg_names(folder, n)
    out = []
    for i, name in enumerate(names):
        fr = cv2.imread(str(folder / name))
        if fr is None:
            continue
        m = re.search(r"(\d+)", name)
        frame_i = int(m.group(1)) if m else i
        out.append((frame_i, float(i), _to_gray3(fr)))
    return out


def load_match_rows() -> list[dict]:
    if not MATCH_CSV.is_file():
        raise FileNotFoundError(f"run build_behavior_manifest.py first: {MATCH_CSV}")
    with MATCH_CSV.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def score_session(row: dict, models, n_frames: int, every_s: float) -> list[dict]:
    import cv2
    sess = row["behavior_session"]
    folder = Path(row["behavior_path"])
    videos = [v for v in (row.get("videos") or "").split(";") if v]
    frames: list[tuple[int, float, np.ndarray]] = []
    source = ""
    if videos:
        source = videos[0]
        frames = sample_video(folder / videos[0], n_frames, every_s)
    elif int(row.get("n_stills") or 0) > 0 or True:
        source = "stills"
        frames = sample_stills(folder, n_frames)
    qc_dir = QC_DIR / f"{row['behavior_date']}_{sess}"
    qc_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    cropper, scorer, qmodel = models
    for k, (fi, t, fr) in enumerate(frames):
        crop, det_c = cropper.crop(fr)
        q = qmodel.score(fr)
        au = [np.nan] * 5
        n_fau = 0
        mgs = np.nan
        used = crop if crop is not None else fr
        if crop is not None:
            pred = scorer.score(crop)
            au = []
            for v in pred:
                if int(v) == 3:
                    au.append(np.nan)
                else:
                    au.append(int(v))
            n_fau = int(np.sum(np.isfinite(au)))
            mgs = float(np.nansum(au)) if n_fau else np.nan
        if k < 8:
            vis = crop if crop is not None else fr
            cv2.imwrite(str(qc_dir / f"frame_{k:03d}_f{fi}.jpg"), vis)
            if k == 0:
                cv2.imwrite(str(qc_dir / "raw_frame.jpg"), fr)
        rec = {
            "behavior_session": sess,
            "behavior_date": row["behavior_date"],
            "assay": row["assay"],
            "mouse_id": row["mouse_id"],
            "epoch": row["epoch"],
            "calcium_mouse": row["calcium_mouse"],
            "source": source,
            "frame": fi,
            "time_s": t,
            "mgs": mgs,
            "n_fau": n_fau,
            "face_conf": det_c,
            "frame_quality": q,
            "face_detected": int(crop is not None),
        }
        for name, val in zip(AU_NAMES, au):
            rec[name] = val
        rows.append(rec)
    return rows


def write_session_csv(rows: list[dict], sess: str) -> Path:
    SCORE_DIR.mkdir(parents=True, exist_ok=True)
    path = SCORE_DIR / f"{sess}_grimace.csv"
    if not rows:
        return path
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return path


def summarize(all_rows: list[dict]) -> list[dict]:
    from collections import defaultdict
    groups = defaultdict(list)
    for r in all_rows:
        groups[(r["behavior_date"], r["mouse_id"], r["epoch"], r["behavior_session"])].append(r)
    out = []
    for key, rs in groups.items():
        det = [r for r in rs if r["face_detected"]]
        mgs = [r["mgs"] for r in det if r["n_fau"] >= 5 and np.isfinite(r["mgs"])]
        out.append({
            "behavior_date": key[0],
            "mouse_id": key[1],
            "epoch": key[2],
            "behavior_session": key[3],
            "n_sampled": len(rs),
            "n_face": len(det),
            "frac_face": len(det) / max(len(rs), 1),
            "n_max_fau": len(mgs),
            "mean_mgs": float(np.mean(mgs)) if mgs else np.nan,
            "mean_quality": float(np.mean([r["frame_quality"] for r in rs])),
        })
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-frames", type=int, default=24)
    p.add_argument("--every-s", type=float, default=10.0)
    p.add_argument("--max-sessions", type=int, default=0)
    p.add_argument("--date", default="", help="restrict to one behavior_date, e.g. 20260408")
    args = p.parse_args()

    download_weights()
    device = _device()
    print(f"device={device}", flush=True)
    cropper = FaceCropper(MODEL_DIR / "face_cropper.pt", device)
    scorer = GrimaceScoreModel(device)
    scorer.load(MODEL_DIR / "grimace_scorer.pt", device)
    qmodel = FrameQualityModel(MODEL_DIR / "frame_quality.pt", device)
    models = (cropper, scorer, qmodel)

    match_rows = load_match_rows()
    if args.date:
        match_rows = [r for r in match_rows if r["behavior_date"] == args.date]
    if args.max_sessions:
        match_rows = match_rows[: args.max_sessions]
    all_rows = []
    for row in match_rows:
        print(f"scoring {row['behavior_date']} {row['behavior_session']} ...", flush=True)
        try:
            recs = score_session(row, models, args.n_frames, args.every_s)
        except Exception as exc:
            print(f"  FAIL {exc}", flush=True)
            continue
        dest = write_session_csv(recs, f"{row['behavior_date']}_{row['behavior_session']}")
        print(f"  {len(recs)} frames -> {dest}", flush=True)
        all_rows.extend(recs)

    if not all_rows:
        print("no frames scored")
        return 1
    all_path = SCORE_DIR / "grimace_all_frames.csv"
    with all_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        w.writeheader()
        w.writerows(all_rows)
    summary = summarize(all_rows)
    sum_path = SCORE_DIR / "grimace_session_summary.csv"
    with sum_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)
    (SCORE_DIR / "grimace_run.json").write_text(json.dumps({
        "device": device,
        "n_frames_per_session": args.n_frames,
        "every_s": args.every_s,
        "weights": "huggingface:rwitz/GrimACE",
        "domain_shift": "GrimACE-box front IR vs mesoscope face camera",
        "n_sessions": len(match_rows),
        "n_scored_frames": len(all_rows),
    }, indent=2), encoding="utf-8")
    print(f"wrote {all_path}")
    print(f"wrote {sum_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
