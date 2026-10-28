"""Pull representative face frames from every behavior session and cache them locally.

Two media types on the share:
  mp4            seek with ffmpeg
  jpg_sequence   ~21,000 or 42,000 loose files named "Image_<idx> HH-MM-SS-mmm.jpg",
                 so the wall clock is in the filename and the index gives the ordering

Frames are sampled evenly across the middle 80% of each session (the first and last
minutes are when the experimenter is usually still at the rig) and written to
frames/<mouse>__<date>__<condition>/.

This does not score anything. It gives the raw material for the methods figure and for a
by-eye baseline vs melittin vs melittin+morphine comparison.
"""
from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "frames")
MANIFEST = os.path.join(HERE, "behavior_calcium_manifest.csv")
FFMPEG = os.path.expanduser(
    "~/Downloads/ffmpeg-master-latest-win64-gpl-shared/"
    "ffmpeg-master-latest-win64-gpl-shared/bin/ffmpeg.exe")
N_FRAMES = 8


def sample_positions(n_total: int, k: int):
    lo, hi = int(0.10 * n_total), int(0.90 * n_total)
    if hi <= lo:
        lo, hi = 0, max(n_total - 1, 0)
    step = max((hi - lo) // max(k, 1), 1)
    return [lo + i * step for i in range(k) if lo + i * step <= hi]


def from_mp4(path: str, fps: float, idxs, dst: str):
    got = []
    for i in idxs:
        t = i / max(fps, 1e-9)
        o = os.path.join(dst, f"f{i:06d}.png")
        if not os.path.exists(o):
            subprocess.run([FFMPEG, "-nostdin", "-v", "error", "-ss", f"{t:.3f}",
                            "-i", path, "-frames:v", "1", "-y", o],
                           capture_output=True, timeout=300)
        if os.path.exists(o):
            got.append(o)
    return got


def from_jpgs(folder: str, idxs, dst: str):
    """Match by the Image_<idx> token rather than sort order, which is lexicographic."""
    try:
        names = os.listdir(folder)
    except Exception:
        return []
    by_idx = {}
    pat = re.compile(r"Image_(\d+)")
    for n in names:
        m = pat.match(n)
        if m:
            by_idx[int(m.group(1))] = n
    got = []
    for i in idxs:
        n = by_idx.get(i)
        if n is None:
            continue
        o = os.path.join(dst, f"f{i:06d}.jpg")
        if not os.path.exists(o):
            try:
                shutil.copyfile(os.path.join(folder, n), o)
            except Exception:
                continue
        got.append(o)
    return got


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = list(csv.DictReader(open(MANIFEST)))
    print(f"{len(rows)} sessions in the manifest\n")
    for r in rows:
        tag = f"{r['mouse']}__{r['date']}__{r['condition']}{r['repeat']}"
        dst = os.path.join(OUT, tag)
        os.makedirs(dst, exist_ok=True)
        have = [f for f in os.listdir(dst) if f.lower().endswith((".png", ".jpg"))]
        if len(have) >= N_FRAMES:
            print(f"  {tag:44s} cached ({len(have)})")
            continue
        n = int(r["n_frames"]) if r["n_frames"] else 0
        if n == 0:
            print(f"  {tag:44s} no frame count, skipped")
            continue
        idxs = sample_positions(n, N_FRAMES)
        if r["media"] == "mp4":
            fps = float(r["fps"]) if r["fps"] else 35.0
            got = from_mp4(r["behavior_path"], fps, idxs, dst)
        else:
            got = from_jpgs(r["behavior_path"], idxs, dst)
        print(f"  {tag:44s} {len(got)}/{len(idxs)} frames  [{r['media']}]")
    print(f"\nframes under {OUT}")


if __name__ == "__main__":
    main()
