"""Convert a face video to the MP4 PainFace expects (1080p-ish, 30 fps, yuv420p)."""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("video")
    p.add_argument("-o", "--out", default="")
    p.add_argument("--fps", type=int, default=30)
    args = p.parse_args()
    src = Path(args.video)
    if not src.is_file():
        sys.stderr.write(f"not found: {src}\n")
        return 2
    if shutil.which("ffmpeg") is None:
        sys.stderr.write("ffmpeg is not on PATH. Install it, then re-run.\n")
        return 2
    out = Path(args.out) if args.out else src.with_name(src.stem + "_painface.mp4")
    cmd = [
        "ffmpeg", "-y", "-i", str(src),
        "-r", str(args.fps),
        "-vf", "scale=1920:1080:force_original_aspect_ratio=decrease,pad=1920:1080:(ow-iw)/2:(oh-ih)/2",
        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-an",
        str(out),
    ]
    subprocess.check_call(cmd)
    print("wrote", out)
    print("Upload that file at http://painface.net then run grimace_load on the export CSV.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
