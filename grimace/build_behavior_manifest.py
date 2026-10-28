"""Inventory face-behavior folders on the crystal-skull share and match them to EXTRACT.

Does not copy zip/mp4. Writes grimace/behavior_calcium_match.csv.
"""
from __future__ import annotations

import csv
import glob
import json
import os
import re
from pathlib import Path

SHARE = Path("//stg-tnr50a1.stanford.edu/crystal_skull")
HAZON = SHARE / "omer_hazon"
OUT_DIR = Path(__file__).resolve().parent
VIDEO_EXT = {".mp4", ".avi", ".mkv", ".mov"}
STILL_EXT = {".jpg", ".jpeg", ".png", ".tif", ".tiff"}

CALCIUM = {
    ("melittin_only", "CS1014-1b"): HAZON / "meltinonly" / "CS1014-1b",
    ("melittin_only", "CS1014-2w"): HAZON / "meltinonly" / "CS1014-2w",
    ("melittin_only", "CS0204-1b"): HAZON / "CS0204-1b",
    ("melittin_only", "CS0204-2w"): HAZON / "CS0204-2w",
    ("melittin_morphine", "1014m01b"): HAZON / "morphine_20260408" / "1014m01b",
    ("melittin_morphine", "1014m02w"): HAZON / "morphine_20260408" / "1014m02w",
}


def _norm_mouse(name: str) -> str:
    s = name.lower().replace(" ", "")
    if "1014m01b" in s or "cs1014-1b" in s or "cs1014_1b" in s:
        return "1014m01b"
    if "1014m02w" in s or "cs1014-2w" in s or "cs1014_2w" in s:
        return "1014m02w"
    if "cs0204-1b" in s or "02041b" in s or "fcs_02041b" in s:
        return "CS0204-1b"
    if "cs0204-2w" in s or "02042w" in s or "fcs_02042w" in s:
        return "CS0204-2w"
    if "mouse1f" in s or re.search(r"\bm1f\b", s):
        return "mouse1f"
    if "mouse2c" in s or re.search(r"\bm2c\b", s):
        return "mouse2c"
    return name


def _epoch(name: str) -> str:
    s = name.lower()
    if "inject3" in s:
        return "Injection3"
    if "inject2" in s:
        return "Injection2"
    if "inject" in s or re.search(r"\binj\b", s):
        return "Injection"
    if "base2" in s:
        return "Baseline2"
    if "base" in s:
        return "Baseline"
    return "unknown"


def _assay(date: str) -> str:
    if date == "20260408":
        return "melittin_morphine"
    if date in {"20260218", "20260219"}:
        return "melittin_only"
    return "unknown"


def _calcium_mouse(canon: str, assay: str) -> str:
    if assay == "melittin_only" and canon == "1014m01b":
        return "CS1014-1b"
    if assay == "melittin_only" and canon == "1014m02w":
        return "CS1014-2w"
    return canon


def _calcium_path(assay: str, ca_mouse: str) -> Path | None:
    return CALCIUM.get((assay, ca_mouse))


def _glob_names(folder: Path, pattern: str) -> list[str]:
    try:
        return sorted(Path(p).name for p in glob.glob(str(folder / pattern)))
    except OSError:
        return []


def _has_stills(folder: Path) -> tuple[int, str, str]:
    """Detect stills without walking tens of thousands of JPEGs."""
    try:
        first = next(glob.iglob(str(folder / "*.jpg")), None)
    except OSError:
        first = None
    if first is None:
        return 0, "", ""
    return -1, Path(first).name, ""  # -1 = present, not counted


def _row(date: str, sess_name: str, sess_path: Path,
         videos: list[str] | None = None, scan_files: bool = True) -> dict:
    assay = _assay(date)
    mouse = _norm_mouse(sess_name)
    epoch = _epoch(sess_name)
    ca_mouse = _calcium_mouse(mouse, assay)
    ca = _calcium_path(assay, ca_mouse)
    if videos is None:
        videos = _glob_names(sess_path, "*.mp4") + _glob_names(sess_path, "*.avi") if scan_files else []
    if scan_files:
        csvs = [n for n in _glob_names(sess_path, "*.csv") if sess_name in n]
        zips = [n for n in _glob_names(sess_path, "*.zip") if sess_name in n]
        n_jpg, jpg0, jpg1 = (0, "", "") if videos else (-1, "", "")
    else:
        csvs, zips = [], []
        n_jpg, jpg0, jpg1 = -1, "", ""
    match_ok = bool(ca and ca.is_dir())
    ca_epoch = "Injection" if epoch.startswith("Injection") else (
        "Baseline" if epoch.startswith("Baseline") else ""
    )
    ca_epoch_path = ""
    if ca and ca_epoch:
        p = ca / ca_epoch
        ca_epoch_path = str(p) if p.is_dir() else str(p)
    return {
        "behavior_date": date,
        "assay": assay,
        "behavior_session": sess_name,
        "behavior_path": str(sess_path),
        "mouse_id": mouse,
        "epoch": epoch,
        "n_stills": n_jpg,
        "still_first": jpg0,
        "still_last": jpg1,
        "videos": ";".join(videos),
        "n_videos": len(videos),
        "timestamp_csvs": ";".join(csvs),
        "zips": ";".join(zips),
        "calcium_mouse": ca_mouse,
        "calcium_path": str(ca) if ca else "",
        "calcium_epoch": ca_epoch,
        "calcium_epoch_path": ca_epoch_path,
        "calcium_exists": int(match_ok),
        "match_status": (
            "matched" if match_ok else
            ("unmatched_no_extract" if ca is None else "extract_missing")
        ),
    }


def inventory_date(date: str) -> list[dict]:
    root = SHARE / date / "behavior"
    if not root.is_dir():
        print(f"  MISSING {root}")
        return []
    try:
        dirs = sorted((Path(e.path) for e in os.scandir(root) if e.is_dir()), key=lambda p: p.name)
    except OSError as exc:
        print(f"  FAIL {root}: {exc}")
        return []
    rows = []
    if dirs:
        for sess in dirs:
            print(f"  dir {sess.name}", flush=True)
            rows.append(_row(date, sess.name, sess, scan_files=False))
    videos = _glob_names(root, "*.mp4") + _glob_names(root, "*.avi")
    if videos:
        groups: dict[str, list[str]] = {}
        for v in videos:
            key = Path(v).stem
            m = re.match(r"(.+?_(?:base2|base|inject\d*|inject|inj).*)", key, re.I)
            g = m.group(1) if m else key
            groups.setdefault(g, []).append(v)
        for g, vids in sorted(groups.items()):
            print(f"  video-group {g}", flush=True)
            rows.append(_row(date, g, root, videos=vids))
    if not dirs and not videos:
        rows.append(_row(date, root.name, root))
    return rows


def peek_timestamp_csv(path: Path, n: int = 3) -> list[str]:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as f:
            return [next(f).rstrip("\n") for _ in range(n + 1)]
    except (OSError, StopIteration) as exc:
        return [f"READ_FAIL {exc}"]


def main() -> int:
    dates = ["20260408", "20260218", "20260219"]
    rows = []
    for d in dates:
        print(f"scanning {d} ...", flush=True)
        part = inventory_date(d)
        print(f"  {len(part)} session(s)", flush=True)
        rows.extend(part)

    fields = list(rows[0].keys()) if rows else []
    inv_path = OUT_DIR / "behavior_inventory.csv"
    match_path = OUT_DIR / "behavior_calcium_match.csv"
    for dest in (inv_path, match_path):
        with dest.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {dest}")

    n_ok = sum(r["calcium_exists"] for r in rows)
    print(f"{n_ok}/{len(rows)} sessions have a matching EXTRACT folder")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
