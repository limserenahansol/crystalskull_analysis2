"""Index the facial behavior recordings and pair them with the calcium sessions.

Three acquisition dates hold face video for the melittin / morphine experiments:

  20260408   melittin + morphine   .mp4 + per-frame timestamp .csv (+ .zip of frames)
  20260218   melittin only         loose .jpg sequences, ~21,000 frames per session
  20260219   melittin only         loose .jpg sequences

Two things this writes:

  behavior_index.csv    one row per behavior session: media type, frame count, duration,
                        wall-clock start/end taken from the timestamp CSV when present
  behavior_calcium_manifest.csv
                        the same rows joined to the calcium session they belong to, so a
                        later analysis can go from a cell to the face video that was
                        running while it was recorded

Scanning the share is slow (~100 s for a 21,000-file folder), so results are cached in
behavior_index_cache.json and reused unless --rescan is passed.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import time

SRV = "//stg-tnr50a1.stanford.edu/crystal_skull"
CAL = f"{SRV}/omer_hazon"
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "behavior_index_cache.json")
DATES = ("20260408", "20260218", "20260219")

FFPROBE = os.path.expanduser(
    "~/Downloads/ffmpeg-master-latest-win64-gpl-shared/"
    "ffmpeg-master-latest-win64-gpl-shared/bin/ffprobe.exe")

# Behavior folder / file prefixes -> the calcium mouse folder they belong to.
# 1014m01b and 1014m02w appear on BOTH 20260218 (melittin only) and 20260408
# (melittin + morphine). If that is right, the melittin vs melittin+morphine comparison
# is within-animal, not between-animal as the current analysis assumes. Confirm with Omer
# before relying on it - see the note printed at the end of this script.
MOUSE_MAP = {
    "1014m01b": ("1014m01b", "morphine_20260408/1014m01b"),
    "1014m02w": ("1014m02w", "morphine_20260408/1014m02w"),
    "fcs_02041b": ("CS0204-1b", "CS0204-1b"),
    "fcs_02042w": ("CS0204-2w", "CS0204-2w"),
    "mouse1f": ("mouse1f", None),      # FOV label, calcium folder not identified
    "mouse2c": ("mouse2c", None),
}
# which drug the session got, by date
DATE_ARM = {"20260408": "melittin_morphine", "20260218": "melittin_only",
            "20260219": "melittin_only"}


def parse_session(name: str):
    """Pull mouse, condition and power out of a behavior folder/file stem."""
    stem = name
    for ext in (".mp4", ".avi", ".mov"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
    low = stem.lower()
    mouse_key = next((k for k in MOUSE_MAP if k in low), None)
    if "inject" in low:
        cond = "Injection"
    elif "base" in low:
        cond = "Baseline"
    else:
        cond = "unknown"
    # base2 / inject1 / inject2 ... keep the repeat number so repeats stay distinct
    rep = ""
    for tag in ("base2", "base3", "inject1", "inject2", "inject3"):
        if tag in low:
            rep = tag[-1]
    pwr = ""
    for tok in stem.split("_"):
        if "mw" in tok.lower():
            pwr = tok
    return mouse_key, cond, rep, pwr


def probe(path: str):
    """Duration / frame count / size of a video, via ffprobe."""
    if not os.path.exists(FFPROBE):
        return {}
    try:
        out = subprocess.run(
            [FFPROBE, "-v", "error", "-select_streams", "v:0",
             "-show_entries", "stream=width,height,r_frame_rate,nb_frames",
             "-show_entries", "format=duration",
             "-of", "default=noprint_wrappers=1:nokey=1", path],
            capture_output=True, text=True, timeout=180).stdout.split()
        if len(out) >= 5:
            w, h, rate, nfr, dur = out[:5]
            num, den = (rate.split("/") + ["1"])[:2]
            return {"width": int(w), "height": int(h),
                    "fps": round(float(num) / float(den), 3),
                    "n_frames": int(nfr), "duration_s": round(float(dur), 2)}
    except Exception:
        pass
    return {}


def read_timestamps(path: str):
    """First/last wall-clock second and row count from a *_behavior_timestamp.csv."""
    try:
        with open(path, newline="") as fh:
            rows = fh.read().splitlines()
        if len(rows) < 3:
            return {}
        first = float(rows[1].split(",")[2])
        last = float(rows[-1].split(",")[2])
        n = len(rows) - 1
        return {"ts_rows": n, "t_start_s": first, "t_end_s": last,
                "wall_span_s": round(last - first, 2)}
    except Exception:
        return {}


def hhmmss(sec: float) -> str:
    sec = float(sec) % 86400
    return f"{int(sec//3600):02d}:{int(sec%3600//60):02d}:{int(sec%60):02d}"


def scan(rescan: bool = False):
    cache = {}
    if os.path.exists(CACHE) and not rescan:
        cache = json.load(open(CACHE))
    sessions = []
    for date in DATES:
        base = f"{SRV}/{date}/behavior"
        try:
            entries = sorted(os.listdir(base))
        except Exception as e:
            print(f"  {date}: cannot list ({e})")
            continue
        # --- flat .mp4 layout (20260408)
        for f in entries:
            if not f.lower().endswith(".mp4"):
                continue
            key = f"{date}/{f}"
            rec = cache.get(key)
            if rec is None:
                p = f"{base}/{f}"
                rec = {"media": "mp4", "path": p, "size_bytes": os.path.getsize(p)}
                rec.update(probe(p))
                ts = f"{base}/{os.path.splitext(f)[0]}_behavior_timestamp.csv"
                if os.path.exists(ts):
                    rec["timestamp_csv"] = ts
                    rec.update(read_timestamps(ts))
                cache[key] = rec
            sessions.append((date, os.path.splitext(f)[0], rec))
        # --- one-folder-per-session layout (20260218 / 20260219)
        for d in entries:
            p = f"{base}/{d}"
            if not os.path.isdir(p):
                continue
            key = f"{date}/{d}"
            rec = cache.get(key)
            if rec is None:
                t0 = time.time()
                try:
                    items = os.listdir(p)
                except Exception as e:
                    print(f"  {key}: {e}")
                    continue
                jpgs = [x for x in items if x.lower().endswith((".jpg", ".jpeg", ".png"))]
                vids = [x for x in items if x.lower().endswith((".mp4", ".avi", ".mov"))]
                csvs = [x for x in items if x.lower().endswith(".csv")]
                rec = {"media": "mp4" if vids else ("jpg_sequence" if jpgs else "unknown"),
                       "path": f"{p}/{vids[0]}" if vids else p,
                       "n_frames": len(vids and [] or jpgs) or None,
                       "n_files": len(items), "scan_s": round(time.time() - t0, 1)}
                if vids:
                    rec.update(probe(f"{p}/{vids[0]}"))
                if csvs:
                    rec["timestamp_csv"] = f"{p}/{csvs[0]}"
                    rec.update(read_timestamps(f"{p}/{csvs[0]}"))
                cache[key] = rec
                print(f"  scanned {key}: {rec['media']} "
                      f"{rec.get('n_frames') or rec.get('n_files')} frames "
                      f"({rec['scan_s']}s)")
            sessions.append((date, d, rec))
    json.dump(cache, open(CACHE, "w"), indent=1)
    return sessions


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rescan", action="store_true")
    args = ap.parse_args()

    print("scanning behavior share (cached unless --rescan)...")
    sessions = scan(args.rescan)
    print(f"\n{len(sessions)} behavior sessions\n")

    rows = []
    for date, name, rec in sessions:
        mouse_key, cond, rep, pwr = parse_session(name)
        mouse, cal_folder = MOUSE_MAP.get(mouse_key, (mouse_key or "?", None))
        cal_path = f"{CAL}/{cal_folder}/{cond}" if cal_folder and cond != "unknown" else ""
        rows.append({
            "date": date,
            "arm": DATE_ARM.get(date, ""),
            "mouse": mouse,
            "condition": cond,
            "repeat": rep,
            "power": pwr,
            "session": name,
            "media": rec.get("media", ""),
            "n_frames": rec.get("n_frames") or "",
            "fps": rec.get("fps", ""),
            "duration_s": rec.get("duration_s", ""),
            "wall_start": hhmmss(rec["t_start_s"]) if "t_start_s" in rec else "",
            "wall_end": hhmmss(rec["t_end_s"]) if "t_end_s" in rec else "",
            "wall_span_s": rec.get("wall_span_s", ""),
            "has_timestamps": "yes" if rec.get("timestamp_csv") else "no",
            "behavior_path": rec.get("path", ""),
            "timestamp_csv": rec.get("timestamp_csv", ""),
            "calcium_folder": cal_path,
            "calcium_exists": ("yes" if cal_path and os.path.isdir(cal_path)
                               else ("no" if cal_path else "")),
        })

    rows.sort(key=lambda r: (r["mouse"], r["date"], r["condition"], r["repeat"]))
    out = os.path.join(HERE, "behavior_calcium_manifest.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    hdr = f"{'mouse':<11}{'date':<10}{'arm':<18}{'cond':<11}{'media':<14}{'frames':>8} {'dur':>8}  {'wall':<18}{'cal'}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        wall = f"{r['wall_start']}-{r['wall_end']}" if r["wall_start"] else "-"
        print(f"{r['mouse']:<11}{r['date']:<10}{r['arm']:<18}"
              f"{r['condition']+r['repeat']:<11}{r['media']:<14}"
              f"{str(r['n_frames']):>8} {str(r['duration_s']):>8}  {wall:<18}"
              f"{r['calcium_exists']}")

    print(f"\nwrote {out}")

    mice = {}
    for r in rows:
        mice.setdefault(r["mouse"], set()).add(r["arm"])
    both = [m for m, a in mice.items() if len(a) > 1]
    if both:
        print("\n" + "!" * 76)
        print("These mice have face video under MORE THAN ONE drug arm:", ", ".join(both))
        print("If the calcium for both arms is also from these same animals, the")
        print("melittin vs melittin+morphine comparison is WITHIN-animal, not between-")
        print("animal as the current analysis treats it. That would replace the n = 4 vs 2")
        print("unpaired test with a paired one on the same mice. Confirm with Omer before")
        print("re-running - the folder names alone are not proof of animal identity.")
        print("!" * 76)


if __name__ == "__main__":
    main()
