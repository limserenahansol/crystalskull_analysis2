import pandas as pd
from pathlib import Path

QC = Path(r"c:\Users\hsollim\behavior_task\crystal_skull_analysis\grimace\grimace_output\qc_crops")
df = pd.read_csv(r"c:\Users\hsollim\behavior_task\crystal_skull_analysis\grimace\grimace_output\grimace_all_frames.csv")


def path_of(r):
    folder = QC / f"{int(r.behavior_date)}_{r.behavior_session}"
    hits = list(folder.glob(f"frame_*_f{int(r.frame)}.jpg"))
    return str(hits[0]) if hits else ""


df["img"] = df.apply(path_of, axis=1)
have = df[df.img != ""].copy()
print("with image", len(have), "/", len(df))


def show(mask, n=3, label=""):
    sub = have[mask].sort_values("frame_quality", ascending=False).head(n)
    print("==", label, "n", int(mask.sum()), "saved", len(have[mask]))
    for _, r in sub.iterrows():
        print(
            f"  {r.mouse_id} {r.behavior_date} {r.epoch} mgs={r.mgs} "
            f"orb={r.orbital} nose={r.nose} cheek={r.cheek} ear={r.ear} wh={r.whisker} "
            f"q={r.frame_quality:.2f}"
        )
        print("   ", r.img)


show(have.orbital == 0, 2, "orbital 0")
show(have.orbital == 1, 2, "orbital 1")
show(have.nose == 1, 5, "nose 1")
show(have.cheek == 1, 5, "cheek 1")
show(have.whisker == 1, 5, "whisker 1")
show(have.ear == 0, 2, "ear 0")
show((have.mgs == 1) & (have.ear == 1) & (have.orbital == 0), 2, "MGS 1 ear only")
show(have.mgs == 4, 2, "MGS 4")
