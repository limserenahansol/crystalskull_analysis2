"""Show what an AU score of 1 means on our own face frames."""
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np

FIG = Path(r"c:\Users\hsollim\behavior_task\crystal_skull_analysis\grimace\grimace_output\figures")
QC = Path(r"c:\Users\hsollim\behavior_task\crystal_skull_analysis\grimace\grimace_output\qc_crops")

panels = [
    (
        QC / "20260219_fcs_02042w_inject1_20260219" / "frame_003_f1372.jpg",
        "MGS = 0",
        "CS0204-2w  melittin injection",
        "orbital 0   nose 0   cheek 0\near 0   whisker 0",
        "Every unit the network scored is absent.",
    ),
    (
        QC / "20260219_fcs_02042w_base_20260219" / "frame_004_f20766.jpg",
        "MGS = 1",
        "CS0204-2w  baseline",
        "orbital 0   nose 0   cheek 0\near 1   whisker 0",
        "The 1 point is the ear only.\nEyes, nose, cheek, whiskers scored 0.",
    ),
    (
        QC / "20260408_1014m01b_16mWpCM2_base_20260408" / "frame_001_f1050.jpg",
        "MGS = 3",
        "1014m01b  morphine-day baseline",
        "orbital 1   nose 0   cheek 0\near 1   whisker 1",
        "Three separate 1-point calls.\n1 + 1 + 1 = 3. No unit scored 2.",
    ),
    (
        QC / "20260218_1014m01b_19mWpCM2_inject1" / "frame_001_f2039.jpg",
        "MGS = 4  (bad crop)",
        "1014m01b  melittin injection",
        "orbital 1   nose 1   cheek 1\near 1   whisker NR",
        "The crop missed the eyes.\nNose/cheek = 1 here is not a full face.",
    ),
]


def main() -> None:
    plt.rcParams.update({
        "font.size": 10,
        "figure.facecolor": "white",
        "savefig.dpi": 160,
        "savefig.bbox": "tight",
    })
    fig, axes = plt.subplots(1, 4, figsize=(14.2, 5.6))
    for ax, (path, title, who, scores, note) in zip(axes, panels):
        im = cv2.cvtColor(cv2.imread(str(path)), cv2.COLOR_BGR2RGB)
        ax.imshow(im)
        ax.set_title(f"{title}\n{who}", fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(scores + "\n" + note, fontsize=8)
    fig.suptitle(
        "What the numbers mean on our crystal-skull face frames\n"
        "Each action unit is 0 (absent), 1 (moderate), or 2 (obvious). "
        "MGS = sum of the units that were rateable. In these 192 frames, no unit was ever scored 2.",
        fontsize=12,
    )
    fig.tight_layout()
    out = FIG / "fig08_what_one_point_means.png"
    fig.savefig(out)
    plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
