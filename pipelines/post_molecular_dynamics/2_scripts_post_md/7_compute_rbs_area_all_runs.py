#!/usr/bin/env python3
"""
compute_rbs_area_all_runs.py

Walk through subdirectories, find *_rbs folders produced by
rbs_geometry_all_directory.sh, and for each protomer (1–3):

  - read:
      d_pX_133_187.dat
      d_pX_187_223.dat
      d_pX_133_223.dat
  - compute triangle area per frame via Heron's formula
  - write: rbs_area_pX.dat (frame, area)
  - collect summary stats into rbs_area_summary.csv

Run this from the parent directory that contains your MD run folders.
"""

import os
from pathlib import Path
import numpy as np
import csv

# Names of distance files per protomer
DIST_FILES = {
    1: ("d_p1_133_187.dat", "d_p1_187_223.dat", "d_p1_133_223.dat"),
    2: ("d_p2_133_187.dat", "d_p2_187_223.dat", "d_p2_133_223.dat"),
    3: ("d_p3_133_187.dat", "d_p3_187_223.dat", "d_p3_133_223.dat"),
}

def load_distances(path):
    """
    Load cpptraj distance file with columns:
        frame_index   distance(Å)
    Returns (frames, distances).
    """
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    frames = data[:, 0].astype(int)
    d = data[:, 1]
    return frames, d

def heron_area(a, b, c):
    """
    Heron's formula: area given side lengths a, b, c.
    Supports numpy arrays.
    """
    s = 0.5 * (a + b + c)
    area_sq = s * (s - a) * (s - b) * (s - c)
    area_sq = np.where(area_sq < 0.0, 0.0, area_sq)  # numerical safety
    return np.sqrt(area_sq)

def process_rbs_dir(rbs_dir: Path):
    """
    Process a single *_rbs directory.
    Returns list of summary dicts (one per protomer).
    """
    summaries = []
    run_name = rbs_dir.parent.name  # MD run directory name

    print(f"[INFO] Processing {rbs_dir}")

    for prot in (1, 2, 3):
        f_a, f_b, f_c = DIST_FILES[prot]
        paths = [rbs_dir / f_a, rbs_dir / f_b, rbs_dir / f_c]
        if not all(p.exists() for p in paths):
            print(f"[WARN] Missing distance files for protomer {prot} in {rbs_dir}, skipping this protomer.")
            continue

        frames_a, d_a = load_distances(paths[0])
        frames_b, d_b = load_distances(paths[1])
        frames_c, d_c = load_distances(paths[2])

        # sanity check frames
        if not (np.array_equal(frames_a, frames_b) and np.array_equal(frames_a, frames_c)):
            raise ValueError(f"Frame indices do not match in {rbs_dir} for protomer {prot}")

        frames = frames_a
        area = heron_area(d_a, d_b, d_c)

        # write per-frame area file
        out_dat = rbs_dir / f"rbs_area_p{prot}.dat"
        np.savetxt(out_dat, np.column_stack([frames, area]),
                   fmt=["%8d", "%12.6f"])
        print(f"  [INFO] Wrote {out_dat.name} (protomer {prot})")

        # summary stats
        summary = {
            "run": run_name,
            "rbs_dir": rbs_dir.name,
            "protomer": prot,
            "nframes": int(area.size),
            "mean_area": float(area.mean()),
            "std_area": float(area.std()),
            "min_area": float(area.min()),
            "max_area": float(area.max()),
        }
        summaries.append(summary)

    return summaries

def main():
    root = Path(".").resolve()
    all_summaries = []

    # Walk the tree and find *_rbs dirs
    for dirpath, dirnames, filenames in os.walk(root):
        d = Path(dirpath)
        if d.name.endswith("_rbs"):
            summaries = process_rbs_dir(d)
            all_summaries.extend(summaries)

    if not all_summaries:
        print("[WARN] No *_rbs directories or no valid data found.")
        return

    # Write global summary CSV
    out_csv = root / "rbs_area_summary.csv"
    fieldnames = ["run", "rbs_dir", "protomer", "nframes",
                  "mean_area", "std_area", "min_area", "max_area"]

    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in all_summaries:
            writer.writerow(row)

    print(f"[INFO] Wrote summary to {out_csv}")

if __name__ == "__main__":
    main()
