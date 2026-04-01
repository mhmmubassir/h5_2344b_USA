#!/usr/bin/env python3
"""
3_prepare_wetlab_benchmark_data.py

Prepare a clean benchmarking folder by:
1) generating experimental_ddg.csv from SPR KD values
2) collecting benchmark-ready protocol CSVs into one place

By default it collects:
- Rosetta final HA-corrected outputs:
    *_vs_HA_final_ddg_avg.csv
    *_vs_HA_final_ddg_min.csv
- Flex ddG summary CSVs:
    *flexddg*.csv

This avoids mixing the benchmarking folder with intermediate HA-only or raw
per-glyco Rosetta tables that are not the final comparison products.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

RT_KCAL = 0.001987 * 298.15

KD_DATA = {
    "WT": {"kd3": 1.38e-7, "kd6": 1e-3},
    "Q226L": {"kd3": 1e-3, "kd6": 2.31e-5},
    "G228S": {"kd3": 1e-3, "kd6": 1e-3},
    "N224K/Q226L": {"kd3": 1e-3, "kd6": 1.18e-6},
    "N224K/Q226L/G228S": {"kd3": 1e-6, "kd6": 5.54e-6},
    "E190D": {"kd3": 1e-3, "kd6": 1e-3},
}


def kd_to_ddg(kd: float, kd_wt: float) -> float:
    return float(RT_KCAL * np.log(kd / kd_wt))


def write_experimental_csv(outfile: Path) -> None:
    wt3 = KD_DATA["WT"]["kd3"]
    wt6 = KD_DATA["WT"]["kd6"]
    rows = []
    for mutation, values in KD_DATA.items():
        rows.append({
            "mutation": mutation,
            "exp_ddG_23": kd_to_ddg(values["kd3"], wt3),
            "exp_ddG_26": kd_to_ddg(values["kd6"], wt6),
        })
    pd.DataFrame(rows).to_csv(outfile, index=False)


def collect_csvs(search_root: Path, out_dir: Path) -> int:
    patterns = [
        "**/*_vs_HA_final_ddg_avg.csv",
        "**/*_vs_HA_final_ddg_min.csv",
        "**/*flexddg*.csv",
    ]
    copied = 0
    for pattern in patterns:
        for src in search_root.glob(pattern):
            if not src.is_file():
                continue
            if out_dir in src.parents:
                continue
            dst = out_dir / src.name
            shutil.copy2(src, dst)
            copied += 1
    return copied


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare wet-lab benchmarking input folder.")
    parser.add_argument("--search-root", default=".", help="Root folder to search for Rosetta/Flex ddG result CSVs.")
    parser.add_argument("--outdir", default="statistical_correlation", help="Output folder for benchmarking input files.")
    args = parser.parse_args()

    search_root = Path(args.search_root).resolve()
    out_dir = Path(args.outdir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    exp_csv = out_dir / "experimental_ddg.csv"
    write_experimental_csv(exp_csv)
    copied = collect_csvs(search_root, out_dir)

    print(f"✅ Wrote experimental file: {exp_csv}")
    print(f"✅ Copied {copied} protocol CSV files into: {out_dir}")


if __name__ == "__main__":
    main()
