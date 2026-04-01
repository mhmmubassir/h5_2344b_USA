#!/usr/bin/env python3
"""
2_summarize_rosetta_ddg.py

Summarize Rosetta .ddg files inside organized groups, build avg/min heatmaps,
and compute final bound-vs-HA-only ΔΔG tables for benchmarking.

Expected grouped layout:
    organized_ddg_groups/
        9dip_<score>_<bb>/
            ..._23_.../
            ..._26_.../
            ..._HAonly_.../
        9dio_<score>_<bb>_ori/
            ..._26_ori_.../
            ..._HAonly_.../
        9dio_<score>_<bb>_rvmut/
            ..._26_rvmut_.../
            ..._HAonly_.../

Outputs are written inside each group root and copied from each glyco subfolder.
A central folder of heatmaps is also created.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import shutil
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm

WT_PATTERNS = [
    re.compile(r"COMPLEX:\s+Round\d+:\s+WT_[^:]*:\s+([-]?\d+\.\d+)"),
    re.compile(r"SCORE:.*WT_.*\s([-]?\d+\.\d+)$"),
]
MUT_PATTERNS = [
    re.compile(r"COMPLEX:\s+Round\d+:\s+MUT_[^:]*:\s+([-]?\d+\.\d+)"),
    re.compile(r"SCORE:.*MUT_.*\s([-]?\d+\.\d+)$"),
]


def nice_cmap():
    return sns.diverging_palette(240, 10, s=95, l=35, as_cmap=True)


def plot_heatmap(df: pd.DataFrame, value_col: str, title: str, outfile: Path, cbar_label: str) -> None:
    if df.empty:
        return
    max_abs = max(abs(float(df[value_col].min())), abs(float(df[value_col].max())))
    max_abs = max(max_abs, 1e-6)
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)

    plt.figure(figsize=(max(8, len(df) * 0.4), 2.5))
    ax = sns.heatmap(
        df.set_index("mutation")[[value_col]].T,
        annot=True,
        fmt=".2f",
        annot_kws={"size": 6},
        cmap=nice_cmap(),
        norm=norm,
        cbar_kws={"label": cbar_label, "extend": "both"},
        linewidths=0.2,
    )
    ax.set_title(title, weight="bold", fontsize=10)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=6)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.collections[0].colorbar.ax.tick_params(labelsize=6)
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()


def extract_scores(score_file: Path) -> Optional[Dict[str, float]]:
    wt: List[float] = []
    mut: List[float] = []
    with score_file.open() as handle:
        for line in handle:
            for pat in WT_PATTERNS:
                m = pat.search(line)
                if m:
                    wt.append(float(m.group(1)))
                    break
            for pat in MUT_PATTERNS:
                m = pat.search(line)
                if m:
                    mut.append(float(m.group(1)))
                    break
    if not wt or not mut:
        return None
    return {
        "avg": sum(mut) / len(mut) - sum(wt) / len(wt),
        "min": min(mut) - min(wt),
    }


def summarize_variant_dir(ddg_dir: Path, prefix: str, group_root: Path) -> None:
    records = {"avg": [], "min": []}
    for root, _, files in os.walk(ddg_dir, followlinks=True):
        for fn in files:
            if not fn.lower().endswith(".ddg"):
                continue
            path = Path(root) / fn
            scores = extract_scores(path)
            if scores is None:
                continue
            mutation = path.stem
            records["avg"].append({"mutation": mutation, "ddG": scores["avg"]})
            records["min"].append({"mutation": mutation, "ddG": scores["min"]})

    for kind in ("avg", "min"):
        data = records[kind]
        if not data:
            continue
        df = pd.DataFrame(data).sort_values("mutation").reset_index(drop=True)
        csv_path = ddg_dir / f"{prefix}_ddg_{kind}_summary.csv"
        txt_path = ddg_dir / f"{prefix}_ddg_{kind}_summary.txt"
        png_path = ddg_dir / f"{prefix}_heatmap_{kind}.png"
        df.to_csv(csv_path, index=False)
        txt_path.write_text(df.to_string(index=False) + "\n")
        plot_heatmap(df, "ddG", f"{prefix} ΔΔG ({kind})", png_path, "ΔΔG")
        copy_outputs_up(ddg_dir, group_root, prefix)


def copy_outputs_up(src_dir: Path, dst_dir: Path, prefix: str) -> None:
    for kind in ("avg", "min"):
        for ext in ("csv", "txt", "png"):
            if ext == "png":
                fn = f"{prefix}_heatmap_{kind}.png"
            else:
                fn = f"{prefix}_ddg_{kind}_summary.{ext}"
            src = src_dir / fn
            if src.exists():
                shutil.copy2(src, dst_dir / fn)


def final_ddg(bound_dir: Path, unbound_dir: Path, bound_prefix: str, unbound_prefix: str, out_prefix: str, out_dir: Path) -> None:
    for kind in ("avg", "min"):
        bound_csv = bound_dir / f"{bound_prefix}_ddg_{kind}_summary.csv"
        unbound_csv = unbound_dir / f"{unbound_prefix}_ddg_{kind}_summary.csv"
        if not (bound_csv.exists() and unbound_csv.exists()):
            continue

        bound = pd.read_csv(bound_csv).set_index("mutation")["ddG"]
        unbound = pd.read_csv(unbound_csv).set_index("mutation")["ddG"]
        merged = (bound - unbound).rename("final_ddG").to_frame().reset_index()

        csv_path = out_dir / f"{out_prefix}_final_ddg_{kind}.csv"
        txt_path = out_dir / f"{out_prefix}_final_ddg_{kind}.txt"
        png_path = out_dir / f"{out_prefix}_final_heatmap_{kind}.png"
        merged.to_csv(csv_path, index=False)
        txt_path.write_text(merged.to_string(index=False) + "\n")
        plot_heatmap(merged, "final_ddG", f"{out_prefix} Final ΔΔG ({kind})", png_path, "Final ΔΔG")


def collect_heatmaps(group_root: Path, central_dir: Path) -> None:
    central_dir.mkdir(parents=True, exist_ok=True)
    for src in group_root.rglob("*_heatmap_*.png"):
        try:
            shutil.copy2(src, central_dir / src.name)
        except shutil.SameFileError:
            pass
    for src in group_root.rglob("*_final_heatmap_*.png"):
        try:
            shutil.copy2(src, central_dir / src.name)
        except shutil.SameFileError:
            pass


def map_group_subdirs(group_path: Path) -> Dict[str, Path]:
    mapping: Dict[str, Path] = {}
    for sub in sorted(group_path.iterdir()):
        if not sub.is_dir():
            continue
        name = sub.name
        if "_23_" in name:
            mapping["23"] = sub
        elif "_26_" in name:
            mapping["26"] = sub
        elif "_HAonly_" in name:
            mapping["HAonly"] = sub
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize Rosetta .ddg files and create benchmark-ready outputs.")
    parser.add_argument("--groups", default="organized_ddg_groups", help="Grouped ddG folder from step 1.")
    parser.add_argument("--heatmaps", default="bound_vs_unbound_heatmaps", help="Central folder for copied heatmaps.")
    args = parser.parse_args()

    groups_root = Path(args.groups).resolve()
    heatmap_root = Path(args.heatmaps).resolve()

    if not groups_root.is_dir():
        raise SystemExit(f"❌ Group folder not found: {groups_root}")

    group_dirs = [d for d in sorted(groups_root.iterdir()) if d.is_dir() and d.name.startswith(("9dip_", "9dio_"))]
    if not group_dirs:
        raise SystemExit("❌ No grouped 9dip_/9dio_ folders found.")

    for group in group_dirs:
        mapping = map_group_subdirs(group)

        for glyco, subdir in mapping.items():
            prefix = f"{group.name}_{glyco}"
            summarize_variant_dir(subdir, prefix, group)

        if "HAonly" in mapping:
            if group.name.startswith("9dip_"):
                for glyco in ("23", "26"):
                    if glyco in mapping:
                        final_ddg(
                            mapping[glyco],
                            mapping["HAonly"],
                            f"{group.name}_{glyco}",
                            f"{group.name}_HAonly",
                            f"{group.name}_{glyco}_vs_HA",
                            group,
                        )
            elif group.name.startswith("9dio_") and "26" in mapping:
                final_ddg(
                    mapping["26"],
                    mapping["HAonly"],
                    f"{group.name}_26",
                    f"{group.name}_HAonly",
                    f"{group.name}_26_vs_HA",
                    group,
                )

    collect_heatmaps(groups_root, heatmap_root)
    print(f"✅ Rosetta summaries finished. Groups: {groups_root}")
    print(f"✅ Heatmaps copied to: {heatmap_root}")


if __name__ == "__main__":
    main()
