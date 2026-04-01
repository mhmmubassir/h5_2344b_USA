#!/usr/bin/env python3
"""
1_group_ddg_outputs.py

Organize collected Rosetta/Flex ddG output folders into a clean benchmarking layout.

Input:
    collected_ddg_outputs/

Output:
    organized_ddg_groups/
        flexddg_all/
        9dip_<score>_<bb>/
        9dio_<score>_<bb>_ori/
        9dio_<score>_<bb>_rvmut/
        others/

Notes
-----
- 9dip folders are grouped by score function and backbone setting.
- 9dio HA-only folders are copied into both ori and rvmut groups so the final
  comparison step can run cleanly for each branch.
- Unmatched directories are copied into organized_ddg_groups/others/ for audit.
"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

STD_RX = re.compile(
    r"""^
        (?P<system>9di[op])_
        (?P<glyco>23|26|HAonly)_
        (?P<score>bn16|cr15)
        (?:_(?P<replicate>[A-Za-z0-9]+))?
        -ddg-5ddg_
        (?P<score2>bn16|cr15)_
        (?P<bb>bb[123])
        (?:-[^/]+)?
        -ddg_outputs$
    """,
    re.VERBOSE,
)

FLEX_RX = re.compile(r".+_flexddg-ddg_outputs$")


def copy_tree_if_missing(src: Path, dest_dir: Path) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / src.name
    if dest.exists():
        return dest
    shutil.copytree(src, dest)
    return dest


def main() -> None:
    parser = argparse.ArgumentParser(description="Group raw ddG output folders for benchmarking.")
    parser.add_argument("--src", default="collected_ddg_outputs", help="Source directory containing raw ddG folders.")
    parser.add_argument("--dst", default="organized_ddg_groups", help="Destination directory for grouped folders.")
    args = parser.parse_args()

    src_root = Path(args.src).resolve()
    dst_root = Path(args.dst).resolve()
    dst_root.mkdir(parents=True, exist_ok=True)

    if not src_root.is_dir():
        raise SystemExit(f"❌ Source folder not found: {src_root}")

    n_total = 0
    n_grouped = 0
    n_other = 0

    for item in sorted(src_root.iterdir()):
        if not item.is_dir():
            continue
        n_total += 1
        name = item.name

        if FLEX_RX.search(name):
            copy_tree_if_missing(item, dst_root / "flexddg_all")
            n_grouped += 1
            continue

        m = STD_RX.match(name)
        if not m or m.group("score") != m.group("score2"):
            copy_tree_if_missing(item, dst_root / "others")
            n_other += 1
            continue

        g = m.groupdict()
        system = g["system"]
        glyco = g["glyco"]
        score = g["score"]
        replicate = g["replicate"]
        bb = g["bb"]

        if system == "9dip":
            if replicate is None and glyco in {"23", "26", "HAonly"}:
                copy_tree_if_missing(item, dst_root / f"9dip_{score}_{bb}")
                n_grouped += 1
            else:
                copy_tree_if_missing(item, dst_root / "others")
                n_other += 1
            continue

        if system == "9dio":
            if glyco == "26" and replicate in {"ori", "rvmut"}:
                copy_tree_if_missing(item, dst_root / f"9dio_{score}_{bb}_{replicate}")
                n_grouped += 1
            elif glyco == "HAonly" and replicate is None:
                for rep in ("ori", "rvmut"):
                    copy_tree_if_missing(item, dst_root / f"9dio_{score}_{bb}_{rep}")
                n_grouped += 1
            else:
                copy_tree_if_missing(item, dst_root / "others")
                n_other += 1
            continue

        copy_tree_if_missing(item, dst_root / "others")
        n_other += 1

    print(f"✅ Grouping finished: {n_grouped} grouped, {n_other} sent to 'others', from {n_total} source folders.")
    print(f"   Output: {dst_root}")


if __name__ == "__main__":
    main()
