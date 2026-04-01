#!/usr/bin/env python3
"""
collect_rmsd.py

Recursively walk through a base directory, find all files named
'rmsd_protein.dat' and 'rmsd_complex.dat', and copy them into a folder 'rmsd'.

Each copied file will be renamed to include its relative path so you can
see where it came from.

Usage:
    python collect_rmsd.py /path/to/base_dir

If no base_dir is given, it uses the current directory.
"""

import os
import shutil
import argparse
from pathlib import Path


def sanitize_rel_path(rel_path: str) -> str:
    """
    Turn a relative path into a compact, safe label.

    - If rel_path is '.' -> 'ROOT'
    - Otherwise use only the last 3 directories of the path.
      e.g. 'big/long/path/run1/rep0' -> 'path__run1__rep0'
    """
    rel_path = rel_path.strip(os.sep)
    if rel_path in ("", "."):
        return "ROOT"

    parts = rel_path.split(os.sep)

    # Keep only the last 3 directories if longer
    if len(parts) > 3:
        parts = parts[-3:]

    safe_parts = []
    for p in parts:
        p = p.replace(" ", "_")
        p = p.replace(":", "_")
        p = p.replace("/", "_")
        p = p.replace("\\", "_")
        safe_parts.append(p)

    return "__".join(safe_parts)


def collect_rmsd(base_dir: Path, out_dir_name: str = "rmsd") -> None:
    base_dir = base_dir.resolve()
    out_dir = base_dir / out_dir_name
    out_dir.mkdir(exist_ok=True)

    print(f"Base directory : {base_dir}")
    print(f"Output directory: {out_dir}\n")

    count = 0
    target_files = ("rmsd_protein.dat", "rmsd_complex.dat")

    for root, dirs, files in os.walk(base_dir):
        # Prune directories:
        #  - don't recurse into the output folder itself
        #  - skip any folder whose name contains 'pdb' or 'png'
        dirs[:] = [
            d for d in dirs
            if d != out_dir_name and "pdb" not in d and "png" not in d
        ]

        present = [f for f in target_files if f in files]
        if not present:
            continue

        root_path = Path(root)
        rel_path = root_path.relative_to(base_dir)
        rel_str = sanitize_rel_path(str(rel_path))

        for fname in present:
            src = root_path / fname
            base = Path(fname).stem  # 'rmsd_protein' or 'rmsd_complex'

            # Example: ROOT__rmsd_protein.dat, run1__rep0__rmsd_complex.dat
            new_name = f"{rel_str}__{base}.dat"
            dst = out_dir / new_name

            # If there is a collision, add a counter suffix
            suffix = 1
            while dst.exists():
                dst = out_dir / f"{rel_str}__{base}_{suffix}.dat"
                suffix += 1

            shutil.copy2(src, dst)
            count += 1
            print(f"[{count:3d}] {src}  -->  {dst.name}")

    print(f"\nDone. Collected {count} file(s) into: {out_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Collect all 'rmsd_protein.dat' and 'rmsd_complex.dat' files into one folder."
    )
    parser.add_argument(
        "base_dir",
        nargs="?",
        default=".",
        help="Base directory to search (default: current directory).",
    )
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    if not base_dir.is_dir():
        raise SystemExit(f"ERROR: {base_dir} is not a directory.")

    collect_rmsd(base_dir)


if __name__ == "__main__":
    main()
