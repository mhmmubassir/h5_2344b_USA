#!/usr/bin/env python3
"""
collect_rmsf_protein.py

Recursively walk through a base directory, find all files named
'rmsf_protein.dat', and copy them into a folder 'all_rmsf_protein'.

Each copied file will be renamed to include its relative path so you can
see where it came from.

Usage:
    python collect_rmsf_protein.py /path/to/base_dir

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

    # Split into components
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



def collect_rmsf(base_dir: Path, out_dir_name: str = "all_rmsf_protein") -> None:
    base_dir = base_dir.resolve()
    out_dir = base_dir / out_dir_name
    out_dir.mkdir(exist_ok=True)

    print(f"Base directory : {base_dir}")
    print(f"Output directory: {out_dir}\n")

    count = 0

    for root, dirs, files in os.walk(base_dir):
        # Don't recurse into the output folder itself
        if out_dir_name in dirs:
            dirs.remove(out_dir_name)

        if "rmsf_protein.dat" not in files:
            continue

        root_path = Path(root)
        rel_path = root_path.relative_to(base_dir)
        rel_str = sanitize_rel_path(str(rel_path))

        src = root_path / "rmsf_protein.dat"

        # New name includes the relative path so you know the origin
        # Example: ROOT__rmsf_protein.dat  or  run1__rep0__rmsf_protein.dat
        new_name = f"{rel_str}__rmsf_protein.dat"
        dst = out_dir / new_name

        # If there is a collision, add a counter suffix
        suffix = 1
        while dst.exists():
            dst = out_dir / f"{rel_str}__rmsf_protein_{suffix}.dat"
            suffix += 1

        shutil.copy2(src, dst)
        count += 1
        print(f"[{count:3d}] {src}  -->  {dst.name}")

    print(f"\nDone. Collected {count} file(s) into: {out_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Collect all 'rmsf_protein.dat' files into one folder."
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

    collect_rmsf(base_dir)


if __name__ == "__main__":
    main()
