#!/usr/bin/env python3
"""
Sort Rosetta score.sc-style files by total_score.

Example:
    python 2_sort_relax_scores.py outputs/cart_relax/scores/my_input_relax.sc
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple


def parse_score_file(path: Path) -> List[Tuple[float, str]]:
    models: List[Tuple[float, str]] = []
    with path.open() as handle:
        for line in handle:
            if not line.startswith("SCORE:"):
                continue
            if line.startswith("SCORE: total_score"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                score = float(parts[1])
            except ValueError:
                continue
            model_name = parts[-1]
            models.append((score, model_name))
    return sorted(models, key=lambda x: x[0])


def main() -> None:
    parser = argparse.ArgumentParser(description="Sort Rosetta score files by total_score.")
    parser.add_argument("score_file", help="Path to score.sc or equivalent Rosetta score file")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output text file. Default: <score_file stem>_sorted.txt",
    )
    args = parser.parse_args()

    score_file = Path(args.score_file)
    if not score_file.is_file():
        raise SystemExit(f"ERROR: score file not found: {score_file}")

    output = Path(args.output) if args.output else score_file.with_name(f"{score_file.stem}_sorted.txt")
    models = parse_score_file(score_file)

    if not models:
        raise SystemExit(f"ERROR: no model scores found in {score_file}")

    with output.open("w") as out:
        for score, name in models:
            out.write(f"{score:.3f}\t{name}\n")

    best_score, best_name = models[0]
    print(f"Sorted {len(models)} models")
    print(f"Best model: {best_name}  total_score={best_score:.3f}")
    print(f"Wrote: {output}")


if __name__ == "__main__":
    main()
