#!/usr/bin/env python3
"""
Pick the lowest-scoring relaxed PDB from a Rosetta score file and copy it to a clean destination.

Example:
    python 3_select_best_relaxed_model.py \
        outputs/cart_relax/scores/my_input_relax.sc \
        outputs/cart_relax/models/my_input \
        -o outputs/cart_relax/best_models/my_input_relaxed_best.pdb
"""

from __future__ import annotations

import argparse
import shutil
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
            models.append((score, parts[-1]))
    return sorted(models, key=lambda x: x[0])


def resolve_model_path(model_dir: Path, description: str) -> Path:
    candidates = [
        model_dir / description,
        model_dir / f"{description}.pdb",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    matches = sorted(model_dir.glob(f"{description}*.pdb"))
    if matches:
        return matches[0]
    raise FileNotFoundError(f"Could not find model '{description}' in {model_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Select the best relaxed Rosetta model.")
    parser.add_argument("score_file", help="Rosetta score file, for example my_relax.sc")
    parser.add_argument("model_dir", help="Directory containing the corresponding relaxed PDB files")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Destination PDB path. Default: <score stem>_best.pdb next to the score file",
    )
    args = parser.parse_args()

    score_file = Path(args.score_file)
    model_dir = Path(args.model_dir)

    if not score_file.is_file():
        raise SystemExit(f"ERROR: score file not found: {score_file}")
    if not model_dir.is_dir():
        raise SystemExit(f"ERROR: model directory not found: {model_dir}")

    models = parse_score_file(score_file)
    if not models:
        raise SystemExit(f"ERROR: no model scores found in {score_file}")

    best_score, best_description = models[0]
    best_model = resolve_model_path(model_dir, best_description)

    output = Path(args.output) if args.output else score_file.with_name(f"{score_file.stem}_best.pdb")
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(best_model, output)

    print(f"Best model description : {best_description}")
    print(f"Best model score       : {best_score:.3f}")
    print(f"Source model           : {best_model}")
    print(f"Copied to              : {output}")


if __name__ == "__main__":
    main()
