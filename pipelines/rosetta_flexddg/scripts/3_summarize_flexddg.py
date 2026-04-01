#!/usr/bin/env python3
"""
Summarize Rosetta Flex ddG database outputs into a single CSV table.

Default input directory : outputs/ddg_db3
Default output CSV      : outputs/flexddg_summary.csv

Example:
    python summarize_flexddg.py
    python summarize_flexddg.py --input-dir outputs/ddg_db3 --output results.csv
"""

from __future__ import annotations

import argparse
import csv
import sqlite3
from pathlib import Path
from statistics import mean
from typing import Dict, List, Optional

REQUIRED_STATES = ("bound_wt", "unbound_wt", "bound_mut", "unbound_mut")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("outputs/ddg_db3"),
        help="Directory containing *_ddG.db3 files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("outputs/flexddg_summary.csv"),
        help="Output CSV file.",
    )
    parser.add_argument(
        "--sort-by",
        choices=("mutation", "ddG_REU"),
        default="ddG_REU",
        help="Sort output rows by mutation tag or ddG_REU.",
    )
    return parser.parse_args()


def extract_scores(db_path: Path) -> Optional[Dict[str, object]]:
    mutation_tag = db_path.name.replace("_ddG.db3", "")

    with sqlite3.connect(str(db_path)) as conn:
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()

        row = cursor.execute(
            "SELECT score_type_id FROM score_types WHERE score_type_name='total_score'"
        ).fetchone()
        if row is None:
            return None
        total_score_id = row[0]

        query = """
            SELECT b.name AS state, s.score_value
            FROM structure_scores s
            JOIN batches b USING(batch_id)
            WHERE s.score_type_id = ?
        """
        rows = cursor.execute(query, (total_score_id,)).fetchall()

    state_to_values: Dict[str, List[float]] = {}
    for row in rows:
        state = row["state"].replace("_dbreport", "")
        state_to_values.setdefault(state, []).append(float(row["score_value"]))

    if not all(state in state_to_values for state in REQUIRED_STATES):
        return None

    averages = {state: mean(values) for state, values in state_to_values.items()}
    ddg_reu = (averages["bound_mut"] - averages["unbound_mut"]) - (
        averages["bound_wt"] - averages["unbound_wt"]
    )

    return {
        "mutation": mutation_tag,
        "ddG_REU": round(ddg_reu, 3),
        "bound_wt_avg": round(averages["bound_wt"], 3),
        "unbound_wt_avg": round(averages["unbound_wt"], 3),
        "bound_mut_avg": round(averages["bound_mut"], 3),
        "unbound_mut_avg": round(averages["unbound_mut"], 3),
        "n_bound_wt": len(state_to_values["bound_wt"]),
        "n_unbound_wt": len(state_to_values["unbound_wt"]),
        "n_bound_mut": len(state_to_values["bound_mut"]),
        "n_unbound_mut": len(state_to_values["unbound_mut"]),
    }


def write_csv(rows: List[Dict[str, object]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "mutation",
        "ddG_REU",
        "bound_wt_avg",
        "unbound_wt_avg",
        "bound_mut_avg",
        "unbound_mut_avg",
        "n_bound_wt",
        "n_unbound_wt",
        "n_bound_mut",
        "n_unbound_mut",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def print_table(rows: List[Dict[str, object]]) -> None:
    headers = ["mutation", "ddG_REU", "bound_wt_avg", "unbound_wt_avg", "bound_mut_avg", "unbound_mut_avg"]
    widths = {header: len(header) for header in headers}

    for row in rows:
        for header in headers:
            widths[header] = max(widths[header], len(str(row[header])))

    line = "  ".join(header.ljust(widths[header]) for header in headers)
    print(line)
    print("  ".join("-" * widths[header] for header in headers))
    for row in rows:
        print("  ".join(str(row[header]).ljust(widths[header]) for header in headers))


def main() -> None:
    args = parse_args()

    if not args.input_dir.is_dir():
        raise SystemExit(f"❌ Input directory not found: {args.input_dir}")

    db_files = sorted(args.input_dir.glob("*_ddG.db3"))
    if not db_files:
        raise SystemExit(f"❌ No *_ddG.db3 files found in: {args.input_dir}")

    rows: List[Dict[str, object]] = []
    skipped: List[Path] = []

    for db_file in db_files:
        result = extract_scores(db_file)
        if result is None:
            skipped.append(db_file)
            continue
        rows.append(result)

    if not rows:
        raise SystemExit("❌ No valid Flex ddG summaries could be extracted.")

    if args.sort_by == "ddG_REU":
        rows.sort(key=lambda row: float(row["ddG_REU"]))
    else:
        rows.sort(key=lambda row: str(row["mutation"]))

    write_csv(rows, args.output)

    print(f"✅ Saved summary CSV: {args.output}")
    print()
    print_table(rows)

    if skipped:
        print()
        print("⚠️ Skipped files with missing required states:")
        for path in skipped:
            print(f"  - {path}")


if __name__ == "__main__":
    main()
