#!/usr/bin/env python3
"""
make_af3_jsons_trimer.py
----------------------------------------
Generate AlphaFold-3 JSON inputs for **homotrimers** (chain IDs A/B/C) and
write `input.lst` for the Slurm workflow. Filenames will use only the first
three underscore-separated parts of the FASTA header (e.g. EPI_ISL_15078246).

Usage examples
  python make_af3_jsons_trimer.py seqs.fa
  python make_af3_jsons_trimer.py seq_dir/
  python make_af3_jsons_trimer.py --fasta seqs.fa
  python make_af3_jsons_trimer.py --fasta_dir seq_dir/ --out_dir jsons --seed 2
"""

import argparse
import json
import re
import sys
from pathlib import Path

TRIMER_IDS = ["A", "B", "C"]    # chain IDs for homotrimer


def safe_stem(name: str) -> str:
    """Filesystem-safe filename stem."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def read_fasta_records(path: Path):
    """Yield (header, sequence) tuples from a FASTA file."""
    header, seq_frag = None, []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq_frag)
                header = line[1:].split()[0]
                seq_frag = []
            else:
                seq_frag.append(line)
        if header:
            yield header, "".join(seq_frag)


def records_from_source(src: Path):
    if src.is_file():
        yield from read_fasta_records(src)
    elif src.is_dir():
        for fa in sorted(src.glob("*.fa*")):
            yield from read_fasta_records(fa)
    else:
        sys.exit(f"ERROR: {src} not found or not FASTA.")


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--fasta", type=Path, help="multi-FASTA file")
    g.add_argument("--fasta_dir", type=Path, help="directory of FASTA files")
    ap.add_argument("path", nargs="?", help="FASTA file or directory")
    ap.add_argument("--out_dir", type=Path, default=Path("jsons"),
                    help="[default: jsons]")
    ap.add_argument("--seed", type=int, default=1,
                    help="[default: 1] modelSeeds value")
    args = ap.parse_args()

    # resolve input source
    if not (args.fasta or args.fasta_dir):
        if not args.path:
            ap.error("Provide a FASTA file/dir or use --fasta / --fasta_dir.")
        p = Path(args.path)
        if p.is_dir():
            args.fasta_dir = p
        else:
            args.fasta = p

    src = args.fasta if args.fasta else args.fasta_dir
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    input_lines = []
    for rec_id, seq in records_from_source(src):
        # Trim everything after a pipe ("|") if present, then grab first three
        base_id = rec_id.split("|")[0]
        parts = base_id.split("_")
        accession = "_".join(parts[:3]) if len(parts) >= 3 else base_id
        stem = safe_stem(accession)
        fname = f"{stem}.json"
        json_path = out_dir / fname

        payload = {
            "name": rec_id,
            "sequences": [
                {
                    "protein": {
                        "id": TRIMER_IDS,
                        "sequence": seq
                    }
                }
            ],
            "modelSeeds": [args.seed],
            "dialect": "alphafold3",
            "version": 1
        }

        with json_path.open("w") as fh:
            json.dump(payload, fh, indent=2)
        input_lines.append(fname)

    # write input.lst in project root
    with open("input.lst", "w") as fh:
        fh.write("\n".join(input_lines) + "\n")

    print(f"✔ Wrote {len(input_lines)} JSON files to {out_dir.resolve()}")
    print("✔ Created input.lst")


if __name__ == "__main__":
    main()
