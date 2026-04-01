#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

H3_LOOPS = {
    "130": (130, 138),
    "150": (150, 158),
    "190": (187, 195),
    "220": (218, 225),
}
H3_TIPS = {
    "130": 135,
    "150": 156,
    "190": 193,
    "220": 222,
}

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pdb-dir", default="best_pdbs")
    p.add_argument("--chain", default="A")
    p.add_argument("--numbering-csv", default="compare_numbering.csv")
    p.add_argument("--reference-substring", default="EPI_ISL_18133029")
    p.add_argument("--output-csv", default="rbs_opening_features.csv")
    return p.parse_args()

def log(msg):
    print(msg, flush=True)

def load_crystal_to_rosetta_map(csv_path: Path, chain: str):
    mapping = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["Crystal_chain"].strip() != chain:
                continue
            try:
                crystal_num = int(row["Crystal_resNum"])
                rosetta_idx = int(row["Rosetta_index"])
            except ValueError:
                continue
            if crystal_num not in mapping:
                mapping[crystal_num] = rosetta_idx
    if not mapping:
        raise SystemExit(f"ERROR: no numbering map loaded from {csv_path} for chain {chain}")
    return mapping

def convert_range(start, end, mapping):
    indices = [mapping[c] for c in range(start, end + 1) if c in mapping]
    if not indices:
        raise SystemExit(f"ERROR: no numbering map for crystal range {start}-{end}")
    return min(indices), max(indices)

def map_tips(h3_tips, mapping):
    out = {}
    for loop_name, h3_pos in h3_tips.items():
        if h3_pos in mapping:
            out[loop_name] = mapping[h3_pos]
        else:
            log(f"WARN: missing numbering map for loop-tip H3 {h3_pos} ({loop_name})")
    return out

def parse_ca_coords(pdb_path: Path, chain_id: str):
    coords = {}
    with open(pdb_path) as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[21].strip() != chain_id:
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                resnum = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            coords[resnum] = np.array([x, y, z], dtype=float)
    return coords

def centroid_for_range(ca_coords, start_idx, end_idx):
    pts = [ca_coords[r] for r in range(start_idx, end_idx + 1) if r in ca_coords]
    if not pts:
        return None
    return np.vstack(pts).mean(axis=0)

def vec_dist(a, b):
    if a is None or b is None:
        return np.nan
    return float(np.linalg.norm(a - b))

def collect_geometry_for_pdb(pdb_path, chain_id, loop_ranges_idx, tips_idx, ref_dists, ref_tip_coords):
    ca_coords = parse_ca_coords(pdb_path, chain_id)
    if not ca_coords:
        return None

    loop_centroids = {name: centroid_for_range(ca_coords, s, e) for name, (s, e) in loop_ranges_idx.items()}
    d130_190 = vec_dist(loop_centroids.get("130"), loop_centroids.get("190"))
    d130_220 = vec_dist(loop_centroids.get("130"), loop_centroids.get("220"))
    d150_220 = vec_dist(loop_centroids.get("150"), loop_centroids.get("220"))

    row = {
        "d130_190": d130_190,
        "d130_220": d130_220,
        "d150_220": d150_220,
        "delta_d130_190": d130_190 - ref_dists["d130_190"] if not np.isnan(d130_190) else np.nan,
        "delta_d130_220": d130_220 - ref_dists["d130_220"] if not np.isnan(d130_220) else np.nan,
        "delta_d150_220": d150_220 - ref_dists["d150_220"] if not np.isnan(d150_220) else np.nan,
    }

    for loop_name, tip_resnum in tips_idx.items():
        ref_tip = ref_tip_coords.get(loop_name)
        tip = ca_coords.get(tip_resnum)
        row[f"tip_disp_{loop_name}"] = vec_dist(tip, ref_tip) if (ref_tip is not None and tip is not None) else np.nan

    return row

def main():
    args = parse_args()
    pdb_dir = Path(args.pdb_dir)
    numbering_csv = Path(args.numbering_csv)
    out_csv = Path(args.output_csv)

    if not numbering_csv.is_file():
        raise SystemExit(f"ERROR: numbering CSV not found: {numbering_csv}")
    if not pdb_dir.is_dir():
        raise SystemExit(f"ERROR: PDB directory not found: {pdb_dir}")

    mapping = load_crystal_to_rosetta_map(numbering_csv, args.chain)
    loop_ranges_idx = {name: convert_range(s, e, mapping) for name, (s, e) in H3_LOOPS.items()}
    tips_idx = map_tips(H3_TIPS, mapping)

    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    if not pdb_files:
        raise SystemExit(f"ERROR: no PDB files found in {pdb_dir}")

    ref_pdb = next((p for p in pdb_files if args.reference_substring in p.name), None)
    if ref_pdb is None:
        raise SystemExit(f"ERROR: no reference PDB containing '{args.reference_substring}' in {pdb_dir}")

    log(f"Reference PDB: {ref_pdb.name}")

    ref_ca = parse_ca_coords(ref_pdb, args.chain)
    if not ref_ca:
        raise SystemExit("ERROR: no CA coordinates found in reference PDB")

    ref_loop_centroids = {name: centroid_for_range(ref_ca, s, e) for name, (s, e) in loop_ranges_idx.items()}
    ref_dists = {
        "d130_190": vec_dist(ref_loop_centroids.get("130"), ref_loop_centroids.get("190")),
        "d130_220": vec_dist(ref_loop_centroids.get("130"), ref_loop_centroids.get("220")),
        "d150_220": vec_dist(ref_loop_centroids.get("150"), ref_loop_centroids.get("220")),
    }
    ref_tip_coords = {name: ref_ca[idx] for name, idx in tips_idx.items() if idx in ref_ca}

    rows = []
    for i, pdb_path in enumerate(pdb_files, start=1):
        feats = collect_geometry_for_pdb(pdb_path, args.chain, loop_ranges_idx, tips_idx, ref_dists, ref_tip_coords)
        if feats is None:
            continue
        row = {"pdb": pdb_path.name}
        row.update(feats)
        rows.append(row)
        if i % 50 == 0:
            log(f"Processed {i} PDBs")

    if not rows:
        raise SystemExit("ERROR: no RBS opening features computed")

    dvals = [r["d130_190"] for r in rows if not np.isnan(r["d130_190"])]
    if dvals:
        q1, q2 = np.quantile(dvals, [1.0/3.0, 2.0/3.0])
        for r in rows:
            d = r["d130_190"]
            if np.isnan(d):
                r["rbs_geom_state"] = "NA"
            elif d <= q1:
                r["rbs_geom_state"] = "closed"
            elif d <= q2:
                r["rbs_geom_state"] = "intermediate"
            else:
                r["rbs_geom_state"] = "open"
    else:
        for r in rows:
            r["rbs_geom_state"] = "NA"

    fieldnames = ["pdb"] + [k for k in rows[0] if k != "pdb"]
    with open(out_csv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    log(f"Wrote {len(rows)} rows to {out_csv}")

if __name__ == "__main__":
    main()
