#!/usr/bin/env python3
"""
Combine MM/GBSA result files from mmgbsa_collect/ into two CSVs:

  - mmgbsa_summary_combined.csv   (one row per run)
  - decom_gbneck2_combined.csv   (one row per residue per run)

Run this script from inside MD_9dip23vs26:
    python collect_mmgbsa_to_csv.py
or:
    python collect_mmgbsa_to_csv.py mmgbsa_collect
"""

import csv
import os
import re
import sys
from typing import Dict, List, Optional


# ---------- helpers to infer metadata from directory names ----------

def infer_linkage(dirname: str) -> str:
    """Return '23' or '26' based on 9dip23/9dip26 prefix, else try _23_/_26_."""
    if dirname.startswith("9dip23"):
        return "23"
    if dirname.startswith("9dip26"):
        return "26"
    m = re.search(r"_(23|26)_", dirname)
    return m.group(1) if m else "NA"


def infer_variant(dirname: str, linkage: str) -> str:
    """
    Infer variant name from the directory name.

      - '2021ref' : if '2021ref' or '2021_ref' present
      - 'V3'...'V8': if pattern '_V3_23_', etc.
      - '9dip'    : anything else containing '9dip'
    """
    d = dirname

    # 2021 reference runs
    if "2021ref" in d or "2021_ref" in d:
        return "2021ref"

    # Top-5 variants: V3/V4/V5/V6/V8
    m = re.search(r"_V(\d+)_" + re.escape(linkage), d)
    if m:
        return f"V{m.group(1)}"

    # 9dip baseline
    if "9dip" in d:
        return "9dip"

    return "UNKNOWN"


def infer_replicate(dirname: str) -> Optional[int]:
    """
    Replicates are encoded as x.y where y is 0/1/2.
    We take the LAST x.y in the dirname and return y as int.
    """
    matches = re.findall(r"(\d+)\.(\d+)", dirname)
    if not matches:
        return None
    _, rep_str = matches[-1]
    try:
        return int(rep_str)
    except ValueError:
        return None


# ---------- parsers for the MMGBSA text files ----------

def parse_mmgbsa_summary(path: str) -> Dict[str, Optional[float]]:
    """
    Parse mmgbsa_gbneck2.dat and extract:
      - delta_g_gas, delta_g_gas_sd, delta_g_gas_sem
      - delta_g_solv, delta_g_solv_sd, delta_g_solv_sem
      - delta_total,  delta_total_sd,  delta_total_sem
    from 'Differences (Complex - Receptor - Ligand)' block.
    """
    result: Dict[str, Optional[float]] = {
        "delta_g_gas": None,
        "delta_g_gas_sd": None,
        "delta_g_gas_sem": None,
        "delta_g_solv": None,
        "delta_g_solv_sd": None,
        "delta_g_solv_sem": None,
        "delta_total": None,
        "delta_total_sd": None,
        "delta_total_sem": None,
    }

    in_diff = False
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not in_diff:
                if s.startswith("Differences (Complex - Receptor - Ligand"):
                    in_diff = True
                continue

            if not s:
                continue

            parts = s.split()
            if len(parts) < 4:
                continue
            try:
                val, std, sem = map(float, parts[-3:])
            except ValueError:
                continue

            if s.startswith("DELTA G gas"):
                result["delta_g_gas"] = val
                result["delta_g_gas_sd"] = std
                result["delta_g_gas_sem"] = sem
            elif s.startswith("DELTA G solv"):
                result["delta_g_solv"] = val
                result["delta_g_solv_sd"] = std
                result["delta_g_solv_sem"] = sem
            elif s.startswith("DELTA TOTAL"):
                result["delta_total"] = val
                result["delta_total_sd"] = std
                result["delta_total_sem"] = sem

    return result


def parse_decomp_file(path: str,
                      variant: str,
                      linkage: str,
                      replicate: Optional[int],
                      dirname: str) -> List[Dict[str, object]]:
    """
    Parse decom_gbneck2.dat and return a list of per-residue dicts.

    Columns:
      variant, linkage, replicate, run_dir,
      residue, resnum, location,
      internal_avg/sd/sem, vdw_avg/sd/sem, elec_avg/sd/sem,
      polar_avg/sd/sem, nonpolar_avg/sd/sem, total_avg/sd/sem
    """
    rows: List[Dict[str, object]] = []
    in_table = False

    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")

            if not in_table:
                if "Residue" in line and "Location" in line and "," in line:
                    in_table = True
                    # skip the "Avg/Std/SEM" line
                    next(f, None)
                continue

            if not line.strip():
                break

            if line.startswith("Residue") or line.startswith("#"):
                continue

            parts = line.split(",")
            if len(parts) < 20:
                continue

            residue_field = parts[0].strip()   # e.g. 'ARG 1481' or 'GLY1006' or '4GB1490'
            location_field = parts[1].strip()  # e.g. 'R ARG 1481' or 'L 4GB1490'

            # --- robust residue / resnum parsing ---
            residue_str = residue_field.strip()
            resname = residue_str
            resnum: Optional[int] = None

            tokens = residue_str.split()
            if len(tokens) >= 2 and tokens[-1].isdigit():
                # e.g. 'ARG 1481' -> 'ARG', 1481 ; 'ROH 1495' -> 'ROH', 1495
                resname = " ".join(tokens[:-1])
                try:
                    resnum = int(tokens[-1])
                except ValueError:
                    resnum = None
            else:
                # No explicit space-separated number; split off trailing digits
                i = len(residue_str) - 1
                while i >= 0 and residue_str[i].isdigit():
                    i -= 1
                if i < len(residue_str) - 1:
                    name_part = residue_str[:i+1]   # up to last non-digit
                    num_part = residue_str[i+1:]    # trailing digits
                    if name_part:
                        resname = name_part
                    try:
                        resnum = int(num_part)
                    except ValueError:
                        resnum = None
                else:
                    # no trailing digits, keep as is
                    resname = residue_str
                    resnum = None

            # location: 'R' or 'L' from first token
            loc_tokens = location_field.split()
            location = loc_tokens[0] if loc_tokens else ""

            def f(idx: int) -> float:
                try:
                    return float(parts[idx])
                except ValueError:
                    return float("nan")

            row = {
                "variant": variant,
                "linkage": linkage,
                "replicate": replicate,
                "run_dir": dirname,
                "residue": resname,
                "resnum": resnum,
                "location": location,
                "internal_avg": f(2),
                "internal_sd": f(3),
                "internal_sem": f(4),
                "vdw_avg": f(5),
                "vdw_sd": f(6),
                "vdw_sem": f(7),
                "elec_avg": f(8),
                "elec_sd": f(9),
                "elec_sem": f(10),
                "polar_avg": f(11),
                "polar_sd": f(12),
                "polar_sem": f(13),
                "nonpolar_avg": f(14),
                "nonpolar_sd": f(15),
                "nonpolar_sem": f(16),
                "total_avg": f(17),
                "total_sd": f(18),
                "total_sem": f(19),
            }
            rows.append(row)

    return rows


# ---------- main driver ----------

def main():
    # root directory with the collected per-run dirs
    if len(sys.argv) > 1:
        mmgbsa_root = sys.argv[1]
    else:
        mmgbsa_root = "mmgbsa_collect"

    if not os.path.isdir(mmgbsa_root):
        raise SystemExit(f"Cannot find directory: {mmgbsa_root}")

    mmgbsa_csv = "mmgbsa_summary_combined.csv"
    decomp_csv = "decom_gbneck2_combined.csv"

    mmgbsa_rows: List[Dict[str, object]] = []
    decomp_rows: List[Dict[str, object]] = []

    for entry in sorted(os.listdir(mmgbsa_root)):
        run_dir = os.path.join(mmgbsa_root, entry)
        if not os.path.isdir(run_dir):
            continue

        dirname = os.path.basename(run_dir)
        linkage = infer_linkage(dirname)
        variant = infer_variant(dirname, linkage)
        replicate = infer_replicate(dirname)

        mmgbsa_file = None
        decomp_file = None

        for fname in os.listdir(run_dir):
            if fname.endswith("mmgbsa_gbneck2.dat"):
                mmgbsa_file = os.path.join(run_dir, fname)
            elif fname.endswith("decom_gbneck2.dat"):
                decomp_file = os.path.join(run_dir, fname)

        if mmgbsa_file is None and decomp_file is None:
            continue

        if mmgbsa_file is not None:
            summary = parse_mmgbsa_summary(mmgbsa_file)
            row = {
                "variant": variant,
                "linkage": linkage,
                "replicate": replicate,
                "run_dir": dirname,
                **summary,
            }
            mmgbsa_rows.append(row)

        if decomp_file is not None:
            decomp_rows.extend(
                parse_decomp_file(
                    decomp_file,
                    variant=variant,
                    linkage=linkage,
                    replicate=replicate,
                    dirname=dirname,
                )
            )

    # ---- write mmgbsa summary CSV ----
    if mmgbsa_rows:
        mm_fields = [
            "variant",
            "linkage",
            "replicate",
            "run_dir",
            "delta_g_gas",
            "delta_g_gas_sd",
            "delta_g_gas_sem",
            "delta_g_solv",
            "delta_g_solv_sd",
            "delta_g_solv_sem",
            "delta_total",
            "delta_total_sd",
            "delta_total_sem",
        ]
        with open(mmgbsa_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=mm_fields)
            writer.writeheader()
            for row in mmgbsa_rows:
                writer.writerow(row)
        print(f"Wrote {len(mmgbsa_rows)} rows to {mmgbsa_csv}")
    else:
        print("No mmgbsa_gbneck2.dat files found.")

    # ---- write decomposition CSV ----
    if decomp_rows:
        decomp_fields = [
            "variant",
            "linkage",
            "replicate",
            "run_dir",
            "residue",
            "resnum",
            "location",
            "internal_avg",
            "internal_sd",
            "internal_sem",
            "vdw_avg",
            "vdw_sd",
            "vdw_sem",
            "elec_avg",
            "elec_sd",
            "elec_sem",
            "polar_avg",
            "polar_sd",
            "polar_sem",
            "nonpolar_avg",
            "nonpolar_sd",
            "nonpolar_sem",
            "total_avg",
            "total_sd",
            "total_sem",
        ]
        with open(decomp_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=decomp_fields)
            writer.writeheader()
            for row in decomp_rows:
                writer.writerow(row)
        print(f"Wrote {len(decomp_rows)} rows to {decomp_csv}")
    else:
        print("No decom_gbneck2.dat files found.")


if __name__ == "__main__":
    main()
