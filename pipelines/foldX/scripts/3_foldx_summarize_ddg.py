#!/usr/bin/env python3
"""
Summarize FoldX BuildModel pH-scan outputs.

Inputs:
  - per-variant output folders containing Dif_*.fxout files, typically:
      outputs/foldx_buildmodel/per_variant/<variant>/T298/fxout/

Outputs:
  - long table for all parsed files
  - per-temperature wide table
  - per-temperature metrics table
  - per-temperature tree value tables and 2-column mapping files
  - progress report

Optional manifest:
  A TSV with at least one mutation-file column and one sequence/taxon column.
  Accepted header names include:
    - Sequence_ID, sequence_id, taxon, tree_key
    - mutfile, mutation_file, variant_folder
"""

import argparse
import csv
import glob
import math
import os
import re
from collections import Counter, defaultdict

PH_RE = re.compile(r"_pH(\d+p\d+)_", re.IGNORECASE)
TEMP_RE = re.compile(r"(?:^|[\\/])T(\d{3})(?:[\\/]|$)", re.IGNORECASE)


def safe_float(value):
    s = str(value).strip().replace(",", ".")
    s = re.sub(r"[^0-9eE\+\-\.]", "", s)
    try:
        return float(s)
    except Exception:
        return float("nan")


def split_fields(line):
    if "\t" in line:
        return [x.strip() for x in line.rstrip("\n").split("\t")]
    return [x.strip() for x in line.rstrip("\n").split()]


def find_table(lines):
    header_i = None
    headers = None
    for i, line in enumerate(lines):
        low = line.lower()
        if ("pdb" in low) and ("total energy" in low) and ("\t" in line or len(line.split()) > 3):
            header_i = i
            headers = split_fields(line)
            break
    if header_i is None:
        return None, []

    rows = []
    for line in lines[header_i + 1 :]:
        text = line.strip()
        if not text:
            continue
        low = text.lower()
        if low.startswith("#") or low.startswith("foldx") or low.startswith("command") or low.startswith("date"):
            continue
        fields = split_fields(line)
        if len(fields) >= 2:
            rows.append(fields)
    return headers, rows


def total_energy_col(headers):
    for i, h in enumerate(headers):
        if h.strip().lower() == "total energy":
            return i
    for i, h in enumerate(headers):
        low = h.strip().lower()
        if "total" in low and "energy" in low:
            return i
    for i, h in enumerate(headers):
        if h.strip().lower() not in ("pdb", "pdbname", "pdb_name"):
            return i
    return None


def parse_ph(name):
    match = PH_RE.search(name)
    if not match:
        return None
    ph = safe_float(match.group(1).replace("p", "."))
    if ph == ph:
        ph = round(ph, 2)
    return ph


def parse_temp(path_or_name):
    match = TEMP_RE.search(path_or_name)
    if not match:
        return "TUNK"
    return f"T{match.group(1)}"


def mean_sd(values):
    vals = [v for v in values if v == v]
    n = len(vals)
    if n == 0:
        return float("nan"), float("nan"), 0
    mu = sum(vals) / n
    if n == 1:
        return mu, 0.0, 1
    var = sum((v - mu) ** 2 for v in vals) / (n - 1)
    return mu, math.sqrt(var), n


def linreg(points):
    xs = [x for x, y in points if x == x and y == y]
    ys = [y for x, y in points if x == x and y == y]
    n = len(xs)
    if n < 2:
        return float("nan"), float("nan"), float("nan")
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    ssx = sum((x - xbar) ** 2 for x in xs)
    if ssx == 0:
        return float("nan"), float("nan"), float("nan")
    sxy = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
    slope = sxy / ssx
    intercept = ybar - slope * xbar
    ss_tot = sum((y - ybar) ** 2 for y in ys)
    ss_res = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
    r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else float("nan")
    return intercept, slope, r2


def ph_col(ph):
    s = f"{ph:.2f}".rstrip("0").rstrip(".")
    if "." in s:
        a, b = s.split(".", 1)
        return f"ddG_pH{a}p{b}"
    return f"ddG_pH{s}"


def load_manifest(manifest_path):
    if not manifest_path or not os.path.isfile(manifest_path):
        return {}

    with open(manifest_path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return {}

        field_map = {name.strip().lower(): name for name in reader.fieldnames}
        seq_field = None
        mut_field = None
        for candidate in ("sequence_id", "taxon", "tree_key"):
            if candidate in field_map:
                seq_field = field_map[candidate]
                break
        for candidate in ("mutfile", "mutation_file", "variant_folder"):
            if candidate in field_map:
                mut_field = field_map[candidate]
                break
        if seq_field is None or mut_field is None:
            return {}

        mapping = {}
        for row in reader:
            seq = (row.get(seq_field) or "").strip()
            mut = (row.get(mut_field) or "").strip()
            if not mut:
                continue
            key = os.path.basename(mut)
            key = os.path.splitext(key)[0]
            if seq:
                mapping[key] = seq
        return mapping


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default="outputs/foldx_buildmodel/per_variant", help="Per-variant FoldX output root")
    parser.add_argument("--outdir", default="outputs/foldx_summary", help="Directory for summary tables")
    parser.add_argument("--manifest", default="foldx_mutfiles/foldx_mutfiles_manifest.tsv", help="Optional TSV mapping mutation files to tree keys")
    parser.add_argument("--acid-phs", default="5.7,5.5,5.3", help="Comma-separated acid pH points for delta metrics")
    parser.add_argument("--neutral-pref", default="6.9,7.0", help="Preferred neutral pH points, first present value is used")
    return parser.parse_args()


def main():
    args = parse_args()
    root = args.root
    outdir = args.outdir
    acid_phs = [round(float(x), 2) for x in args.acid_phs.split(",") if x.strip()]
    neutral_pref = [round(float(x), 2) for x in args.neutral_pref.split(",") if x.strip()]

    os.makedirs(outdir, exist_ok=True)
    folder_to_tree = load_manifest(args.manifest)

    variant_dirs = sorted([p for p in glob.glob(os.path.join(root, "*")) if os.path.isdir(p)])
    if not variant_dirs:
        raise SystemExit(f"ERROR: No variant folders found under {root}")

    long_rows = []
    wide = defaultdict(dict)
    progress = Counter()
    temps_seen = set()

    for vdir in variant_dirs:
        variant_folder = os.path.basename(vdir)
        tree_key = folder_to_tree.get(variant_folder, variant_folder)

        difs = sorted(glob.glob(os.path.join(vdir, "**", "Dif_*.fxout"), recursive=True))
        if not difs:
            continue

        ph_seen_by_temp = defaultdict(set)

        for dif in difs:
            base = os.path.basename(dif)
            ph = parse_ph(base)
            if ph is None or ph != ph:
                continue
            temp = parse_temp(dif)
            temps_seen.add(temp)

            with open(dif, "r", errors="replace") as handle:
                lines = handle.read().splitlines()
            headers, rows = find_table(lines)
            if not headers or not rows:
                long_rows.append((variant_folder, tree_key, temp, ph, math.nan, math.nan, 0, base, "no_table"))
                continue

            col = total_energy_col(headers)
            if col is None:
                long_rows.append((variant_folder, tree_key, temp, ph, math.nan, math.nan, 0, base, "no_total_energy_col"))
                continue

            values = []
            for row in rows:
                if col < len(row):
                    values.append(safe_float(row[col]))

            mean_val, sd_val, n_models = mean_sd(values)
            ph_seen_by_temp[temp].add(ph)
            long_rows.append((variant_folder, tree_key, temp, ph, mean_val, sd_val, n_models, base, ""))
            wide[(tree_key, temp)][ph] = mean_val

        for temp, seen_phs in ph_seen_by_temp.items():
            progress[(temp, len(seen_phs))] += 1

    if not long_rows:
        raise SystemExit("ERROR: Found variant folders but no parsable Dif_*.fxout tables.")

    long_path = os.path.join(outdir, "foldx_ddg_all_long.tsv")
    with open(long_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["variant_folder", "tree_key", "temp", "pH", "ddG_total_mean", "ddG_total_sd", "n_models", "dif_file", "note"])
        for row in sorted(long_rows, key=lambda x: (x[2], x[1], x[3])):
            writer.writerow(row)

    for temp in sorted(temps_seen):
        keys = sorted([k for (k, t) in wide.keys() if t == temp])
        if not keys:
            continue

        all_phs = sorted({ph for (k, t) in wide.keys() if t == temp for ph in wide[(k, t)].keys()})
        neutral_ph = None
        for candidate in neutral_pref:
            if candidate in all_phs:
                neutral_ph = candidate
                break
        if neutral_ph is None and all_phs:
            neutral_ph = max(all_phs)

        wide_path = os.path.join(outdir, f"foldx_ddg_wide_{temp}.tsv")
        with open(wide_path, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["tree_key", "temp"] + [ph_col(ph) for ph in all_phs])
            for key in keys:
                writer.writerow([key, temp] + [wide[(key, temp)].get(ph, "") for ph in all_phs])

        tree_values_path = os.path.join(outdir, f"foldx_tree_values_{temp}.tsv")
        metrics_path = os.path.join(outdir, f"foldx_metrics_{temp}.tsv")
        tree_maps_dir = os.path.join(outdir, f"tree_maps_{temp}")
        os.makedirs(tree_maps_dir, exist_ok=True)

        map_buffers = {"slope_ddG_per_pH": [], "acid_mean_minus_neutral": []}
        for acid_ph in acid_phs:
            map_buffers[f"delta_{ph_col(acid_ph)}minus_neutral"] = []

        with open(tree_values_path, "w", newline="") as tv_handle, open(metrics_path, "w", newline="") as m_handle:
            tv_writer = csv.writer(tv_handle, delimiter="\t")
            m_writer = csv.writer(m_handle, delimiter="\t")

            tree_cols = ["tree_key", "temp", "neutral_pH_used"] + [ph_col(ph) for ph in all_phs]
            tree_cols += [f"delta_{ph_col(ph)}minus_neutral" for ph in acid_phs]
            tree_cols += ["acid_mean_minus_neutral", "slope_ddG_per_pH", "r2"]
            tv_writer.writerow(tree_cols)

            m_writer.writerow([
                "tree_key", "temp", "neutral_pH_used", "n_pH_points",
                "ddG_min", "ddG_max", "ddG_range",
                "ddG_neutral", "ddG_lowest_pH", "ddG_highest_pH",
                "slope_ddG_per_pH", "intercept", "r2",
            ])

            for key in keys:
                points = sorted((ph, wide[(key, temp)][ph]) for ph in wide[(key, temp)].keys())
                values = [val for _, val in points if val == val]
                n_points = len(points)

                ddg_min = min(values) if values else math.nan
                ddg_max = max(values) if values else math.nan
                ddg_range = (ddg_max - ddg_min) if (ddg_min == ddg_min and ddg_max == ddg_max) else math.nan
                ddg_low = points[0][1] if points else math.nan
                ddg_high = points[-1][1] if points else math.nan
                ddg_neutral = wide[(key, temp)].get(neutral_ph, math.nan) if neutral_ph is not None else math.nan

                intercept, slope, r2 = linreg(points)

                deltas = []
                acid_delta_vals = []
                for acid_ph in acid_phs:
                    acid_val = wide[(key, temp)].get(acid_ph, math.nan)
                    if acid_val == acid_val and ddg_neutral == ddg_neutral:
                        delta = acid_val - ddg_neutral
                        deltas.append(delta)
                        acid_delta_vals.append(delta)
                    else:
                        deltas.append(math.nan)

                acid_mean = sum(acid_delta_vals) / len(acid_delta_vals) if acid_delta_vals else math.nan

                row = [key, temp, neutral_ph]
                row += [wide[(key, temp)].get(ph, "") for ph in all_phs]
                row += deltas
                row += [acid_mean, slope, r2]
                tv_writer.writerow(row)

                m_writer.writerow([
                    key, temp, neutral_ph, n_points,
                    ddg_min, ddg_max, ddg_range,
                    ddg_neutral, ddg_low, ddg_high,
                    slope, intercept, r2,
                ])

                map_buffers["slope_ddG_per_pH"].append((key, slope))
                map_buffers["acid_mean_minus_neutral"].append((key, acid_mean))
                for acid_ph, delta in zip(acid_phs, deltas):
                    map_buffers[f"delta_{ph_col(acid_ph)}minus_neutral"].append((key, delta))

        for metric, pairs in map_buffers.items():
            metric_path = os.path.join(tree_maps_dir, f"{metric}.tsv")
            with open(metric_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["taxon", metric])
                for taxon, value in pairs:
                    if value == value:
                        writer.writerow([taxon, value])

    progress_path = os.path.join(outdir, "progress.txt")
    with open(progress_path, "w") as handle:
        handle.write(f"Total variant folders found: {len(variant_dirs)}\n")
        handle.write(f"Temperatures detected: {', '.join(sorted(temps_seen)) if temps_seen else 'None'}\n")
        handle.write("Variants by number of pH points present (per temp):\n")
        for (temp, count) in sorted(progress.keys(), key=lambda x: (x[0], x[1])):
            handle.write(f"  {temp} | {count} pH points: {progress[(temp, count)]}\n")

    print(f"[Wrote] {long_path}")
    print(f"[Wrote] per-temperature files under: {outdir}")
    print(f"[Wrote] {progress_path}")


if __name__ == "__main__":
    main()
