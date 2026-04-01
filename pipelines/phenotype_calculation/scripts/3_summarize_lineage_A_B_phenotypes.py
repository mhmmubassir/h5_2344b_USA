#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
summarize_lineage_phenotypes_noargs.py

Creates lineage-wise phenotype summary tables (median and IQR) for:
  - Lineage A (all variants)
  - Lineage B (all variants)
  - Lineage B "cattle-associated" variants (operational definition: domestic-mammal enriched)
  - Lineage B non-cattle variants

No command-line arguments are required.
"""

import os
import numpy as np
import pandas as pd

# ─── USER SETTINGS ─────────────────────────────────────────────────────────────
REP_VARIANT_CSV = "master_v3_2195_rep_phenotypes_with_variants_sources_with_lineage.csv"
ALLSEQ_CSV      = "master_v3_13000_all_phenotypes_with_variants_sources_with_lineage_working_rows_only.csv"

OUT_DIR         = "phenotype_summary_outputs"   # output folder
THRESH_DM_PROP  = 0.50   # domestic mammal fraction threshold per variant
THRESH_DM_COUNT = 10     # minimum domestic mammal sequences per variant
# ──────────────────────────────────────────────────────────────────────────────

def compute_variant_dm_stats(seq_df, variant_col="variant_num_rep"):
    """Compute per-variant domestic mammal counts and fraction using host == 'DM'."""
    if "host" not in seq_df.columns:
        raise ValueError("ALLSEQ_CSV must contain column: host")
    g = seq_df.groupby(variant_col, dropna=False)
    dm_count = g.apply(lambda x: (x["host"] == "DM").sum())
    total = g.size()
    dm_prop = (dm_count / total).fillna(0)
    out = pd.DataFrame({"DM_count": dm_count, "total_count": total, "DM_prop": dm_prop}).reset_index()
    return out

def summarize_metrics(df, mask, metric_cols):
    """Return per-metric median and IQR (q1, q3) within a mask."""
    rows = []
    sub = df.loc[mask, metric_cols]
    for col in metric_cols:
        s = sub[col].dropna()
        if len(s) == 0:
            rows.append((col, np.nan, np.nan, np.nan, 0))
            continue
        q1 = float(s.quantile(0.25))
        q3 = float(s.quantile(0.75))
        med = float(s.median())
        rows.append((col, med, q1, q3, int(s.shape[0])))
    return pd.DataFrame(rows, columns=["metric", "median", "q1", "q3", "n"]).set_index("metric")

def fmt_iqr(med, q1, q3, decimals=3):
    if pd.isna(med):
        return ""
    return f"{med:.{decimals}f} ({q1:.{decimals}f} to {q3:.{decimals}f})"

def main():
    rep = pd.read_csv(REP_VARIANT_CSV)
    seq = pd.read_csv(ALLSEQ_CSV)

    # Required columns
    for col in ["variant_num", "lineage"]:
        if col not in rep.columns:
            raise ValueError(f"REP_VARIANT_CSV must contain column: {col}")
    for col in ["variant_num_rep", "host"]:
        if col not in seq.columns:
            raise ValueError(f"ALLSEQ_CSV must contain column: {col}")

    # DM enrichment (proxy for cattle association) computed from all sequences
    dm_stats = compute_variant_dm_stats(seq, "variant_num_rep").rename(columns={"variant_num_rep": "variant_num"})
    rep2 = rep.merge(dm_stats, on="variant_num", how="left")
    rep2[["DM_count", "total_count", "DM_prop"]] = rep2[["DM_count", "total_count", "DM_prop"]].fillna(0)

    lineage = rep2["lineage"].astype(str).str.upper()
    mask_A = lineage.eq("A")
    mask_B = lineage.eq("B")
    mask_B_cattle = mask_B & (rep2["DM_prop"] >= THRESH_DM_PROP) & (rep2["DM_count"] >= THRESH_DM_COUNT)
    mask_B_non = mask_B & (~mask_B_cattle)

    # Phenotype columns = numeric columns excluding obvious metadata
    exclude = {"variant_num", "variant_count", "DM_count", "total_count", "DM_prop"}
    exclude |= {c for c in rep2.columns if c.startswith("imputed_")}
    metric_cols = [c for c in rep2.columns if pd.api.types.is_numeric_dtype(rep2[c]) and c not in exclude]

    sum_A = summarize_metrics(rep2, mask_A, metric_cols)
    sum_B = summarize_metrics(rep2, mask_B, metric_cols)
    sum_B_cat = summarize_metrics(rep2, mask_B_cattle, metric_cols)
    sum_B_non = summarize_metrics(rep2, mask_B_non, metric_cols)

    # Wide formatted output
    wide = pd.DataFrame(index=sum_A.index)
    wide["Lineage A (median, IQR)"] = [fmt_iqr(*row) for row in zip(sum_A["median"], sum_A["q1"], sum_A["q3"])]
    wide["Lineage B (median, IQR)"] = [fmt_iqr(*row) for row in zip(sum_B["median"], sum_B["q1"], sum_B["q3"])]
    wide["B cattle-associated variants (median, IQR)"] = [fmt_iqr(*row) for row in zip(sum_B_cat["median"], sum_B_cat["q1"], sum_B_cat["q3"])]
    wide["B non-cattle variants (median, IQR)"] = [fmt_iqr(*row) for row in zip(sum_B_non["median"], sum_B_non["q1"], sum_B_non["q3"])]
    wide["n_A"] = sum_A["n"].astype(int)
    wide["n_B"] = sum_B["n"].astype(int)
    wide["n_B_cattle"] = sum_B_cat["n"].astype(int)
    wide["n_B_non"] = sum_B_non["n"].astype(int)
    wide = wide.reset_index().rename(columns={"index": "metric"})

    # Long numeric output
    def to_long(sum_df, group):
        df = sum_df.reset_index()
        df["group"] = group
        return df[["metric", "group", "median", "q1", "q3", "n"]]

    long = pd.concat([
        to_long(sum_A, "Lineage A"),
        to_long(sum_B, "Lineage B"),
        to_long(sum_B_cat, f"Lineage B cattle-associated (DM_prop>={THRESH_DM_PROP}, DM_count>={THRESH_DM_COUNT})"),
        to_long(sum_B_non, "Lineage B non-cattle"),
    ], ignore_index=True)

    # Variant-level transparency table
    variant_stats = rep2[["variant_num", "lineage", "variant_count", "DM_count", "total_count", "DM_prop"]].copy()
    variant_stats["B_cattle_associated"] = mask_B_cattle.astype(int)

    os.makedirs(OUT_DIR, exist_ok=True)
    wide.to_csv(os.path.join(OUT_DIR, "phenotype_summary_lineageAB_wide.csv"), index=False)
    long.to_csv(os.path.join(OUT_DIR, "phenotype_summary_lineageAB_long.csv"), index=False)
    variant_stats.to_csv(os.path.join(OUT_DIR, "variant_domestic_mammal_enrichment_stats.csv"), index=False)

    print("Wrote outputs to:", OUT_DIR)

if __name__ == "__main__":
    main()
