#!/usr/bin/env python3
"""
4_benchmark_rosetta_vs_wetlab.py

Benchmark Rosetta and Flex ddG outputs against experimental ΔΔG values.

Input folder should contain:
- experimental_ddg.csv
- Rosetta final CSVs such as:
    9dip_bn16_bb1_23_vs_HA_final_ddg_avg.csv
    9dip_bn16_bb1_23_vs_HA_final_ddg_min.csv
    9dio_bn16_bb1_ori_26_vs_HA_final_ddg_avg.csv
- Flex ddG summary CSVs such as:
    flexddg_wCrelax_summary_9dip23.csv
    flexddg_woCrelax_summary_9dip26.csv

Outputs:
    results/
        all_predictions_merged.csv
        stats_summary.tsv
        stats_overall.csv
        scatter_23.png
        scatter_26.png
        correlation_23.png
        correlation_26.png
        performance_overall.png
        scatter_<protocol>_<linkage>.png
"""

from __future__ import annotations

import argparse
import re
import warnings
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st
from matplotlib.ticker import MaxNLocator

warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)


def norm_mut(mut: str) -> str:
    return re.sub(r"__+", "_", str(mut).replace("/", "_").strip().upper())


def first_numeric(df: pd.DataFrame) -> Optional[str]:
    cols = df.select_dtypes(include="number").columns.tolist()
    return cols[0] if cols else None


def fisher_ci(r: float, n: int) -> Tuple[float, float]:
    if n < 4 or not np.isfinite(r) or abs(r) >= 1:
        return np.nan, np.nan
    z = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    return np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)


def parse_filename(path: Path) -> Tuple[str, Optional[str]]:
    stem = path.stem
    lower = stem.lower()

    if "flexddg" in lower:
        if "wcrelax" in lower:
            protocol = "flexddg_wCrelax"
        elif "wocrelax" in lower:
            protocol = "flexddg_woCrelax"
        else:
            protocol = "flexddg"
        m = re.search(r"(23|26)", lower)
        linkage = m.group(1) if m else None
        return protocol, linkage

    tokens = stem.split("_")
    system = tokens[0] if tokens else "unknown"
    score = next((t for t in tokens if t in {"bn16", "cr15"}), "")
    bb = next((t for t in tokens if t in {"bb1", "bb2", "bb3"}), "")
    replicate = next((t for t in tokens if t in {"ori", "rvmut"}), None)
    linkage = next((t for t in tokens if t in {"23", "26"}), None)

    if "min" in tokens:
        flavor = "min"
    elif "avg" in tokens:
        flavor = "avg"
    else:
        flavor = ""

    parts = [p for p in [system, score, bb, replicate, flavor] if p]
    protocol = "_".join(parts)
    return protocol, linkage


def read_prediction(csv_file: Path, exp_csv: Path) -> Optional[pd.DataFrame]:
    if csv_file.name == exp_csv.name:
        return None
    lower = csv_file.name.lower()
    if "haonly" in lower:
        return None

    try:
        df = pd.read_csv(csv_file)
    except Exception:
        return None

    protocol, linkage = parse_filename(csv_file)

    score_col = None
    for candidate in ("ddG_REU", "final_ddG", "ddG", "ddg"):
        if candidate in df.columns:
            score_col = candidate
            break
    if score_col is None:
        score_col = first_numeric(df)
    if score_col is None or "mutation" not in df.columns:
        return None

    out = df[["mutation", score_col]].rename(columns={score_col: "pred_ddG"}).copy()
    out["mutation"] = out["mutation"].map(norm_mut)
    out["pred_ddG"] = pd.to_numeric(out["pred_ddG"], errors="coerce")
    out = out[np.isfinite(out["pred_ddG"])].copy()
    if out.empty:
        return None

    out["protocol"] = protocol
    out["linkage"] = linkage
    return out


def get_color(protocol: str) -> str:
    p = protocol.lower()
    if p == "flexddg_wcrelax":
        return "#C0392B"
    if p == "flexddg_wocrelax":
        return "#E74C3C"
    if p.startswith("9dip"):
        if "_bb1" in p:
            return "#1F77B4"
        if "_bb2" in p:
            return "#2CA02C"
        if "_bb3" in p:
            return "#17BECF"
    if p.startswith("9dio"):
        if "_bb1" in p:
            return "#FF7F0E"
        if "_bb2" in p:
            return "#9467BD"
        if "_bb3" in p:
            return "#8C564B"
    return "#7F7F7F"


def overview_scatter(preds: pd.DataFrame, exp_long: pd.DataFrame, outdir: Path, linkage: str) -> None:
    sub = preds[preds["linkage"] == linkage]
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    ax.plot([-6, 6], [-6, 6], "--", lw=0.6, color="gray", zorder=0)
    for protocol, grp in sub.groupby("protocol"):
        merged = grp.merge(exp_long[exp_long["linkage"] == linkage], on=["mutation", "linkage"], how="inner")
        ax.scatter(merged["pred_ddG"], merged["exp_ddG"], s=12, alpha=0.8, label=protocol)
    ax.set_xlabel("Predicted ΔΔG")
    ax.set_ylabel("Experimental ΔΔG (kcal/mol)")
    ax.set_title(f"Overview scatter, linkage {linkage}")
    ax.legend(fontsize=4, ncol=2, loc="best")
    fig.tight_layout()
    fig.savefig(outdir / f"scatter_{linkage}.png", dpi=300)
    plt.close(fig)


def correlation_bar(stats: pd.DataFrame, outdir: Path, linkage: str) -> None:
    sub = stats[stats["linkage"] == linkage].sort_values("pearson_r", ascending=False).reset_index(drop=True)
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 3), dpi=300)
    ax.bar(
        sub.index,
        sub["pearson_r"],
        yerr=[sub["pearson_r"] - sub["r_CI_lo"], sub["r_CI_hi"] - sub["pearson_r"]],
        capsize=2,
    )
    ax.set_xticks(sub.index)
    ax.set_xticklabels(sub["protocol"], rotation=90)
    ax.set_ylabel("Pearson r")
    ax.set_title(f"Correlation vs wet lab, linkage {linkage}")
    ax.yaxis.set_major_locator(MaxNLocator(5))
    fig.tight_layout()
    fig.savefig(outdir / f"correlation_{linkage}.png", dpi=300)
    plt.close(fig)


def per_protocol_plots(preds: pd.DataFrame, exp_long: pd.DataFrame, outdir: Path) -> None:
    for (protocol, linkage), grp in preds.groupby(["protocol", "linkage"]):
        merged = grp.merge(exp_long, on=["mutation", "linkage"], how="inner").dropna()
        if len(merged) < 2:
            continue
        fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
        ax.plot([-6, 6], [-6, 6], "--", lw=0.6, color="gray")
        ax.scatter(merged["pred_ddG"], merged["exp_ddG"], s=18)
        r, _ = st.pearsonr(merged["pred_ddG"], merged["exp_ddG"])
        ax.text(0.05, 0.9, f"r = {r:.2f}", transform=ax.transAxes)
        ax.set_xlabel("Predicted ΔΔG")
        ax.set_ylabel("Experimental ΔΔG (kcal/mol)")
        ax.set_title(f"{protocol} ({linkage})")
        fig.tight_layout()
        fig.savefig(outdir / f"scatter_{protocol}_{linkage}.png", dpi=300)
        plt.close(fig)


def performance_overall(overall: pd.DataFrame, outdir: Path) -> None:
    if overall.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 4), dpi=300)
    for _, row in overall.iterrows():
        x = row["rmse"]
        y = row["pearson_r"]
        color = get_color(row["protocol"])
        ax.scatter(x, y, s=40, color=color, edgecolor="k", linewidth=0.3, zorder=2)
        ax.text(
            x,
            y,
            row["protocol"],
            fontsize=4,
            color=color,
            ha="left",
            va="top",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=0.2),
            zorder=3,
        )
    ax.set_xlabel("RMSE (kcal/mol)")
    ax.set_ylabel("Pearson r")
    ax.set_title("Overall protocol performance")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(float(overall["rmse"].min()) - 0.5, float(overall["rmse"].max()) + 0.5)
    ax.set_ylim(float(overall["pearson_r"].min()) - 0.05, float(overall["pearson_r"].max()) + 0.05)
    fig.tight_layout()
    fig.savefig(outdir / "performance_overall.png", dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark Rosetta and Flex ddG against wet-lab ΔΔG values.")
    parser.add_argument("--indir", default="statistical_correlation", help="Input folder containing experimental_ddg.csv and protocol CSVs.")
    parser.add_argument("--outdir", default=None, help="Output folder for stats and figures. Default: <indir>/results")
    args = parser.parse_args()

    indir = Path(args.indir).resolve()
    outdir = Path(args.outdir).resolve() if args.outdir else indir / "results"
    outdir.mkdir(parents=True, exist_ok=True)

    exp_csv = indir / "experimental_ddg.csv"
    if not exp_csv.exists():
        raise SystemExit(f"❌ Missing experimental file: {exp_csv}")

    exp = pd.read_csv(exp_csv).rename(columns={"exp_ddG_23": "23", "exp_ddG_26": "26"})
    exp["mutation"] = exp["mutation"].map(norm_mut)
    exp_long = exp.melt(id_vars="mutation", var_name="linkage", value_name="exp_ddG")

    frames = []
    for csv_file in sorted(indir.glob("*.csv")):
        pred = read_prediction(csv_file, exp_csv)
        if pred is not None:
            frames.append(pred)

    if not frames:
        raise SystemExit("❌ No valid prediction CSVs found in the benchmarking folder.")

    preds = pd.concat(frames, ignore_index=True)
    merged_all = preds.merge(exp_long, on=["mutation", "linkage"], how="inner")
    merged_all.to_csv(outdir / "all_predictions_merged.csv", index=False)

    rows_link = []
    rows_all = []

    for (protocol, linkage), grp in preds.groupby(["protocol", "linkage"]):
        merged = grp.merge(exp_long, on=["mutation", "linkage"], how="inner").dropna()
        if len(merged) < 2:
            continue
        x = merged["pred_ddG"]
        y = merged["exp_ddG"]
        r, p = st.pearsonr(x, y)
        rho, p2 = st.spearmanr(x, y)
        rmse = float(np.sqrt(((x - y) ** 2).mean()))
        lo, hi = fisher_ci(r, len(merged))
        rows_link.append([protocol, linkage, len(merged), r, p, rho, p2, rmse, lo, hi])

    for protocol, grp in preds.groupby("protocol"):
        merged = grp.merge(exp_long, on=["mutation", "linkage"], how="inner").dropna()
        if len(merged) < 2:
            continue
        x = merged["pred_ddG"]
        y = merged["exp_ddG"]
        if len(set(x)) > 1:
            r, p = st.pearsonr(x, y)
            rho, p2 = st.spearmanr(x, y)
        else:
            r, p, rho, p2 = np.nan, np.nan, np.nan, np.nan
        rmse = float(np.sqrt(((x - y) ** 2).mean()))
        lo, hi = fisher_ci(r, len(merged))
        rows_all.append([protocol, "all", len(merged), r, p, rho, p2, rmse, lo, hi])

    stats = pd.DataFrame(
        rows_link + rows_all,
        columns=["protocol", "linkage", "N", "pearson_r", "pearson_p", "spearman_rho", "spearman_p", "rmse", "r_CI_lo", "r_CI_hi"],
    )
    stats.to_csv(outdir / "stats_summary.tsv", sep="\t", index=False)

    overall = stats.query("linkage == 'all'").sort_values("rmse").reset_index(drop=True)
    overall.to_csv(outdir / "stats_overall.csv", index=False)

    plt.rcParams.update({
        "font.size": 6,
        "font.family": "sans-serif",
        "axes.linewidth": 0.6,
        "legend.frameon": False,
    })

    for linkage in ("23", "26"):
        overview_scatter(preds, exp_long, outdir, linkage)
        correlation_bar(stats, outdir, linkage)
    per_protocol_plots(preds, exp_long, outdir)
    performance_overall(overall, outdir)

    print(f"✅ Benchmarking finished. Results: {outdir}")
    if not overall.empty:
        print("\nTop overall performers:")
        print(overall[["protocol", "N", "pearson_r", "rmse"]].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
