#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib import cm
from matplotlib.colorbar import ColorbarBase

APEX_COL = "delta_d130_190"
FLANK_COL = "delta_d150_220"

VARIANT_COLORS = {
    "2021ref": (0.133, 0.467, 0.698),
    "9dip":    (0.906, 0.529, 0.169),
    "V3":      (0.267, 0.627, 0.263),
    "V4":      (0.839, 0.153, 0.157),
    "V5":      (0.569, 0.388, 0.757),
    "V6":      (0.769, 0.612, 0.580),
    "V8":      (0.898, 0.898, 0.898),
}
RANK_LABEL = {1: "2021ref", 2: "9dip", 3: "V3", 4: "V4", 5: "V5", 6: "V6", 8: "V8"}
KEY_RANKS_PANEL_B = [1, 2, 3, 4, 5, 6, 8]
LINEAGE_MARKER = {"A": "o", "B": "s"}
YEAR_PALETTE = "tab10"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--rbs-csv", default="rbs_opening_features.csv")
    p.add_argument("--top-csv", default="top_variants_membership.csv")
    p.add_argument("--lineage-csv", default="tips_lineageA_lineageB_mapping.csv")
    p.add_argument("--output-dir", default="fig_af3_rbs_geometry")
    return p.parse_args()

def acc_from_pdb(pdb_name):
    return str(pdb_name).replace(".pdb", "").strip().upper()

def year_from_date(date_str):
    try:
        return int(str(date_str).split("-")[0])
    except Exception:
        return np.nan

def save_legend(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=300, transparent=True, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rbs = pd.read_csv(args.rbs_csv)
    top = pd.read_csv(args.top_csv)
    linmap = pd.read_csv(args.lineage_csv)

    if "pdb" not in rbs.columns:
        raise SystemExit("ERROR: rbs_opening_features.csv must contain a 'pdb' column")

    rbs["acc"] = rbs["pdb"].map(acc_from_pdb)

    rep_meta = top[["variant_rank", "variant_rep_acc", "variant_rep_header"]].drop_duplicates()
    freq = top.groupby("variant_rank")["sequence_acc"].nunique().rename("n_seq")
    rep_meta = rep_meta.merge(freq, on="variant_rank", how="left")
    rep_meta["variant_rep_acc"] = rep_meta["variant_rep_acc"].astype(str).str.upper()

    rbs_rep = rbs.merge(rep_meta, left_on="acc", right_on="variant_rep_acc", how="inner")

    linmap["accession"] = linmap["accession"].astype(str).str.upper()
    linmap["lineage"] = linmap["lineage"].astype(str).str.upper()
    acc_to_lineage = dict(zip(linmap["accession"], linmap["lineage"]))
    acc_to_date = dict(zip(linmap["accession"], linmap["collection_date"]))

    rbs["lineage"] = rbs["acc"].map(acc_to_lineage)
    rbs["collection_date"] = rbs["acc"].map(acc_to_date)
    rbs["year"] = rbs["collection_date"].map(year_from_date)

    rbs_rep["lineage"] = rbs_rep["variant_rep_acc"].map(acc_to_lineage)
    sel_rep = rbs_rep[rbs_rep["variant_rank"].isin(KEY_RANKS_PANEL_B)].copy()
    sel_rep["variant_label"] = sel_rep["variant_rank"].map(RANK_LABEL)

    rbs[["pdb", "acc", "lineage", "collection_date", "year", APEX_COL, FLANK_COL]].to_csv(
        out_dir / "af3_panelA_rbs_geom_points.csv", index=False
    )
    sel_rep[["variant_rank", "variant_label", "variant_rep_acc", "variant_rep_header", "n_seq", "lineage", APEX_COL, FLANK_COL]].to_csv(
        out_dir / "af3_panelB_apex_flank.csv", index=False
    )

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "font.size": 16,
    })

    years = rbs["year"].dropna().astype(int)
    unique_years = sorted(years.unique().tolist()) if not years.empty else [2022, 2023, 2024, 2025]
    base = cm.get_cmap(YEAR_PALETTE, len(unique_years))
    colors = [base(i) for i in range(len(unique_years))]
    cmap_year = ListedColormap(colors)
    year_to_idx = {y: i for i, y in enumerate(unique_years)}
    rbs["year_idx"] = rbs["year"].map(lambda y: year_to_idx.get(int(y), np.nan) if np.isfinite(y) else np.nan)
    bounds = np.arange(-0.5, len(unique_years) + 0.5, 1.0)
    norm_year = BoundaryNorm(bounds, cmap_year.N)

    figA, axA = plt.subplots(figsize=(4.2, 4.2))
    for lin in ["A", "B"]:
        sub = rbs[rbs["lineage"] == lin]
        if sub.empty:
            continue
        axA.scatter(
            sub[APEX_COL], sub[FLANK_COL],
            s=10, alpha=0.70,
            c=sub["year_idx"].values.astype(float),
            cmap=cmap_year, norm=norm_year,
            marker=LINEAGE_MARKER[lin],
            edgecolors="none", zorder=2
        )

    sub_unknown = rbs[~rbs["lineage"].isin(["A", "B"])]
    if not sub_unknown.empty:
        axA.scatter(sub_unknown[APEX_COL], sub_unknown[FLANK_COL], s=10, alpha=0.35, c="0.8", marker="o", edgecolors="none", zorder=1)

    axA.axvline(0.0, color="0.7", lw=0.7, ls="--", zorder=0)
    axA.axhline(0.0, color="0.7", lw=0.7, ls="--", zorder=0)
    axA.set_xlabel("Δd130–190 (Å)")
    axA.set_ylabel("Δd150–220 (Å)")
    axA.set_title("AF3 RBS geometry map")
    for spine in ("top", "right"):
        axA.spines[spine].set_visible(False)
    figA.tight_layout()
    figA.savefig(out_dir / "af3_panelA_rbs_geom.pdf", transparent=True)
    figA.savefig(out_dir / "af3_panelA_rbs_geom.png", dpi=300, transparent=True)
    plt.close(figA)

    sel_rep["lineage_sort"] = sel_rep["lineage"].map({"A": 0, "B": 1}).fillna(2).astype(int)
    sel_rep = sel_rep.sort_values(["lineage_sort", "variant_rank"]).copy()
    labels = sel_rep["variant_label"].tolist()
    x = np.arange(len(labels))
    d_apex = sel_rep[APEX_COL].values
    d_flank = sel_rep[FLANK_COL].values

    figB, axB = plt.subplots(figsize=(4.6, 3.9))
    axB.plot(x, d_apex, color="0.55", lw=1.6, label="Δd130–190 (apex)")
    axB.plot(x, d_flank, color="0.35", lw=1.6, linestyle="--", label="Δd150–220 (flank)")

    for i, label in enumerate(labels):
        color = VARIANT_COLORS.get(label, (0.2, 0.2, 0.2))
        axB.scatter(x[i], d_apex[i], marker="o", s=60, facecolors=color, edgecolors="black", linewidths=0.6, zorder=5)
        axB.scatter(x[i], d_flank[i], marker="s", s=60, facecolors=color, edgecolors="black", linewidths=0.6, zorder=5)

    lineages = sel_rep["lineage"].tolist()
    if "A" in lineages and "B" in lineages:
        last_A = max(i for i, l in enumerate(lineages) if l == "A")
        first_B = min(i for i, l in enumerate(lineages) if l == "B")
        if first_B == last_A + 1:
            axB.axvline(last_A + 0.5, color="0.6", lw=0.8, ls=(0, (3, 3)), zorder=1)

    axB.axhline(0.0, color="0.7", lw=0.7, ls="--")
    axB.set_xticks(x)
    axB.set_xticklabels(labels, rotation=90)
    axB.set_ylabel("Distance change (Å)")
    axB.set_title("Apex and flank shifts in key variants")
    for spine in ("top", "right"):
        axB.spines[spine].set_visible(False)
    figB.tight_layout()
    figB.savefig(out_dir / "af3_panelB_apex_flank.pdf", transparent=True)
    figB.savefig(out_dir / "af3_panelB_apex_flank.png", dpi=300, transparent=True)
    plt.close(figB)

    figL1, axL1 = plt.subplots(figsize=(2.2, 1.1))
    axL1.axis("off")
    handles_lin = [
        Line2D([0], [0], marker=LINEAGE_MARKER["A"], linestyle="none", markersize=7, markerfacecolor="0.3", markeredgecolor="none", label="Lineage A (circle)"),
        Line2D([0], [0], marker=LINEAGE_MARKER["B"], linestyle="none", markersize=7, markerfacecolor="0.3", markeredgecolor="none", label="Lineage B (square)"),
    ]
    leg1 = axL1.legend(handles=handles_lin, frameon=False, loc="center")
    for t in leg1.get_texts():
        t.set_fontsize(14)
    save_legend(figL1, out_dir / "af3_panelA_lineage_legend.png")

    figL2, axL2 = plt.subplots(figsize=(3.4, 0.7))
    axL2.axis("off")
    cax = figL2.add_axes([0.10, 0.45, 0.80, 0.30])
    cb = ColorbarBase(cax, cmap=cmap_year, norm=norm_year, orientation="horizontal")
    cb.set_ticks(np.arange(len(unique_years)))
    cb.set_ticklabels([str(y) for y in unique_years])
    cb.ax.tick_params(labelsize=12, pad=1.5, length=0)
    save_legend(figL2, out_dir / "af3_panelA_year_legend.png")

    figL3, axL3 = plt.subplots(figsize=(4.4, 1.8))
    axL3.axis("off")
    handles = [
        Line2D([0], [0], color="0.55", lw=1.6, label="Δd130–190 (apex)"),
        Line2D([0], [0], color="0.35", lw=1.6, linestyle="--", label="Δd150–220 (flank)"),
    ] + [
        Line2D([0], [0], marker="o", linestyle="none", markersize=7, markerfacecolor=VARIANT_COLORS[label], markeredgecolor="black", label=label)
        for label in ["2021ref", "9dip", "V3", "V4", "V5", "V6", "V8"]
    ]
    leg3 = axL3.legend(handles=handles, frameon=False, loc="center", ncol=4)
    for t in leg3.get_texts():
        t.set_fontsize(12)
    save_legend(figL3, out_dir / "af3_panelB_legend.png")

    print(f"Wrote AF3 RBS geometry outputs to {out_dir}", flush=True)

if __name__ == "__main__":
    main()
