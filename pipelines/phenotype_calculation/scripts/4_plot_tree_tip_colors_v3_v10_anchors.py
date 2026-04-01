#!/usr/bin/env python3
"""
5_tree_colorTips_overall_and_eachPhenotype_V3V10Anchors_py36.py

Tree tip coloring using the V3/V10-anchor scoring outputs:
  - variant2195_V3V10Anchors_centroid_direction_overall.csv
  - variant2195_V3V10Anchors_centroid_direction_long.csv

Generates:
  (1) Tip-colored tree using overall_signed_mean_lvl
  (2) Tip-colored tree using overall_signed_mean_lvl_nonzero
  (3) Tip-colored trees for EACH individual phenotype using per-phenotype signed_lvl (-4..+4)

Color meaning:
  - Negative values  -> blue (closer to V3 "avian" anchor)
  - Positive values  -> red  (closer to V10 "cattle" anchor)
  - 0                -> white/neutral (no-change band)

REQUIRED INPUTS (place in the same directory where you run this script):
  1) timetree.nexus
  2) variant2195_V3V10Anchors_centroid_direction_overall.csv
  3) variant2195_V3V10Anchors_centroid_direction_long.csv
  4) A full→rep mapping source (any ONE of these works; needs columns full_acc, rep_acc):
       - master_v3_13000_all_phenotypes_with_variants_sources_with_lineage_working_rows_only.csv
       - master_v3_13000_all_phenotypes_with_variants_sources_with_lineage.csv
       - master_v3_13000_all_phenotypes_with_variants_sources.csv

OPTIONAL fallback mapping (only if mapping CSV is missing/unusable):
  - 2344b_US_h5n1_1jan21to4may25_TR2_IQsnDNA_13000_woo_aligned.aa.fasta
  - 2344b_US_h5n1_1jan21to4may25_TR2_IQsnDNA_13000_woo_aligned.aa_U_2195.fasta

Run:
  python 5_tree_colorTips_overall_and_eachPhenotype_V3V10Anchors_py36.py

Outputs (transparent PNG + PDF):
  ./supp_tree_tipColor_scalar_V3V10Anchors/overall_signed_mean_lvl/tree_overall_signed_mean_lvl.png|pdf
  ./supp_tree_tipColor_scalar_V3V10Anchors/overall_signed_mean_lvl_nonzero/tree_overall_signed_mean_lvl_nonzero.png|pdf
  ./supp_tree_tipColor_scalar_V3V10Anchors/phenotypes/tree_<phenotype>.png|pdf
"""

from pathlib import Path
import datetime
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm, rcParams
from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap

import baltic as bt
from Bio import SeqIO


# ── CONFIG ────────────────────────────────────────────────────────────────
TREE_FILE = Path("timetree.nexus")

MASTER_FULL_CSV_CANDIDATES = [
    Path("master_v3_13000_all_phenotypes_with_variants_sources_with_lineage_working_rows_only.csv"),
    Path("master_v3_13000_all_phenotypes_with_variants_sources_with_lineage.csv"),
    Path("master_v3_13000_all_phenotypes_with_variants_sources.csv"),
]

FULL_AA_FASTA = Path("2344b_US_h5n1_1jan21to4may25_TR2_IQsnDNA_13000_woo_aligned.aa.fasta")
REP_FASTA     = Path("2344b_US_h5n1_1jan21to4may25_TR2_IQsnDNA_13000_woo_aligned.aa_U_2195.fasta")

# V3/V10 anchor scoring outputs:
OVERALL_CSV = Path("variant2195_V3V10Anchors_centroid_direction_overall.csv")
LONG_CSV    = Path("variant2195_V3V10Anchors_centroid_direction_long.csv")

OUT_BASE = Path("supp_tree_tipColor_scalar_V3V10Anchors")

# Visual
BASE_FONT           = 25
YEAR_LABEL_FONT     = 25
COLORBAR_TITLE_FONT = 25
COLORBAR_TICK_FONT  = 25

BRANCH_COLOUR, BRANCH_LW = "#666666", 0.65
DOT_SIZE, TIP_ALPHA      = 38, 0.85
EDGE_COLOUR, EDGE_LW     = "#444444", 0.15

MISSING_COLOUR      = (0.90, 0.90, 0.90, 0.35)
MISSING_EDGE_COLOUR = "#AAAAAA"
MISSING_EDGE_LW     = 0.12
GREY_SCALE          = 0.35

# Reference marker (only 2021 ref)
SHOW_REF_CIRCLE = True
TARGET_ACC      = "EPI_ISL_18133029"
SHOW_REF_LABEL  = False
REFERENCE_LABEL = "A/AmericanWigeon/SouthCarolina/22-000345-001/2021"

# Time axis (keep consistent with your figures)
AXIS_YEAR_START = 2021.5
AXIS_YEAR_END   = 2025.5
YEAR_LABELS     = [2022, 2023, 2024, 2025]
ENABLE_GRID     = True
GRID_LINE_WIDTH = 0.10

# Scale choices:
# For signed levels, a fixed [-4, +4] is ideal and comparable.
USE_FIXED_RANGE = True
FIXED_VMIN = -4.0
FIXED_VMAX =  4.0
COLOR_SCALE_PERCENTILES = (2, 98)  # used only if USE_FIXED_RANGE=False

# Colormap: blue ↔ white(0) ↔ red
DIVERGING_CMAP_NAME  = "RdBu_r"
ZERO_EPS = 0.0  # exact zero -> pure white

SHOW_COLORBAR = True
PNG_DPI = 300
# ──────────────────────────────────────────────────────────────────────────


# ── STYLE ─────────────────────────────────────────────────────────────────
rcParams.update({
    "font.size":       BASE_FONT,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "font.family":     "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.linewidth":  0,
})


def log(message):
    ts = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("[{0}] {1}".format(ts, message))


def resolve_path(p):
    """Try CWD first, then script folder."""
    if p.exists():
        return p
    base = Path(__file__).resolve().parent
    alt = base / p
    return alt if alt.exists() else p


def load_tree(path):
    path = resolve_path(path)
    if not path.exists():
        raise FileNotFoundError("Tree file not found: {0}".format(path))
    try:
        tree = bt.loadNexus(str(path))
    except AssertionError:
        tree = bt.loadNewick(str(path))
    return tree[0] if isinstance(tree, (list, tuple)) else tree


def is_tip(node):
    kids = getattr(node, "Children", None) or getattr(node, "children", None)
    return not bool(kids)


def label_of(node):
    return getattr(node, "Name", getattr(node, "name", ""))


def decimal_year(label):
    """
    tip label format:
      ACC|strain|state|region|host|YYYY-MM-DD
    """
    parts = label.split("|")
    if len(parts) < 2:
        return None
    date_str = parts[-1]
    d = date_str.split("-")
    if not d or not d[0].isdigit():
        return None
    y = int(d[0])
    m = int(d[1]) if len(d) > 1 and d[1].isdigit() else 1
    day = int(d[2]) if len(d) > 2 and d[2].isdigit() else 1
    from datetime import datetime as dt
    curr = dt(y, m, day)
    start, end = dt(y, 1, 1), dt(y + 1, 1, 1)
    return y + (curr - start).total_seconds() / (end - start).total_seconds()


def build_full_to_rep_map_from_master(candidates):
    for p in candidates:
        p = resolve_path(p)
        if not p.exists():
            continue
        try:
            dfm = pd.read_csv(p, usecols=["full_acc", "rep_acc"])
            dfm["full_acc"] = dfm["full_acc"].astype(str).str.upper()
            dfm["rep_acc"] = dfm["rep_acc"].astype(str).str.upper()
            mp = dict(zip(dfm["full_acc"], dfm["rep_acc"]))
            if mp:
                log("Loaded full→rep map from: {0} ({1} entries)".format(p.name, len(mp)))
                return mp
        except Exception as e:
            log("Mapping CSV unusable ({0}): {1}".format(p.name, str(e)))
            continue
    return {}


def build_full_to_rep_map_from_fasta(full_fasta, rep_fasta):
    full_fasta = resolve_path(full_fasta)
    rep_fasta = resolve_path(rep_fasta)
    if not full_fasta.exists() or not rep_fasta.exists():
        return {}
    seq_to_rep = {}
    for r in SeqIO.parse(str(rep_fasta), "fasta"):
        seq_to_rep[str(r.seq)] = r.id.split("|")[0].upper()

    full_map = {}
    n_total = 0
    n_mapped = 0
    for rec in SeqIO.parse(str(full_fasta), "fasta"):
        n_total += 1
        acc_full = rec.id.split("|")[0].upper()
        rep_acc = seq_to_rep.get(str(rec.seq))
        full_map[acc_full] = rep_acc
        if rep_acc is not None:
            n_mapped += 1
    frac = (float(n_mapped) / float(n_total)) if n_total else 0.0
    log("full→rep mapping (FASTA): {0}/{1} tips mapped ({2:.1%})".format(n_mapped, n_total, frac))
    return full_map


def make_centered_white_cmap(base_name):
    """Force PURE WHITE exactly at midpoint so 0 maps to white with TwoSlopeNorm(vcenter=0)."""
    base = cm.get_cmap(base_name, 256)
    x = np.linspace(0, 1, 256)
    colors = base(x)
    mid = 128
    colors[mid] = [1.0, 1.0, 1.0, 1.0]
    for k in (mid - 1, mid + 1):
        if 0 <= k < 256:
            colors[k] = 0.5 * colors[k] + 0.5 * np.array([1.0, 1.0, 1.0, 1.0])
    return LinearSegmentedColormap.from_list("{0}_whitecenter".format(base_name), colors)


CENTERED_WHITE_CMAP = make_centered_white_cmap(DIVERGING_CMAP_NAME)


def _fmt_tick(v):
    av = abs(float(v))
    if np.isclose(v, 0.0):
        return "0"
    if av >= 100:
        return "{0:.0f}".format(v)
    if av >= 10:
        return "{0:.1f}".format(v)
    if av >= 1:
        return "{0:.2f}".format(v)
    if av >= 0.01:
        return "{0:.3f}".format(v)
    return "{0:.2g}".format(v)


def choose_norm(values):
    if USE_FIXED_RANGE:
        vmin, vmax = float(FIXED_VMIN), float(FIXED_VMAX)
    else:
        p_lo, p_hi = COLOR_SCALE_PERCENTILES
        lo, hi = np.percentile(values, [p_lo, p_hi])
        vmin, vmax = float(min(lo, 0.0)), float(max(hi, 0.0))
        if np.isclose(vmin, vmax):
            vmin, vmax = float(np.min(values)), float(np.max(values))
            if np.isclose(vmin, vmax):
                vmin, vmax = -1.0, 1.0

    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    cmap = CENTERED_WHITE_CMAP
    ticks = [vmin, 0.0, vmax]
    ticklabels = [_fmt_tick(vmin), "0", _fmt_tick(vmax)]
    return norm, cmap, ticks, ticklabels, vmin, vmax


def sanitize_filename(s):
    s = str(s)
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def series_from_overall(overall_csv, value_col):
    overall_csv = resolve_path(overall_csv)
    if not overall_csv.exists():
        raise FileNotFoundError("Overall CSV not found: {0}".format(overall_csv))
    df = pd.read_csv(overall_csv)
    if "rep_acc" in df.columns:
        df["rep_acc"] = df["rep_acc"].astype(str).str.upper()
        df = df.set_index("rep_acc")
    elif "accession" in df.columns:
        df["accession"] = df["accession"].astype(str).str.upper()
        df = df.set_index("accession")
    else:
        raise ValueError("Overall CSV must contain 'rep_acc' or 'accession'.")

    if value_col not in df.columns:
        raise ValueError("Overall CSV missing column: {0}".format(value_col))

    ser = pd.to_numeric(df[value_col], errors="coerce")
    return ser


def series_by_phenotype(long_csv):
    long_csv = resolve_path(long_csv)
    if not long_csv.exists():
        raise FileNotFoundError("Long CSV not found: {0}".format(long_csv))
    df = pd.read_csv(long_csv)
    need = {"rep_acc", "phenotype", "signed_lvl"}
    if not need.issubset(set(df.columns)):
        raise ValueError("Long CSV must contain columns: {0}".format(", ".join(sorted(list(need)))))

    df["rep_acc"] = df["rep_acc"].astype(str).str.upper()
    df["phenotype"] = df["phenotype"].astype(str)

    out = {}
    for pheno, sub in df.groupby("phenotype"):
        ser = pd.to_numeric(sub.set_index("rep_acc")["signed_lvl"], errors="coerce")
        out[pheno] = ser
    return out


def plot_tree_scalar(tree, tips, values_by_rep, out_dir, title, full_to_rep, filename_stub):
    # map per tip
    present_vals = []
    tip_values = []
    for tip in tips:
        lbl = label_of(tip)
        acc = lbl.split("|")[0].upper() if lbl else ""
        rep = full_to_rep.get(acc)
        v = values_by_rep.get(rep, np.nan) if rep is not None else np.nan
        tip_values.append(v)
        if not pd.isna(v):
            present_vals.append(float(v))

    if len(present_vals) == 0:
        raise RuntimeError("No tip values found for: {0}".format(title))

    norm, cmap, ticks, ticklabels, vmin, vmax = choose_norm(np.array(present_vals, dtype=float))
    log("{0}: vmin={1}, vmax={2}, n_present={3}".format(title, vmin, vmax, len(present_vals)))

    # collect coordinates
    all_x = [node.x for node in tree.Objects]
    missing_x, missing_y = [], []
    present_x, present_y, colors = [], [], []

    ref_x = ref_y = None

    for tip, v in zip(tips, tip_values):
        lbl = label_of(tip)
        acc = lbl.split("|")[0].upper() if lbl else ""

        if acc == TARGET_ACC:
            ref_x, ref_y = tip.x, tip.y

        if pd.isna(v):
            missing_x.append(tip.x); missing_y.append(tip.y)
            continue

        v = float(v)
        present_x.append(tip.x); present_y.append(tip.y)

        if abs(v) <= ZERO_EPS:
            rgba = (1.0, 1.0, 1.0, TIP_ALPHA)
        else:
            rgba = list(cmap(norm(v)))
            rgba[3] = TIP_ALPHA
            rgba = tuple(rgba)
        colors.append(rgba)

    fig_height = min(max(6.0, tree.ySpan * 0.0018), 14.0)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    # branches
    tree.plotTree(ax, colour=BRANCH_COLOUR, linewidth=BRANCH_LW, alpha=0.9)

    # tips
    if missing_x:
        ax.scatter(missing_x, missing_y, s=DOT_SIZE * GREY_SCALE, c=[MISSING_COLOUR],
                   edgecolors=MISSING_EDGE_COLOUR, lw=MISSING_EDGE_LW, zorder=5)

    if present_x:
        ax.scatter(present_x, present_y, s=DOT_SIZE, c=colors,
                   edgecolors=EDGE_COLOUR, lw=EDGE_LW, zorder=7)

    # reference marker
    if SHOW_REF_CIRCLE and ref_x is not None:
        ax.scatter([ref_x], [ref_y], s=240, c="#FFC107",
                   edgecolors="black", lw=0.8, zorder=9)
        if SHOW_REF_LABEL:
            ax.text(ref_x, ref_y + 0.01 * tree.ySpan, REFERENCE_LABEL,
                    ha="center", va="bottom", fontsize=BASE_FONT, clip_on=False, zorder=10)

    # limits
    min_x, max_x = min(all_x), max(all_x)
    pad = 0.05 * (max_x - min_x) if max_x > min_x else 1.0
    ax.set_xlim(left=min_x - pad, right=max_x + pad)
    ax.margins(y=0.02)

    # time axis mapping x→year by regression on tip dates
    reg = []
    for t in tips:
        lbl = label_of(t)
        y_dec = getattr(t, "numdate", None) or decimal_year(lbl)
        if y_dec is not None:
            reg.append((t.x, y_dec))

    if reg:
        xs, ys = zip(*reg)
        m, b = np.polyfit(xs, ys, 1)

        years = YEAR_LABELS
        x_years = [(yr - b) / m for yr in years]

        y_all_tips = []
        if present_y: y_all_tips.extend(present_y)
        if missing_y: y_all_tips.extend(missing_y)

        ymin_tree = min(min(ys), min(y_all_tips))
        ymax_tree = max(max(ys), max(y_all_tips))

        x_start = (AXIS_YEAR_START - b) / m
        x_end   = (AXIS_YEAR_END - b) / m
        ybase = ymin_tree - 0.01 * tree.ySpan

        ax.hlines(ybase, x_start, x_end, lw=2.2, color="#666666", clip_on=False, zorder=12)
        long_tick = 0.0065 * tree.ySpan

        for x_pos, yr in zip(x_years, years):
            ax.vlines(x_pos, ybase, ybase - long_tick, lw=2.2, color="#666666", clip_on=False, zorder=12)
            ax.text(x_pos, ybase - 0.016 * tree.ySpan, str(yr),
                    ha="center", va="top", fontsize=YEAR_LABEL_FONT)

        if ENABLE_GRID:
            for yr in years:
                for month in (1, 4, 7, 10):
                    frac = yr + (month - 1) / 12.0
                    if AXIS_YEAR_START <= frac <= AXIS_YEAR_END:
                        xg = (frac - b) / m
                        ax.vlines(xg, ymin_tree, ymax_tree, color="#CCCCCC",
                                  lw=GRID_LINE_WIDTH, ls=(0, (2, 2)),
                                  clip_on=False, zorder=1)

    # cosmetics
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_visible(True)

    # colorbar
    if SHOW_COLORBAR:
        pos = ax.get_position()
        cax = fig.add_axes([pos.x0 + (pos.width - 0.32) / 2.0, pos.y1 + 0.085, 0.32, 0.018])

        sm = cm.ScalarMappable(norm=norm, cmap=cmap)
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_ticks(ticks)
        cb.set_ticklabels(ticklabels)
        cb.set_label(title, fontsize=COLORBAR_TITLE_FONT, labelpad=4)

        cb.ax.xaxis.set_label_position("top")
        cb.ax.xaxis.set_ticks_position("bottom")
        cb.ax.tick_params(labelsize=COLORBAR_TICK_FONT, pad=1.5, length=3,
                          bottom=True, labelbottom=True, top=False, labeltop=False)

    # save
    out_dir.mkdir(exist_ok=True, parents=True)
    for ext in (".pdf", ".png"):
        fn = out_dir / "{0}{1}".format(filename_stub, ext)
        fig.savefig(str(fn),
                    dpi=(PNG_DPI if ext == ".png" else None),
                    bbox_inches="tight", transparent=True, pad_inches=1)

    plt.close(fig)
    log("Saved: {0}".format(out_dir / (filename_stub + ".[pdf|png]")))


def main():
    # mapping
    full_to_rep = build_full_to_rep_map_from_master(MASTER_FULL_CSV_CANDIDATES)
    if not full_to_rep:
        log("No usable mapping CSV found; falling back to FASTA mapping …")
        full_to_rep = build_full_to_rep_map_from_fasta(FULL_AA_FASTA, REP_FASTA)
    if not full_to_rep:
        raise RuntimeError("Could not build full→rep mapping (need mapping CSV or FASTA mapping files).")

    # tree
    tree = load_tree(TREE_FILE)
    tips = [n for n in tree.Objects if is_tip(n)]
    log("Tree tips: {0}".format(len(tips)))

    # ── (A) Overall trees (both) ───────────────────────────────────────────
    overall_cols = [
        ("overall_signed_mean_lvl", "Overall lineage-likeness (V3 avian ↔ V10 cattle), mean incl. zeros"),
        ("overall_signed_mean_lvl_nonzero", "Overall lineage-likeness (V3 avian ↔ V10 cattle), mean nonzero only"),
    ]
    for col, title in overall_cols:
        ser = series_from_overall(OVERALL_CSV, col)
        subdir = OUT_BASE / col
        plot_tree_scalar(
            tree=tree,
            tips=tips,
            values_by_rep=ser,
            out_dir=subdir,
            title=title,
            full_to_rep=full_to_rep,
            filename_stub="tree_{0}".format(col),
        )

    # ── (B) Per-phenotype signed_lvl trees ─────────────────────────────────
    pheno_map = series_by_phenotype(LONG_CSV)
    log("Phenotypes found: {0}".format(len(pheno_map)))

    pheno_dir = OUT_BASE / "phenotypes"
    for pheno, ser in sorted(pheno_map.items(), key=lambda x: x[0].lower()):
        fname = "tree_{0}".format(sanitize_filename(pheno))
        title = "{0} (signed level: V3 avian ↔ V10 cattle)".format(pheno)
        plot_tree_scalar(
            tree=tree,
            tips=tips,
            values_by_rep=ser,
            out_dir=pheno_dir,
            title=title,
            full_to_rep=full_to_rep,
            filename_stub=fname,
        )

    log("DONE.")


if __name__ == "__main__":
    main()
