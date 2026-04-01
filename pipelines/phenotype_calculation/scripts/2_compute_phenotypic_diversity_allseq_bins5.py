#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diversity_allSeq_SEQ_bins5_individualPlusCombined_groupedColors_pngTransparent_noTitles.py

What it does:
- Uses ALL sequences (no host types, no lineages, no subsets)
- SEQ basis only
- Exactly 5 bins (quantile bins)
- Computes effective bins (Shannon effective number of bins) per time period
- Makes:
    1) Individual figure for each phenotype (metric)
    2) Combined all-phenotypes figure + legend-only PNG
- No titles at the top of figures
- PNG only, all transparent
- NEW output folder only (does not touch old results)

Python 3.6+ compatible.
"""

import math
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator


# =========================
# USER SETTINGS
# =========================
MASTER_CSV = "master_v3_13000_all_phenotypes_with_variants_sources.csv"

OUTDIR = "rep_figures_SEQ_bins5_allSeq_individualPlusCombined_groupedPalette_v1"

VIEW = "quarterly"   # "quarterly" / "monthly" / "yearly"
BINS = 5             # fixed

OUTBREAK_DATE = "2024-03-25"  # vertical dashed line (set None to disable)

# sample-size rules
MIN_PERIOD_N = 25
MIN_COUNT_PER_BIN = 5

# bootstrap / rarefaction
DO_BOOTSTRAP = True
BOOTSTRAP_ITERS = 600
SEED = 7

# Outputs
MAKE_COMBINED = True
MAKE_INDIVIDUAL = True
WRITE_COLOR_KEY_CSV = True

# Appearance
DPI = 600
FIG_TRANSPARENT = True

FIGSIZE_SINGLE = (18, 6.2)
FIGSIZE_COMBINED = (18, 6.2)
FIGSIZE_LEGEND = (10, 7)

# sample-size bar styling
N_BAR_FRAC = 0.028
N_BAR_ALPHA = 0.10
N_BAR_COLOR = "0.35"

# low-sample shade styling
LOW_SAMPLE_SHADE_COLOR = "#94a3b8"
LOW_SAMPLE_SHADE_ALPHA = 0.18
LOW_SAMPLE_HATCH = ".."
LOW_SAMPLE_EDGE_COLOR = "#64748b"

# raw low-sample points styling
RAW_DOT_SIZE = 60
RAW_DOT_LW = 1.6

# Typography (kept similar to your pipeline vibe)
TITLE_FONTSIZE = 46   # used for year labels below axis, NOT top titles
TICK_FONTSIZE  = 36
LABEL_FONTSIZE = 40
SMALL_FONTSIZE = 18

# draw faint raw line (usually False)
RAW_FALLBACK = False


# =========================
# STYLE
# =========================
def set_style():
    plt.rcParams.update({
        "figure.dpi": DPI,
        "savefig.dpi": DPI,
        "figure.facecolor": "none" if FIG_TRANSPARENT else "white",
        "axes.facecolor": "none" if FIG_TRANSPARENT else "white",
        "axes.linewidth": 1.0,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.width": 0.9,
        "ytick.major.width": 0.9,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    })


def save_png(fig, out_png):
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_png), dpi=DPI, transparent=FIG_TRANSPARENT,
                bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


# =========================
# TIME HELPERS
# =========================
def period_label(dt, view):
    if view == "quarterly":
        q = (dt.month - 1) // 3 + 1
        return "%dQ%d" % (dt.year, q)
    if view == "monthly":
        return dt.strftime("%Y-%m")
    if view == "yearly":
        return str(dt.year)
    raise ValueError("Unknown view: %s" % view)


def period_start_date(period, view):
    if view == "quarterly":
        y = int(period[:4])
        q = int(period[-1])
        m = 1 + (q - 1) * 3
        return datetime(y, m, 1)
    if view == "monthly":
        y, m = period.split("-")
        return datetime(int(y), int(m), 1)
    if view == "yearly":
        return datetime(int(period), 1, 1)
    raise ValueError(view)


def sort_periods(periods, view):
    return sorted(periods, key=lambda p: period_start_date(p, view))


def compute_outbreak_x(periods, view, outbreak_dt):
    if outbreak_dt is None or not periods:
        return None
    starts = [period_start_date(p, view) for p in periods]
    start_nums = np.array([s.toordinal() for s in starts], dtype=float)
    out_num = float(outbreak_dt.toordinal())

    if out_num <= start_nums[0]:
        return 0.0
    if out_num >= start_nums[-1]:
        return float(len(periods) - 1)

    i = int(np.searchsorted(start_nums, out_num, side="right") - 1)
    i = max(0, min(i, len(start_nums) - 2))
    frac = (out_num - start_nums[i]) / (start_nums[i + 1] - start_nums[i])
    return float(i + frac)


def year_boundaries(periods):
    years = [int(p[:4]) for p in periods]
    blocks = []
    if not years:
        return blocks
    cur = years[0]
    s = 0
    for i in range(1, len(years)):
        if years[i] != cur:
            blocks.append((cur, s, i - 1))
            cur = years[i]
            s = i
    blocks.append((cur, s, len(periods) - 1))
    return blocks


# =========================
# DIVERSITY + BINNING
# =========================
def quantile_edges(series, k):
    s = pd.to_numeric(series, errors="coerce").dropna().astype(float)
    if s.empty:
        return np.array([0.0, 1.0], dtype=float)

    uniq = int(s.nunique())
    k_use = max(1, min(int(k), uniq))

    try:
        _, edges = pd.qcut(s, q=k_use, retbins=True, duplicates="drop")
        edges = np.unique(edges)
        if len(edges) < 2:
            edges = np.array([float(s.min()), float(s.max())], dtype=float)
    except Exception:
        edges = np.array([float(s.min()), float(s.max())], dtype=float)

    if len(edges) == 2 and edges[0] == edges[1]:
        edges = np.array([edges[0] - 1.0, edges[1] + 1.0], dtype=float)
    return edges.astype(float)


def assign_bins(values, edges):
    vals = pd.to_numeric(values, errors="coerce").astype(float).values
    out = np.full(vals.shape[0], -1, dtype=int)
    good = np.isfinite(vals)
    if not good.any():
        return out
    k = len(edges) - 1
    idx = np.searchsorted(edges, vals[good], side="right") - 1
    idx = np.clip(idx, 0, max(0, k - 1))
    out[good] = idx.astype(int)
    return out


def shannon_effective_bins(counts):
    total = float(np.sum(counts))
    if total <= 0:
        return np.nan
    p = counts[counts > 0].astype(float) / total
    H = -np.sum(p * np.log2(p))
    return float(2.0 ** H)


def supported_bins(counts, min_count):
    return int(np.sum(counts >= int(min_count)))


def bootstrap_rarefy(bin_idx, k, rarefy_n, iters, rng, min_count):
    n = int(len(bin_idx))
    if rarefy_n <= 0 or n < rarefy_n:
        return (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    sh = np.empty(iters, dtype=float)
    sup = np.empty(iters, dtype=float)

    for i in range(iters):
        pick = rng.choice(n, size=rarefy_n, replace=False)
        counts = np.bincount(bin_idx[pick], minlength=k)
        sh[i] = shannon_effective_bins(counts)
        sup[i] = float(supported_bins(counts, min_count))

    def ci(arr):
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            return (np.nan, np.nan, np.nan)
        return (float(np.mean(arr)),
                float(np.percentile(arr, 2.5)),
                float(np.percentile(arr, 97.5)))

    sh_m, sh_lo, sh_hi = ci(sh)
    su_m, su_lo, su_hi = ci(sup)
    return (sh_m, sh_lo, sh_hi, su_m, su_lo, su_hi)


def compute_rarefy_n(n_by_period, min_period_n):
    eligible = [int(v) for v in n_by_period.values() if int(v) >= int(min_period_n)]
    return int(min(eligible)) if eligible else 0


# =========================
# PLOT DECORATIONS
# =========================
def add_sample_size_bars(ax, x, nvals):
    if nvals is None or len(nvals) == 0:
        return
    nvals = np.asarray(nvals, dtype=float)
    finite = np.isfinite(nvals)
    if not finite.any():
        return
    max_n = float(np.nanmax(nvals))
    if max_n <= 0:
        return

    heights = (nvals / max_n) * float(N_BAR_FRAC)
    ax.bar(
        x, heights, width=0.90, bottom=0.0,
        color=N_BAR_COLOR, alpha=float(N_BAR_ALPHA), linewidth=0.0,
        transform=ax.get_xaxis_transform(),
        zorder=0.70, clip_on=False
    )

    fs = int(SMALL_FONTSIZE)
    yoff = float(N_BAR_FRAC) * 0.06
    for i, h in enumerate(heights):
        if not np.isfinite(h):
            continue
        ax.text(
            x[i], h + yoff, "%d" % int(nvals[i]),
            transform=ax.get_xaxis_transform(),
            ha="center", va="bottom",
            fontsize=fs, color="0.35",
            alpha=0.92, zorder=0.75, clip_on=False
        )


def shade_low_sample_periods(ax, low_mask):
    if low_mask is None:
        return
    for i, is_low in enumerate(low_mask):
        if bool(is_low):
            ax.axvspan(
                i - 0.5, i + 0.5,
                facecolor=LOW_SAMPLE_SHADE_COLOR,
                alpha=float(LOW_SAMPLE_SHADE_ALPHA),
                hatch=LOW_SAMPLE_HATCH,
                edgecolor=LOW_SAMPLE_EDGE_COLOR,
                linewidth=0.0,
                zorder=0.50
            )


def add_quarter_separators(ax, n):
    for i in range(1, n):
        ax.axvline(i - 0.5, color="0.65", lw=0.6, ls=(0, (2, 3)), alpha=0.35, zorder=0)


def add_year_lines_and_labels(ax, periods):
    blocks = year_boundaries(periods)
    for (y, s, e) in blocks:
        if s > 0:
            ax.axvline(s - 0.5, color="0.30", lw=1.0, alpha=0.55, zorder=1)
    for (y, s, e) in blocks:
        xc = (s + e) / 2.0
        ax.text(
            xc, -0.14, str(y),
            transform=ax.get_xaxis_transform(),
            ha="center", va="top",
            fontsize=TITLE_FONTSIZE, color="0.12",
            clip_on=False
        )


def add_outbreak_line(ax, outbreak_x):
    if outbreak_x is None:
        return
    ax.axvline(outbreak_x, color="#0ea5a4", lw=2.0, ls=(0, (6, 6)), alpha=0.95, zorder=5)


# =========================
# GROUPED COLOR PALETTE
# =========================
def metric_group(metric):
    m = str(metric).lower()

    if "flexddg" in m or "bind" in m:
        return "binding"
    if "cartddg" in m or m.startswith("ph_") or "foldx" in m or "stabil" in m:
        return "stability"
    if "af3_rmsd" in m or ("rmsd" in m and "af3" in m):
        return "structure"
    if m.startswith("af3_d") or ("rbs" in m and "d" in m):
        return "geometry"
    if "helix_frac" in m or "sheet_frac" in m or "turn_frac" in m:
        return "secstruct"
    if ("net_charge" in m or "gravy" in m or "pos_frac" in m or "neg_frac" in m or
        "pi_" in m or m == "pi" or "aliph" in m or "aromat" in m):
        return "composition"
    return "other"


GROUP_BASE = {
    "composition": "#2ca02c",
    "secstruct":   "#17becf",
    "structure":   "#1f77b4",
    "geometry":    "#9467bd",
    "stability":   "#d62728",
    "binding":     "#ff7f0e",
    "other":       "#7f7f7f",
}


def blend_to_white(rgb, frac):
    return tuple((1.0 - frac) * c + frac * 1.0 for c in rgb)


def blend_to_black(rgb, frac):
    return tuple((1.0 - frac) * c + frac * 0.0 for c in rgb)


def make_group_shades(base_hex, n):
    base = mcolors.to_rgb(base_hex)
    if n <= 1:
        return [mcolors.to_hex(base)]
    span = 0.28
    ws = np.linspace(-span, span, n)
    out = []
    for w in ws:
        if w < 0:
            c = blend_to_black(base, abs(w))
        else:
            c = blend_to_white(base, w)
        out.append(mcolors.to_hex(c))
    return out


def palette_for_metrics_grouped(metrics):
    groups = {}
    for m in metrics:
        g = metric_group(m)
        groups.setdefault(g, []).append(m)

    for g in groups:
        groups[g] = sorted(groups[g])

    colors = {}
    for g, ms in groups.items():
        base = GROUP_BASE.get(g, GROUP_BASE["other"])
        shades = make_group_shades(base, len(ms))
        for m, col in zip(ms, shades):
            colors[m] = col

    for m in metrics:
        if m not in colors:
            colors[m] = GROUP_BASE["other"]
    return colors


# =========================
# METRICS (auto-detect)
# =========================
def list_numeric_metrics(df):
    meta = set([
        "full_acc", "strain", "state", "region", "host", "date",
        "rep_acc", "variant_num", "variant_count", "period", "lineage"
    ])
    metrics = []
    for c in df.columns:
        if c in meta:
            continue
        if str(c).startswith("source_"):
            continue
        if str(c).endswith("_std"):
            continue
        if pd.api.types.is_numeric_dtype(df[c]):
            metrics.append(c)
    return sorted(metrics)


# =========================
# SUMMARY BUILD (ALL SEQ)
# =========================
def build_summary_allseq(df, metrics, edges_by_metric, view,
                         min_period_n, min_count_per_bin,
                         do_bootstrap, iters, seed):

    periods = sort_periods(df["period"].dropna().unique().tolist(), view)
    if not periods:
        return pd.DataFrame(), [], {}

    n_by_period = {p: int((df["period"] == p).sum()) for p in periods}
    rarefy_n = compute_rarefy_n(n_by_period, min_period_n)

    rng = np.random.RandomState(int(seed))
    rows = []

    for metric in metrics:
        edges = edges_by_metric[metric]
        k = len(edges) - 1
        bin_idx_all = assign_bins(df[metric], edges)

        for p in periods:
            mask = (df["period"].values == p)
            idx = bin_idx_all[mask]
            idx = idx[idx >= 0]
            n_items = int(mask.sum())
            counts = np.bincount(idx, minlength=k) if idx.size else np.zeros(k, dtype=int)

            raw_sh = shannon_effective_bins(counts)
            raw_sup = float(supported_bins(counts, min_count_per_bin))

            eligible = (do_bootstrap and rarefy_n > 0 and n_items >= max(min_period_n, rarefy_n) and idx.size >= rarefy_n)
            if eligible:
                sh_m, sh_lo, sh_hi, su_m, su_lo, su_hi = bootstrap_rarefy(
                    idx.astype(int), k, rarefy_n, iters, rng, min_count_per_bin
                )
            else:
                sh_m = sh_lo = sh_hi = np.nan
                su_m = su_lo = su_hi = np.nan

            rows.append({
                "view": view,
                "basis": "seq",
                "period": p,
                "metric": metric,
                "bins": int(BINS),
                "k_bins": int(k),
                "n_items": int(n_items),
                "rarefy_n": int(rarefy_n),
                "eligible": bool(eligible),
                "raw_shannon_q1": raw_sh,
                "raw_supported_bins": raw_sup,
                "mean_shannon": sh_m,
                "lo_shannon": sh_lo,
                "hi_shannon": sh_hi,
                "mean_supported": su_m,
                "lo_supported": su_lo,
                "hi_supported": su_hi,
            })

    return pd.DataFrame(rows), periods, n_by_period


# =========================
# PLOTS
# =========================
def stylize_axes(ax):
    ax.tick_params(axis="y", labelsize=TICK_FONTSIZE)
    ax.set_yticks(np.arange(1, BINS + 1))
    ax.set_ylim(0.9, BINS + 0.2)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.grid(True, linewidth=0.9, alpha=0.12)
    ax.set_axisbelow(True)
    ax.margins(y=0.12)


def plot_single_metric(summary, periods, n_by_period, outbreak_x, metric, color_hex, out_png):
    x = np.arange(len(periods))

    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE)
    fig.patch.set_facecolor("none" if FIG_TRANSPARENT else "white")
    ax.set_facecolor("none" if FIG_TRANSPARENT else "white")

    ax.set_ylabel("Effective Bins", fontsize=LABEL_FONTSIZE)

    nvals = np.array([n_by_period.get(p, 0) for p in periods], dtype=int)
    add_sample_size_bars(ax, x, nvals)
    low_mask = np.asarray(nvals, dtype=float) < float(MIN_PERIOD_N)
    shade_low_sample_periods(ax, low_mask)

    add_quarter_separators(ax, len(periods))
    add_year_lines_and_labels(ax, periods)
    add_outbreak_line(ax, outbreak_x)

    sub = summary[summary["metric"] == metric].set_index("period").reindex(periods)

    y = pd.to_numeric(sub["mean_shannon"], errors="coerce").values.astype(float)
    lo = pd.to_numeric(sub["lo_shannon"], errors="coerce").values.astype(float)
    hi = pd.to_numeric(sub["hi_shannon"], errors="coerce").values.astype(float)
    raw = pd.to_numeric(sub["raw_shannon_q1"], errors="coerce").values.astype(float)

    c = color_hex

    if RAW_FALLBACK:
        ax.plot(x, raw, color=c, lw=2.0, alpha=0.18, zorder=2)

    good = np.isfinite(y) & np.isfinite(lo) & np.isfinite(hi)
    if np.any(good):
        ax.fill_between(x[good], lo[good], hi[good], color=c, alpha=0.12, linewidth=0, zorder=3)
    ax.plot(x, y, color=c, lw=4.0, alpha=0.94, zorder=4)

    raw_mask = low_mask & np.isfinite(raw) & (~np.isfinite(y))
    if np.any(raw_mask):
        ax.scatter(
            x[raw_mask], raw[raw_mask],
            s=RAW_DOT_SIZE,
            facecolors="none",
            edgecolors=c,
            linewidths=RAW_DOT_LW,
            alpha=0.98,
            zorder=5,
        )

    # no top title by request
    ax.set_xticks(x)
    ax.set_xticklabels([""] * len(periods))

    stylize_axes(ax)
    save_png(fig, out_png)


def plot_combined(summary, periods, metrics, n_by_period, outbreak_x, colors,
                  out_png, legend_out_png):
    x = np.arange(len(periods))

    fig, ax = plt.subplots(figsize=FIGSIZE_COMBINED)
    fig.patch.set_facecolor("none" if FIG_TRANSPARENT else "white")
    ax.set_facecolor("none" if FIG_TRANSPARENT else "white")

    ax.set_ylabel("Effective Bins", fontsize=LABEL_FONTSIZE)

    nvals = np.array([n_by_period.get(p, 0) for p in periods], dtype=int)
    add_sample_size_bars(ax, x, nvals)
    low_mask = np.asarray(nvals, dtype=float) < float(MIN_PERIOD_N)
    shade_low_sample_periods(ax, low_mask)

    add_quarter_separators(ax, len(periods))
    add_year_lines_and_labels(ax, periods)
    add_outbreak_line(ax, outbreak_x)

    mean_lines = []

    for m in metrics:
        sub = summary[summary["metric"] == m].set_index("period").reindex(periods)

        y = pd.to_numeric(sub["mean_shannon"], errors="coerce").values.astype(float)
        lo = pd.to_numeric(sub["lo_shannon"], errors="coerce").values.astype(float)
        hi = pd.to_numeric(sub["hi_shannon"], errors="coerce").values.astype(float)
        raw = pd.to_numeric(sub["raw_shannon_q1"], errors="coerce").values.astype(float)

        c = colors[m]

        if RAW_FALLBACK:
            ax.plot(x, raw, color=c, lw=2.0, alpha=0.14, zorder=2)

        good = np.isfinite(y) & np.isfinite(lo) & np.isfinite(hi)
        if np.any(good):
            ax.fill_between(x[good], lo[good], hi[good], color=c, alpha=0.08, linewidth=0, zorder=3)
        ax.plot(x, y, color=c, lw=3.0, alpha=0.90, zorder=4)

        raw_mask = low_mask & np.isfinite(raw) & (~np.isfinite(y))
        if np.any(raw_mask):
            ax.scatter(
                x[raw_mask], raw[raw_mask],
                s=RAW_DOT_SIZE,
                facecolors="none",
                edgecolors=c,
                linewidths=RAW_DOT_LW,
                alpha=0.98,
                zorder=5,
            )

        mean_lines.append(y)

    # mean across metrics (dashed black)
    if mean_lines:
        mat = np.vstack([np.array(v, dtype=float) for v in mean_lines])
        finite = np.isfinite(mat)
        col_n = finite.sum(axis=0)
        col_sum = np.where(finite, mat, 0.0).sum(axis=0)

        ymean = np.full(mat.shape[1], np.nan)
        ok = col_n >= max(2, int(math.ceil(len(metrics) * 0.4)))
        ymean[ok] = col_sum[ok] / col_n[ok]
        ax.plot(x, ymean, color="0.10", lw=5.0, ls=(0, (2, 6)), alpha=0.95, zorder=6)

    # no top title by request
    ax.set_xticks(x)
    ax.set_xticklabels([""] * len(periods))

    stylize_axes(ax)
    save_png(fig, out_png)

    # legend-only PNG (transparent)
    handles = []
    for m in metrics:
        handles.append(Line2D([0], [0], color=colors[m], lw=5.0, label=m))
    handles.append(Line2D([0], [0], color="0.10", lw=6.0, ls=(0, (2, 6)), label="Mean across metrics"))
    handles.append(Patch(facecolor="0.25", alpha=0.10, edgecolor="none", label="95% CI (shaded)"))
    handles.append(Line2D([0], [0], marker="o", markersize=10, linestyle="None",
                          markerfacecolor="none", markeredgewidth=2.0,
                          markeredgecolor="0.25", label="Raw (low-sample)"))

    fig2, ax2 = plt.subplots(figsize=FIGSIZE_LEGEND)
    fig2.patch.set_facecolor("none" if FIG_TRANSPARENT else "white")
    ax2.set_facecolor("none" if FIG_TRANSPARENT else "white")
    ax2.axis("off")
    ax2.legend(handles=handles, frameon=False, loc="center", fontsize=26, ncol=1)
    save_png(fig2, legend_out_png)


# =========================
# MAIN
# =========================
def main():
    set_style()

    master_fp = Path(MASTER_CSV)
    if not master_fp.exists():
        raise SystemExit("ERROR: master csv not found: %s" % str(master_fp))

    df = pd.read_csv(str(master_fp))

    # Date -> period
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df[df["date"].notna()].copy()
    df["period"] = df["date"].apply(lambda d: period_label(d.to_pydatetime(), VIEW))

    # Metrics
    metrics = list_numeric_metrics(df)
    if not metrics:
        raise SystemExit("ERROR: no numeric phenotype columns detected.")

    # Bin edges (IMPORTANT: based on ALL sequences, as you requested)
    edges_by_metric = {}
    for m in metrics:
        edges_by_metric[m] = quantile_edges(df[m], BINS)

    # Grouped palette (unique per metric, families by type)
    colors = palette_for_metrics_grouped(metrics)

    # Outbreak marker
    outbreak_dt = None
    if OUTBREAK_DATE:
        outbreak_dt = datetime.strptime(OUTBREAK_DATE, "%Y-%m-%d")

    # Build summary
    summary, periods, n_by_period = build_summary_allseq(
        df, metrics, edges_by_metric,
        view=VIEW,
        min_period_n=MIN_PERIOD_N,
        min_count_per_bin=MIN_COUNT_PER_BIN,
        do_bootstrap=DO_BOOTSTRAP,
        iters=BOOTSTRAP_ITERS,
        seed=SEED
    )

    if summary.empty or not periods:
        raise SystemExit("ERROR: no valid periods after filtering.")

    out_base = Path(OUTDIR) / VIEW / "allSeq" / "seq" / ("bins%d" % int(BINS))
    fig_dir_single = out_base / "figures_individual"
    fig_dir_combo = out_base / "figures_combined"
    out_base.mkdir(parents=True, exist_ok=True)

    summary.to_csv(str(out_base / "Summary.csv"), index=False)

    if WRITE_COLOR_KEY_CSV:
        pd.DataFrame({
            "metric": metrics,
            "group": [metric_group(m) for m in metrics],
            "color_hex": [colors[m] for m in metrics],
        }).to_csv(str(out_base / "MetricColors.csv"), index=False)

    outbreak_x = compute_outbreak_x(periods, VIEW, outbreak_dt) if outbreak_dt else None

    # Individual figures
    if MAKE_INDIVIDUAL:
        for m in metrics:
            safe = "".join([ch if ch.isalnum() or ch in ("_", "-", ".") else "_" for ch in str(m)])
            out_png = fig_dir_single / ("%s_SEQ_bins5.png" % safe)
            plot_single_metric(
                summary=summary,
                periods=periods,
                n_by_period=n_by_period,
                outbreak_x=outbreak_x,
                metric=m,
                color_hex=colors[m],
                out_png=out_png
            )

    # Combined
    if MAKE_COMBINED:
        plot_combined(
            summary=summary,
            periods=periods,
            metrics=metrics,
            n_by_period=n_by_period,
            outbreak_x=outbreak_x,
            colors=colors,
            out_png=fig_dir_combo / "AllPhenotypes_Shannon_SEQ_bins5.png",
            legend_out_png=fig_dir_combo / "AllPhenotypes_Legend.png"
        )

    print("[Done] Outputs written to:", str(out_base))


if __name__ == "__main__":
    main()
