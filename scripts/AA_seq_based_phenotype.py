#!/usr/bin/env python3
"""
ha2344b_metrics_clean_manuscript_only.py

Sequence-derived HA phenotype panel restricted to the metrics described in the
manuscript.

Included metrics
----------------
Full HA sequence:
  - pI
  - net charge at pH 7
  - GRAVY
  - aliphatic index
  - fraction of positively charged residues
  - fraction of negatively charged residues
  - secondary-structure propensities: helix_frac, sheet_frac, turn_frac

HA1 and HA2 domains:
  - pI
  - net charge at pH 7
  - GRAVY
  - aliphatic index
  - fraction of positively charged residues
  - fraction of negatively charged residues

Notes
-----
- Metrics are computed on gap-stripped amino-acid sequences.
- HA1/HA2 domain boundaries are read from site_numbering_map.csv.
- A delta table relative to the chosen reference sequence is also written,
  because reference-relative phenotype shifts are used downstream in the
  manuscript figures.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ─── USER SETTINGS ────────────────────────────────────────────────────────
FASTA_PATH = "2344b_US_h5n1_1jan21to4may25_TR2_IQsnDNA_13000_woo_aligned.aa_U_2195.fasta"
REF_ID_SUBSTR = "EPI_ISL_18133029"
SITE_MAP_CSV = "site_numbering_map.csv"
OUT_DIR = "aa_metrics_out"
# ──────────────────────────────────────────────────────────────────────────

STD_AA = set("ACDEFGHIKLMNPQRSTVWY")
OUT_PATH = Path(OUT_DIR)
OUT_PATH.mkdir(parents=True, exist_ok=True)


# Fallback for older Biopython versions without ProteinAnalysis.aliphatic_index
if not hasattr(ProteinAnalysis, "aliphatic_index"):
    def _aliphatic_index(seq: str) -> float:
        counts = Counter(seq)
        n = float(len(seq)) if seq else 1.0
        return (counts["A"] / n + 2.9 * counts["V"] / n + 3.9 * (counts["I"] + counts["L"]) / n) * 100.0

    ProteinAnalysis.aliphatic_index = lambda self: _aliphatic_index(self.sequence)  # type: ignore[attr-defined]


def load_site_map(path: str | Path) -> tuple[list[int], list[int]]:
    """Return HA1 and HA2 sequential-site lists from the site map."""
    df = pd.read_csv(path, dtype=str)
    df["sequential_site"] = pd.to_numeric(df["sequential_site"], errors="coerce")

    ha1_sites = sorted(
        df.loc[df["region"] == "HA1", "sequential_site"]
        .dropna()
        .astype(int)
        .tolist()
    )
    ha2_sites = sorted(
        df.loc[df["region"] == "HA2", "sequential_site"]
        .dropna()
        .astype(int)
        .tolist()
    )

    if not ha1_sites or not ha2_sites:
        raise ValueError("Failed to recover HA1/HA2 sites from site_numbering_map.csv")

    return ha1_sites, ha2_sites


def seq_to_aln_map(ref_aln_seq: str) -> dict[int, int]:
    """Map ungapped reference position (1-based) to alignment index (0-based)."""
    mapping: dict[int, int] = {}
    ungapped_pos = 0
    for aln_idx, aa in enumerate(ref_aln_seq):
        if aa != "-":
            ungapped_pos += 1
            mapping[ungapped_pos] = aln_idx
    return mapping


def clean_aa(seq: str) -> str:
    """Keep only standard amino acids."""
    return "".join(aa for aa in seq if aa in STD_AA)


def extract_region_from_alignment(aln_seq: str, seq2aln: dict[int, int], seq_sites: list[int]) -> str:
    """Extract a domain sequence using sequential-site coordinates."""
    residues: list[str] = []
    for site in seq_sites:
        aln_idx = seq2aln.get(site)
        if aln_idx is None:
            continue
        aa = aln_seq[aln_idx]
        if aa in STD_AA:
            residues.append(aa)
    return "".join(residues)


def basic_composition_metrics(seq: str) -> dict[str, float]:
    """Return manuscript-used ProtParam metrics for one amino-acid sequence."""
    if not seq:
        return {
            "pI": np.nan,
            "net_charge_pH7": np.nan,
            "gravy": np.nan,
            "aliphatic_index": np.nan,
            "pos_frac": np.nan,
            "neg_frac": np.nan,
        }

    pa = ProteinAnalysis(seq)
    n = float(len(seq))
    return {
        "pI": round(pa.isoelectric_point(), 3),
        "net_charge_pH7": round(pa.charge_at_pH(7.0), 3),
        "gravy": round(pa.gravy(), 3),
        "aliphatic_index": round(pa.aliphatic_index(), 3),
        "pos_frac": round((seq.count("K") + seq.count("R")) / n, 4),
        "neg_frac": round((seq.count("D") + seq.count("E")) / n, 4),
    }



def secondary_structure_metrics(seq: str) -> dict[str, float]:
    """Return ProtParam secondary-structure propensity fractions for full HA."""
    if not seq:
        return {
            "helix_frac": np.nan,
            "sheet_frac": np.nan,
            "turn_frac": np.nan,
        }

    pa = ProteinAnalysis(seq)
    helix, turn, sheet = pa.secondary_structure_fraction()
    return {
        "helix_frac": round(helix, 3),
        "sheet_frac": round(sheet, 3),
        "turn_frac": round(turn, 3),
    }



def prefix_keys(d: dict[str, float], prefix: str) -> dict[str, float]:
    return {f"{key}_{prefix}": value for key, value in d.items()}



def calc_metrics(aln_seq: str, seq2aln: dict[int, int], ha1_sites: list[int], ha2_sites: list[int]) -> dict[str, float]:
    """Compute manuscript-only sequence phenotypes for one aligned HA sequence."""
    full_seq = clean_aa(aln_seq.replace("-", ""))
    ha1_seq = extract_region_from_alignment(aln_seq, seq2aln, ha1_sites)
    ha2_seq = extract_region_from_alignment(aln_seq, seq2aln, ha2_sites)

    metrics: dict[str, float] = {}

    # Full HA metrics
    metrics.update(prefix_keys(basic_composition_metrics(full_seq), "HA"))
    metrics.update(secondary_structure_metrics(full_seq))

    # Domain metrics
    metrics.update(prefix_keys(basic_composition_metrics(ha1_seq), "HA1"))
    metrics.update(prefix_keys(basic_composition_metrics(ha2_seq), "HA2"))

    return metrics



def ordered_columns() -> list[str]:
    return [
        "pI_HA",
        "pI_HA1",
        "pI_HA2",
        "net_charge_pH7_HA",
        "net_charge_pH7_HA1",
        "net_charge_pH7_HA2",
        "gravy_HA",
        "gravy_HA1",
        "gravy_HA2",
        "aliphatic_index_HA",
        "aliphatic_index_HA1",
        "aliphatic_index_HA2",
        "pos_frac_HA",
        "pos_frac_HA1",
        "pos_frac_HA2",
        "neg_frac_HA",
        "neg_frac_HA1",
        "neg_frac_HA2",
        "helix_frac",
        "sheet_frac",
        "turn_frac",
    ]



def main() -> None:
    ha1_sites, ha2_sites = load_site_map(SITE_MAP_CSV)

    records = list(SeqIO.parse(FASTA_PATH, "fasta"))
    if not records:
        raise SystemExit("[ERR] No sequences found in FASTA")

    try:
        ref_record = next(rec for rec in records if REF_ID_SUBSTR in rec.id)
    except StopIteration as exc:
        raise SystemExit(f"[ERR] No reference sequence with id containing '{REF_ID_SUBSTR}'") from exc

    ref_aln_seq = str(ref_record.seq)
    seq2aln = seq_to_aln_map(ref_aln_seq)

    rows: list[dict[str, float | str]] = []
    for rec in records:
        aln_seq = str(rec.seq)
        row = calc_metrics(aln_seq, seq2aln, ha1_sites, ha2_sites)
        row["seq_id"] = rec.id
        rows.append(row)

    df = pd.DataFrame(rows).set_index("seq_id")
    cols = [c for c in ordered_columns() if c in df.columns]
    df = df[cols]

    metrics_csv = OUT_PATH / "ha_metrics_manuscript_only.csv"
    df.to_csv(metrics_csv)
    print(f"[✓] wrote metrics -> {metrics_csv}")

    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    ref_vals = df.loc[ref_record.id, numeric_cols]
    delta_df = df[numeric_cols].subtract(ref_vals, axis=1)

    delta_csv = OUT_PATH / "ha_metrics_manuscript_only_delta.csv"
    delta_df.to_csv(delta_csv)
    print(f"[✓] wrote delta metrics -> {delta_csv}")


if __name__ == "__main__":
    main()
