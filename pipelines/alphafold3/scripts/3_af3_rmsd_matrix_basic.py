#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path
from typing import Dict, Tuple, Optional

import numpy as np
import pandas as pd
import MDAnalysis as mda

RBS_BLOCKS = [
    ("130-loop", 143, 150),
    ("150-loop", 162, 172),
    ("190-helix", 199, 209),
    ("220-loop", 232, 240),
]
RBS_RESIDS = {r for _, s, e in RBS_BLOCKS for r in range(s, e + 1)}

def chain_resid_key(atom) -> Tuple[str, int]:
    chain = (atom.segid or atom.chainID or "").strip()
    return chain, atom.resid

def ca_coord_map(universe: mda.Universe) -> Dict[Tuple[str, int], np.ndarray]:
    return {chain_resid_key(atom): atom.position.copy() for atom in universe.select_atoms("name CA and protein")}

def shared_coords(coords1, coords2, subset_resids: Optional[set] = None):
    keys = set(coords1) & set(coords2)
    if subset_resids is not None:
        keys = {k for k in keys if k[1] in subset_resids}
    if not keys:
        return None, None
    keys = sorted(keys)
    x1 = np.vstack([coords1[k] for k in keys])
    x2 = np.vstack([coords2[k] for k in keys])
    return x1, x2

def kabsch_fit(P, Q):
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    C = Pc.T @ Qc
    V, _, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V) * np.linalg.det(Wt))
    U = V @ np.diag([1.0, 1.0, d]) @ Wt
    Pr = Pc @ U + Q.mean(axis=0)
    rmsd = np.sqrt(((Pr - Q) ** 2).mean())
    return rmsd, U, P.mean(axis=0), Q.mean(axis=0)

def rmsd_after_global_fit(P_subset, U, P_centroid, Q_centroid, Q_subset):
    Pr = (P_subset - P_centroid) @ U + Q_centroid
    return float(np.sqrt(((Pr - Q_subset) ** 2).mean()))

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pdb-dir", default="best_pdbs")
    p.add_argument("--output-dir", default="results")
    p.add_argument("--reference-substring", default="EPI_ISL_18133029")
    return p.parse_args()

def main():
    args = parse_args()
    pdb_dir = Path(args.pdb_dir)
    output_dir = Path(args.output_dir)

    warnings.filterwarnings("ignore")

    if not pdb_dir.is_dir():
        sys.exit(f"ERROR: missing PDB directory: {pdb_dir}")

    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    if len(pdb_files) < 2:
        sys.exit("ERROR: need at least two PDB files")

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {len(pdb_files)} PDBs from {pdb_dir}", flush=True)
    universes = [mda.Universe(str(p), format="PDB", in_memory=True) for p in pdb_files]
    names = [p.stem for p in pdb_files]
    coord_maps = [ca_coord_map(u) for u in universes]

    try:
        ref_idx = next(i for i, name in enumerate(names) if args.reference_substring in name)
    except StopIteration:
        sys.exit(f"ERROR: reference substring '{args.reference_substring}' not found in PDB names")

    print(f"Reference model: {names[ref_idx]}", flush=True)

    n = len(names)
    rmsd_all = np.full((n, n), np.nan, dtype=np.float32)
    rmsd_rbs = np.full((n, n), np.nan, dtype=np.float32)

    for i in range(n):
        rmsd_all[i, i] = 0.0
        rmsd_rbs[i, i] = 0.0
        for j in range(i + 1, n):
            P, Q = shared_coords(coord_maps[i], coord_maps[j])
            if P is None:
                continue
            rmsd_global, U, Pc, Qc = kabsch_fit(P, Q)
            rmsd_all[i, j] = rmsd_all[j, i] = rmsd_global

            Psub, Qsub = shared_coords(coord_maps[i], coord_maps[j], RBS_RESIDS)
            if Psub is not None:
                rmsd_local = rmsd_after_global_fit(Psub, U, Pc, Qc, Qsub)
                rmsd_rbs[i, j] = rmsd_rbs[j, i] = rmsd_local

    pd.DataFrame(rmsd_all, index=names, columns=names).to_csv(output_dir / "rmsd_all.csv")
    pd.DataFrame(rmsd_rbs, index=names, columns=names).to_csv(output_dir / "rmsd_rbs.csv")

    ref_all = np.full(n, np.nan, dtype=np.float32)
    ref_rbs = np.full(n, np.nan, dtype=np.float32)
    for j in range(n):
        P, Q = shared_coords(coord_maps[ref_idx], coord_maps[j])
        if P is None:
            continue
        rmsd_global, U, Pc, Qc = kabsch_fit(P, Q)
        ref_all[j] = rmsd_global
        Psub, Qsub = shared_coords(coord_maps[ref_idx], coord_maps[j], RBS_RESIDS)
        if Psub is not None:
            ref_rbs[j] = rmsd_after_global_fit(Psub, U, Pc, Qc, Qsub)

    pd.DataFrame({"name": names, "rmsd_to_ref_all": ref_all, "rmsd_to_ref_rbs": ref_rbs}).to_csv(
        output_dir / "rmsd_to_ref.csv", index=False
    )

    print(f"Wrote RMSD outputs to {output_dir}", flush=True)

if __name__ == "__main__":
    main()
