#!/usr/bin/env python3
"""
Convert Rosetta-style glycoprotein PDB files into Amber/GLYCAM-style PDBs.

Key behavior:
- keeps protein residue numbering unchanged
- renumbers inserted ROH and glycan residues after the protein
- inserts ROH O1 at each glycan-tree root
- converts common sugar residue names to GLYCAM names
- renames selected sugar atoms to GLYCAM-compatible names
- writes ATOM records and removes hydrogens
- writes a per-file change log

Examples:
    python rosetta_to_amber.py
    python rosetta_to_amber.py --ref reference.pdb
    python rosetta_to_amber.py model1.pdb model2.pdb --out-dir amber_glycan_pdbs
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import OrderedDict
from math import sqrt
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

STD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

GLC_CODES = {"GLC", "Glc"}
NAG_CODES = {"NAG", "Nag"}
GAL_CODES = {"GAL", "Gal"}
SIA_CODES = {"SIA", "NEU", "Neu", "0SA"}

MAP_4YA = {
    "C10": "C5N", "O10": "O5N", "C11": "CME", "O11": "OME",
    "C7": "C2N", "O7": "O2N", "C8": "CME",
    "CN2": "C2N", "CAN2": "CME", "OCN2": "O2N",
}
MAP_GAL = dict(MAP_4YA)
MAP_0SA = {
    "C10": "C5N", "O10": "O5N", "C11": "CME",
    "1O1": "O1A", "2O1": "O1B",
    "CN5": "C5N", "OCN5": "O5N", "CAN5": "CME",
}
ATOM_REMAP = {
    "4YA": MAP_4YA,
    "4YB": MAP_4YA,
    "3YB": MAP_4YA,
    "3LB": MAP_GAL,
    "6LB": MAP_GAL,
    "0SA": MAP_0SA,
}


class Atom:
    __slots__ = (
        "raw", "rec", "serial", "name", "alt", "res", "chain", "resnum", "icode",
        "x", "y", "z", "xyz", "occ", "beta", "el", "charge"
    )

    def __init__(self, line: str) -> None:
        self.raw = line.rstrip("\n")
        self.rec = line[:6]
        self.serial = int(line[6:11])
        self.name = line[12:16]
        self.alt = line[16]
        self.res = line[17:20].strip()
        self.chain = line[21]
        self.resnum = int(line[22:26])
        self.icode = line[26]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.xyz = line[30:54]
        self.occ = line[54:60]
        self.beta = line[60:66]
        self.el = line[76:78] if len(line) >= 78 else "  "
        self.charge = line[78:80] if len(line) >= 80 else "  "

    def format_pdb(self, serial: int, resnum: int, resname: str, atom_name: str) -> str:
        return (
            f"ATOM  {serial:>5} {atom_name:<4}{self.alt}{resname:>3} "
            f"{self.chain}{resnum:>4}{self.icode}   {self.xyz}{self.occ}{self.beta}          "
            f"{self.el:<2}{self.charge}\n"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdbs", nargs="*", help="Input PDB files. If omitted, all .pdb files in the current directory are used.")
    parser.add_argument("--ref", help="Reference PDB with desired ligand residue codes.")
    parser.add_argument("--out-dir", default="amber_glycan_pdbs", help="Output directory for converted PDBs.")
    parser.add_argument("--log-dir", default="amber_glycan_logs", help="Output directory for change logs.")
    return parser.parse_args()


def is_ligand(resname: str) -> bool:
    return resname not in STD_AA


def is_hydrogen(atom: Atom) -> bool:
    return atom.el.strip().upper() == "H" or atom.name.strip().upper().startswith("H")


def distance(atom_a: Atom, atom_b: Atom) -> float:
    return sqrt((atom_a.x - atom_b.x) ** 2 + (atom_a.y - atom_b.y) ** 2 + (atom_a.z - atom_b.z) ** 2)


def is_hidden_glcnac(atoms: Iterable[Atom]) -> bool:
    acetamide_atoms = {"CN2", "CAN2", "OCN2", "C2N", "O2N", "C5N"}
    return any(atom.name.strip() in acetamide_atoms for atom in atoms)


def true_resname(old_resname: str, atoms: List[Atom]) -> str:
    if old_resname in GLC_CODES and is_hidden_glcnac(atoms):
        return "NAG"
    return old_resname


def load_reference_map(ref_pdb: str | None) -> Dict[Tuple[str, str], str]:
    reference_map: Dict[Tuple[str, str], str] = {}
    if not ref_pdb:
        return reference_map

    ref_path = Path(ref_pdb)
    if not ref_path.is_file():
        raise SystemExit(f"❌ Reference PDB not found: {ref_pdb}")

    with ref_path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21]
            resid = line[22:26].strip()
            resname = line[17:20].strip().upper()
            if is_ligand(resname):
                reference_map[(chain, resid)] = resname
    return reference_map


def map_resname(chain: str, resid: str, old_resname: str, atoms: List[Atom], reference_map: Dict[Tuple[str, str], str]) -> str:
    if (chain, resid) in reference_map:
        return reference_map[(chain, resid)]
    return true_resname(old_resname, atoms).upper()


def collect_residues(pdb_path: Path) -> Tuple["OrderedDict[Tuple[str, int], List[Atom]]", List[Tuple[str, int]]]:
    residues: "OrderedDict[Tuple[str, int], List[Atom]]" = OrderedDict()
    residue_order: List[Tuple[str, int]] = []

    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom = Atom(line)
            key = (atom.chain, atom.resnum)
            residues.setdefault(key, []).append(atom)
            if len(residues[key]) == 1:
                residue_order.append(key)

    return residues, residue_order


def choose_gal_code(current_atoms: List[Atom], later_glyc_keys: List[Tuple[str, int]], residues: OrderedDict, chain: str, reference_map: Dict[Tuple[str, str], str]) -> str:
    neu_key = None
    for next_key in later_glyc_keys:
        next_chain, next_resid = next_key
        next_atoms = residues[next_key]
        next_resname = map_resname(next_chain, str(next_resid), next_atoms[0].res, next_atoms, reference_map)
        if next_chain == chain and next_resname in SIA_CODES:
            neu_key = next_key
            break

    if neu_key is None:
        return "3LB"

    c2_atoms = [atom for atom in residues[neu_key] if atom.name.strip() == "C2"]
    o3_atoms = [atom for atom in current_atoms if atom.name.strip() == "O3"]
    o6_atoms = [atom for atom in current_atoms if atom.name.strip() == "O6"]

    if not c2_atoms:
        return "3LB"
    if o3_atoms and o6_atoms:
        return "3LB" if distance(c2_atoms[0], o3_atoms[0]) < distance(c2_atoms[0], o6_atoms[0]) else "6LB"
    return "3LB"


def convert_one_pdb(pdb_path: Path, reference_map: Dict[Tuple[str, str], str], out_dir: Path, log_dir: Path) -> Path:
    residues, residue_order = collect_residues(pdb_path)

    protein_keys: List[Tuple[str, int]] = []
    glycan_keys: List[Tuple[str, int]] = []
    for key in residue_order:
        chain, resid = key
        atoms = residues[key]
        mapped_name = map_resname(chain, str(resid), atoms[0].res, atoms, reference_map)
        if is_ligand(mapped_name):
            glycan_keys.append(key)
        else:
            protein_keys.append(key)

    max_protein_resnum = max((resid for _, resid in protein_keys), default=0)
    next_resnum = max_protein_resnum
    next_serial = 1
    output_lines: List[str] = []
    change_log: List[str] = []

    for chain, resid in protein_keys:
        atoms = residues[(chain, resid)]
        mapped_resname = map_resname(chain, str(resid), atoms[0].res, atoms, reference_map)
        for atom in atoms:
            if is_hydrogen(atom):
                continue
            line = atom.format_pdb(next_serial, atom.resnum, mapped_resname, atom.name.strip())
            output_lines.append(line)
            if atom.res != mapped_resname:
                change_log.append(f"{atom.raw}  ->  {line.strip()}")
            next_serial += 1
    output_lines.append("TER\n")

    previous_chain = None
    for index, key in enumerate(glycan_keys):
        chain, resid = key
        atoms = residues[key]
        original_name = map_resname(chain, str(resid), atoms[0].res, atoms, reference_map)

        if chain != previous_chain:
            root_o1 = [atom for atom in atoms if atom.name.strip() == "O1"]
            target_atom = root_o1[0] if root_o1 else atoms[0]
            next_resnum += 1
            roh_line = target_atom.format_pdb(next_serial, next_resnum, "ROH", "O1")
            output_lines.extend(["TER\n", roh_line, "TER\n"])
            change_log.append(f"{target_atom.raw}  ->  {roh_line.strip()}")
            next_serial += 1
        previous_chain = chain

        if original_name in GLC_CODES:
            target_resname = "4GB"
        elif original_name in NAG_CODES:
            target_resname = "3YB"
        elif original_name in GAL_CODES:
            target_resname = choose_gal_code(atoms, glycan_keys[index + 1 :], residues, chain, reference_map)
        elif original_name in SIA_CODES:
            target_resname = "0SA"
        else:
            target_resname = original_name

        atom_remap = ATOM_REMAP.get(target_resname, {})
        next_resnum += 1
        for atom in atoms:
            atom_name = atom.name.strip()
            if target_resname == "0SA" and atom_name == "O1":
                continue
            if is_hydrogen(atom):
                continue
            new_atom_name = atom_remap.get(atom_name, atom_name)
            line = atom.format_pdb(next_serial, next_resnum, target_resname, new_atom_name)
            output_lines.append(line)
            if atom.rec.startswith("HETATM") or atom_name != new_atom_name or original_name != target_resname:
                change_log.append(f"{atom.raw}  ->  {line.strip()}")
            next_serial += 1
        output_lines.append("TER\n")

    cleaned_lines: List[str] = []
    previous_was_ter = False
    for line in output_lines:
        if line == "TER\n":
            if previous_was_ter:
                continue
            previous_was_ter = True
        else:
            previous_was_ter = False
        cleaned_lines.append(line)

    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    output_pdb = out_dir / f"{pdb_path.stem}_amber.pdb"
    output_log = log_dir / f"{pdb_path.stem}_amber.changes.txt"
    output_pdb.write_text("".join(cleaned_lines))
    output_log.write_text("\n".join(change_log))

    return output_pdb


def main() -> None:
    args = parse_args()
    reference_map = load_reference_map(args.ref)

    if args.pdbs:
        pdb_paths = [Path(path) for path in args.pdbs]
    else:
        pdb_paths = sorted(Path(".").glob("*.pdb"))

    if not pdb_paths:
        raise SystemExit("❌ No input PDB files found.")

    out_dir = Path(args.out_dir)
    log_dir = Path(args.log_dir)

    for pdb_path in pdb_paths:
        if args.ref and os.path.abspath(str(pdb_path)) == os.path.abspath(args.ref):
            continue
        if not pdb_path.is_file():
            raise SystemExit(f"❌ Input PDB not found: {pdb_path}")
        output_pdb = convert_one_pdb(pdb_path, reference_map, out_dir, log_dir)
        print(f"✔ {pdb_path} -> {output_pdb}")

    print("All done.")


if __name__ == "__main__":
    main()
