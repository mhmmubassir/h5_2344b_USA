#!/usr/bin/env python3
"""
rosetta2amber_glycan.py • v8 (03‑May‑2025)

• Converts Rosetta‑style glycan residue/atom names to Amber / GLYCAM.
• Inserts ROH with anomeric atom “O1” and a TER record before each glycan tree.
• Residue mapping implemented in this version
    Neu5Ac / SIA / 0SA  ➜ 0SA
    Gal                 ➜ 3LB (α2‑3) or 6LB (α2‑6) depending on linkage geometry
    Glc                 ➜ 4GB (plain β‑D‑Glcp)
    GlcNAc (NAG)        ➜ 3YB (β‑D‑GlcpNAc)
• Multi‑chain aware; original chain IDs preserved.
• HETATM → ATOM, atom names remapped where GLYCAM differs.
• Serial numbers and residue indices are rebuilt so protein and each glycan tree
  continue in order.
• <output>.changes.txt lists only the lines that changed for easy diffing.

Usage
-----
$ python rosetta2amber_glycan.py input.pdb output.pdb
"""

import sys
from collections import OrderedDict, defaultdict
from math import sqrt

# ─── residue categories ──────────────────────────────────────────────
GLC_CODES = {"Glc", "GLC"}          # → 4GB (plain glucose)
NAG_CODES = {"NAG"}                 # → 3YB (GlcNAc)

GAL_CODES = {"Gal", "GAL", "3LB", "6LB"}
SIA_CODES = {"Neu", "NEU", "SIA", "0SA"}
SUGAR_CODES = GLC_CODES | NAG_CODES | GAL_CODES | SIA_CODES

# ─── atom‑rename templates ───────────────────────────────────────────
# NB: 4GB matches GLYCAM names already – no mapping needed.
MAP_4YA = {
    "C10": "C5N", "C11": "CME", "O10": "O5N", "O11": "OME",
    "C7": "C2N", "O7": "O2N", "C8": "CME",
    "CN2": "C2N", "CAN2": "CME", "OCN2": "O2N",
}
MAP_GAL = {
    "C10": "C5N", "C11": "CME", "O10": "O5N", "O11": "OME",
    "C7": "C2N", "O7": "O2N", "C8": "CME",
}
MAP_0SA = {
    "C10": "C5N", "O10": "O5N", "C11": "CME",
    "1O1": "O1A", "2O1": "O1B",
}

ATOM_REMAP = {
    "4YA": MAP_4YA, "4YB": MAP_4YA,
    "3YB": MAP_4YA,             # GlcNAc (β)
    "3LB": MAP_GAL, "6LB": MAP_GAL,
    "0SA": MAP_0SA,
}

# ─────────────────────────────────────────────────────────────────────

def parse_atom(line):
    """Return a dict of the PDB atom line fields we care about."""
    return {
        "record": line[:6],
        "serial": int(line[6:11]),
        "atom":   line[12:16],
        "alt":    line[16],
        "resname": line[17:20],
        "chain":  line[21],
        "resnum": int(line[22:26]),
        "icode":  line[26],
        "x": float(line[30:38]),
        "y": float(line[38:46]),
        "z": float(line[46:54]),
        "xyz": line[30:54],
        "occ": line[54:60],
        "beta": line[60:66],
        "element": line[76:78] if len(line) >= 78 else "  ",
        "charge":  line[78:80] if len(line) >= 80 else "  ",
        "raw": line.rstrip("\n"),
    }


def fmt_atom(f, serial, resnum, resname, atom):
    """Return a formatted PDB ATOM line based on field template *f*."""
    return (f"ATOM  {serial:>5} {atom:<4}{f['alt']}{resname:>3} "
            f"{f['chain']}{resnum:>4}{f['icode']}   {f['xyz']}{f['occ']}{f['beta']}          "
            f"{f['element']:<2}{f['charge']}\n")


def dist(a, b):
    return sqrt((a['x']-b['x'])**2 + (a['y']-b['y'])**2 + (a['z']-b['z'])**2)


def collapse_ter(lines):
    """Remove duplicate sequential TER records."""
    out, prev = [], False
    for ln in lines:
        if ln.startswith("TER"):
            if not prev:
                out.append(ln)
            prev = True
        else:
            out.append(ln)
            prev = False
    return out


def pick_roh_atom(res_atoms):
    """Return the atom dictionary to become ROH‑O1 (heuristic)."""
    for tag in ("O1", " O1", "O5", "O11", " O "):
        for a in res_atoms:
            if a["atom"].strip() == tag.strip():
                return a
    # Fallback: O closest in serial to C1
    c1 = next((x for x in res_atoms if x["atom"].strip() == "C1"), None)
    os = [x for x in res_atoms if x["element"].strip() == "O"]
    if c1 and os:
        return min(os, key=lambda o: abs(o["serial"]-c1["serial"]))
    return None


# ─── main conversion routine ─────────────────────────────────────────

def main(inp, out):
    with open(inp) as fh:
        raw = fh.readlines()

    residues = OrderedDict()      # (chain, resnum) → [atom, …]
    order = []                    # preserve original order
    chain_break_flags = {}        # mark first glycan residue of each tree
    last_key, last_was_sugar, last_chain = None, False, None
    prot_last_resnum = 0

    # Pass 1 – gather atoms per residue and detect tree starts
    for ln in raw:
        tag = ln[:6].strip()
        if tag not in ("ATOM", "HETATM"):
            last_was_sugar = False
            continue
        a = parse_atom(ln)
        key = (a["chain"], a["resnum"])
        is_sugar = a["resname"].strip() in SUGAR_CODES
        residues.setdefault(key, []).append(a)

        if key != last_key:
            order.append(key)
            if is_sugar and ((not last_was_sugar) or a["chain"] != last_chain):
                chain_break_flags[key] = True
        if not is_sugar and a["resnum"] > prot_last_resnum:
            prot_last_resnum = a["resnum"]
        last_key, last_was_sugar, last_chain = key, is_sugar, a["chain"]

    out_lines, log = [], []
    serial = 1
    resnum = prot_last_resnum
    wrote_prot_ter = False

    # Pass 2 – write updated records
    for key in order:
        chain, old_rn = key
        atoms = residues[key]
        is_sugar = atoms[0]["resname"].strip() in SUGAR_CODES

        # ───── Protein residues: copy with serial‑renumbering ─────
        if not is_sugar:
            for a in atoms:
                out_lines.append(fmt_atom(a, serial, a["resnum"], a["resname"], a["atom"].strip()))
                serial += 1
            wrote_prot_ter = False
            continue

        # ───── First residue of glycan tree ─────
        if chain_break_flags.get(key, False):
            if not wrote_prot_ter:
                out_lines.append("TER\n")
                wrote_prot_ter = True
            # insert ROH
            roh_atom = pick_roh_atom(atoms)
            if roh_atom is None:
                sys.exit(f"❌  Anomeric O not found in residue {atoms[0]['resname']} {chain}{old_rn}")
            atoms.remove(roh_atom)
            resnum += 1
            roh_line = fmt_atom(roh_atom, serial, resnum, "ROH", "O1  ")
            out_lines += [roh_line, "TER\n"]
            log.append(f"{roh_atom['raw']}  →  {roh_line.strip()}")
            serial += 1

        # ───── Determine target residue name ─────
        old_name = atoms[0]["resname"].strip()

        if old_name in GLC_CODES:          # plain glucose
            tgt = "4GB"
        elif old_name in NAG_CODES:        # GlcNAc
            tgt = "3YB"
        elif old_name in GAL_CODES:        # Gal – decide linkage
            next_idx = order.index(key) + 1
            if next_idx < len(order):
                next_chain, _ = order[next_idx]
                next_atoms = residues[order[next_idx]]
                if next_chain == chain and next_atoms[0]["resname"].strip() in SIA_CODES:
                    neu_C2 = next(a for a in next_atoms if a["atom"].strip() == "C2")
                    O3 = next((x for x in atoms if x["atom"].strip() == "O3"), None)
                    O6 = next((x for x in atoms if x["atom"].strip() == "O6"), None)
                    tgt = "3LB" if (O3 and O6 and dist(neu_C2, O3) < dist(neu_C2, O6)) else "6LB"
                else:
                    tgt = "3LB"
        elif old_name in SIA_CODES:
            tgt = "0SA"
        else:
            tgt = old_name  # should not happen

        amap = ATOM_REMAP.get(tgt, {})
        resnum += 1
        for a in atoms:
            old_atom = a["atom"].strip()
            if chain_break_flags.get(key, False) and old_atom == "O1":
                continue  # anomeric O was moved to ROH
            new_atom = amap.get(old_atom, old_atom)
            line = fmt_atom(a, serial, resnum, tgt, new_atom)
            if (old_atom != new_atom) or (old_name != tgt) or a["record"].startswith("HETATM"):
                log.append(f"{a['raw']}  →  {line.strip()}")
            out_lines.append(line)
            serial += 1
        out_lines.append("TER\n")
        wrote_prot_ter = True

    # ─── Finish ──────────────────────────────────────────────────────
    out_lines = collapse_ter(out_lines)
    with open(out, "w") as fh:
        fh.writelines(out_lines)
    with open(out + ".changes.txt", "w") as fh:
        fh.write("\n".join(log))

    print(f"✔  Wrote '{out}' and '{out}.changes.txt' (modified {len(log)} atoms).")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rosetta2amber_glycan.py  input.pdb  output.pdb")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
