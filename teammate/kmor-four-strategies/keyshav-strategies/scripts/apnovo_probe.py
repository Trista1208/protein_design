#!/usr/bin/env python
"""Run the organisers' ancestral scorer on the supplied theozyme.

The addendum says one call settles the angle2/dihedral2 question. This is it.

Builds an alphafold3 Structure directly from arrays (from_atom_arrays) rather
than via mmCIF, because gemmi's CIF writes '.' into _entity_poly_seq and AF3's
parser rejects it. Building the arrays also lets us control atom_b_factor,
which is exactly the field the motif-inference trap keys on.

Run with .venv-apnovo/bin/python
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from alphafold3 import structure
from alphaprotein_novo.metrics import kemp_metrics

PDB = Path(sys.argv[1] if len(sys.argv) > 1
           else "theozymes-discovered/Theozyme_withORI_bfix.pdb")

PROTEIN = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
           "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


def read_pdb(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        rows.append(dict(
            chain=line[21:22].strip() or "A",
            res_id=int(line[22:26]),
            res_name=line[17:20].strip(),
            atom_name=line[12:16].strip(),
            element=(line[76:78].strip() or line[12:16].strip()[0]),
            x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54]),
            b=float(line[60:66] or 0.0), occ=float(line[54:60] or 1.0),
        ))
    return rows


def build(rows):
    # protein chain first, ligand as its own chain - the ligand MUST be a
    # separate entity or there is no interface to measure
    prot = [r for r in rows if r["res_name"] in PROTEIN]
    lig = [r for r in rows if r["res_name"] not in PROTEIN]
    # PRESERVE the original protein numbering. Renumbering 1..N breaks
    # motif_residues: callers refer to residues by their real numbers
    # (Asp127/Gln50 from the crystal, or A95/A39 from RFD3's
    # diffused_index_map), and a renumbered structure silently has no such
    # residue - APNovo then reports
    #   "No valid catalytic role combinations found for residues [127, 50]".
    seen, order = {}, []
    for r in prot:
        key = (r["chain"], r["res_id"])
        if key not in seen:
            seen[key] = r["res_id"]
            order.append(key)
        r["_new"] = seen[key]
    # Each distinct non-protein residue needs its OWN res_id. Assigning 1 to
    # all of them collides as soon as the structure has more than one
    # heteroatom residue (a crystal chain usually does), and AF3 then raises
    # "Inconsistent res_name, res_id and all_residues".
    lseen = {}
    for r in lig:
        key = (r["chain"], r["res_id"], r["res_name"])
        if key not in lseen:
            lseen[key] = len(lseen) + 1
        r["_new"] = lseen[key]

    allr = prot + lig
    kw = dict(
        res_id=np.array([r["_new"] for r in allr], dtype=np.int32),
        # AF3 requires object dtype for string arrays, not numpy <U
        chain_id=np.array(["A"] * len(prot) + ["L"] * len(lig), dtype=object),
        chain_type=np.array(["polypeptide(L)"] * len(prot)
                            + ["non-polymer"] * len(lig), dtype=object),
        res_name=np.array([r["res_name"] for r in allr], dtype=object),
        atom_name=np.array([r["atom_name"] for r in allr], dtype=object),
        atom_element=np.array([r["element"] for r in allr], dtype=object),
        atom_x=np.array([r["x"] for r in allr]),
        atom_y=np.array([r["y"] for r in allr]),
        atom_z=np.array([r["z"] for r in allr]),
        atom_b_factor=np.array([r["b"] for r in allr]),
        atom_occupancy=np.array([r["occ"] for r in allr]),
        name="theozyme",
    )
    return structure.from_atom_arrays(**kw), order


def main() -> int:
    rows = read_pdb(PDB)
    st, order = build(rows)
    print(f"  {PDB.name}: {len(rows)} atoms")
    print(f"  protein residue map (orig -> new): "
          f"{[(f'{c}{i}', n + 1) for n, (c, i) in enumerate(order)]}")
    print(f"  res_name set: {sorted(set(st.res_name))}")
    nz = int((np.asarray(st.atom_b_factor) > 0).sum())
    print(f"  atoms with B > 0: {nz}/{len(st.atom_b_factor)}"
          f"   <- motif inference keys on this")

    motif = list(range(1, len(order) + 1))        # all protein residues are motif
    for spec in (motif, None):
        label = "explicit motif_residues=" + str(spec) if spec else "motif_residues=None (B-factor inference)"
        try:
            m = kemp_metrics.kemp_structure_metrics(
                st, motif_residues=spec, ligand_res_name="6NT")
        except Exception as exc:                      # noqa: BLE001
            print(f"\n  === {label} -> {type(exc).__name__}: {str(exc)[:200]}")
            continue
        print(f"\n  === {label} -> OK ===")
        for k in sorted(m):
            v = m[k]
            if isinstance(v, (int, float, np.floating)):
                print(f"    {k:46} {float(v):10.4f}")
            else:
                print(f"    {k:46} {str(v)[:60]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
