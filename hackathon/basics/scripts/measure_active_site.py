#!/usr/bin/env python3
"""Report catalytic geometry. Default residue IDs are 4BS0's.

Prints distances only. Does not invent a pass/fail cutoff — the organizers
did not publish numeric geometry thresholds in the slides.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pdb_io import dist, one, parse_atoms  # noqa: E402

CRYSTAL = {
    "ASP.OD2__6NT.N3": 2.53,
    "GLN.NE2__6NT.N1": 3.10,
}


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("pdb", type=Path)
    p.add_argument("--chain", default="A")
    p.add_argument("--base-res", type=int, default=127, help="Asp residue number")
    p.add_argument("--donor-res", type=int, default=50, help="Gln residue number")
    p.add_argument("--ligand-resname", default="6NT")
    args = p.parse_args()

    atoms = parse_atoms(args.pdb, chain=args.chain)
    try:
        od2 = one(atoms, ri=args.base_res, an="OD2", rn="ASP")
        ne2 = one(atoms, ri=args.donor_res, an="NE2", rn="GLN")
        n3 = one(atoms, rn=args.ligand_resname, an="N3")
        n1 = one(atoms, rn=args.ligand_resname, an="N1")
    except KeyError as e:
        print(f"ERROR {e}", file=sys.stderr)
        print("Pass --base-res / --donor-res / --ligand-resname for your design numbering.", file=sys.stderr)
        return 1

    measured = {
        "ASP.OD2__6NT.N3": round(dist(od2, n3), 3),
        "GLN.NE2__6NT.N1": round(dist(ne2, n1), 3),
        "crystal_4BS0_A": CRYSTAL,
        "delta_vs_4BS0": {
            "ASP.OD2__6NT.N3": round(dist(od2, n3) - CRYSTAL["ASP.OD2__6NT.N3"], 3),
            "GLN.NE2__6NT.N1": round(dist(ne2, n1) - CRYSTAL["GLN.NE2__6NT.N1"], 3),
        },
    }
    json.dump(measured, sys.stdout, indent=2)
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
