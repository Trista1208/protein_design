#!/usr/bin/env python3
"""Extract the 4BS0 chain-A Kemp theozyme: GLN50 + ASP127 + ligand 6NT.

Residue numbers and ligand id were read from the RCSB PDB file, not guessed.
Chain B 6NT is disordered (altloc A/B) — we use chain A only.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pdb_io import dist, find, one, parse_atoms, write_pdb  # noqa: E402

REF = HERE.parent / "references"
OUT_DIR = HERE.parent / "theozyme"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--pdb", type=Path, default=REF / "4BS0.pdb")
    p.add_argument("--out", type=Path, default=OUT_DIR / "4bs0_A_GLN50_ASP127_6NT.pdb")
    args = p.parse_args()
    if not args.pdb.exists():
        print(f"missing {args.pdb} — run scripts/download_references.sh", file=sys.stderr)
        return 1

    atoms = parse_atoms(args.pdb, chain="A")
    gln = [a for a in find(atoms, ri=50, rn="GLN", het=False) if a["heavy"]]
    asp = [a for a in find(atoms, ri=127, rn="ASP", het=False) if a["heavy"]]
    lig = [a for a in find(atoms, rn="6NT", het=True) if a["heavy"]]
    if len(gln) < 5 or len(asp) < 5 or len(lig) < 8:
        print("extract failed: residue identities do not match 4BS0 A (GLN50, ASP127, 6NT)", file=sys.stderr)
        return 1

    geom = {
        "source": "PDB 4BS0 chain A (HG3.17 + 6-nitrobenzotriazole)",
        "base": "A:ASP127",
        "donor": "A:GLN50",
        "ligand": "A:6NT:1305",
        "distances_A": {
            "ASP127.OD2__6NT.N3": round(dist(one(asp, an="OD2"), one(lig, an="N3")), 3),
            "ASP127.OD2__6NT.N2": round(dist(one(asp, an="OD2"), one(lig, an="N2")), 3),
            "GLN50.NE2__6NT.N1": round(dist(one(gln, an="NE2"), one(lig, an="N1")), 3),
        },
        "note": (
            "OD2–N3 ~2.5 Å is the crystal contact of the catalytic carboxylate "
            "to the triazole. NE2–N1 ~3.1 Å is the Gln donor to 6NT. "
            "These are experimental HG3.17 numbers, not DFT transition-state coordinates. "
            "Do not submit 4BS0 as a design (TM / identity gates)."
        ),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_pdb(
        args.out,
        gln + asp + lig,
        remarks=[
            "4BS0 A GLN50 ASP127 6NT heavy atoms",
            "teacher motif only — generate a NEW scaffold around this",
        ],
    )
    json_path = args.out.with_suffix(".json")
    json_path.write_text(json.dumps(geom, indent=2) + "\n")
    print(f"wrote {args.out}")
    print(f"wrote {json_path}")
    print(json.dumps(geom["distances_A"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
