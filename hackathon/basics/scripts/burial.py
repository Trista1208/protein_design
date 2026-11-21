#!/usr/bin/env python3
"""Catalytic-base burial as defined on the slides:

  Number of CA atoms close to the catalytic base (up to 12 Å)

We count protein CA atoms within 12 Å of the base residue's CA.
The slides did not give a numeric pass cutoff, so this prints the count.
On 4BS0 A ASP127 the crystal count (excluding self) is 47.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from pdb_io import dist, find, one, parse_atoms  # noqa: E402

CUTOFF = 12.0


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("pdb", type=Path)
    p.add_argument("--chain", default="A")
    p.add_argument("--base-res", type=int, required=True)
    p.add_argument("--cutoff", type=float, default=CUTOFF)
    args = p.parse_args()

    atoms = parse_atoms(args.pdb, chain=args.chain)
    base_ca = one(atoms, ri=args.base_res, an="CA", het=False)
    cas = find(atoms, an="CA", het=False)
    n_incl = sum(1 for c in cas if dist(c, base_ca) <= args.cutoff)
    report = {
        "base": f"{args.chain}:{base_ca['rn']}{args.base_res}",
        "cutoff_A": args.cutoff,
        "n_CA_within_cutoff_including_self": n_incl,
        "n_CA_within_cutoff_excluding_self": n_incl - 1,
        "crystal_4BS0_A_ASP127_excluding_self": 47,
        "note": "Higher = more buried. No official minimum was published; compare to 4BS0 (~47) and reject obviously exposed bases (very small counts).",
    }
    json.dump(report, sys.stdout, indent=2)
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
