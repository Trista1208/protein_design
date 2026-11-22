"""Emit a Rosetta enzdes constraint file for the Asp/Gln/6NT theozyme.

WHY THIS IS THE TOOL THAT MATTERS
---------------------------------
The 84-variant screen showed every hit fixing ONE term and breaking the other:

    F26Y    d_oxy 9.419 (bad)   dihedral  16.4 (IN BAND)
    H103F   d_oxy 2.899 (GOOD)  dihedral  93.0 (dead zone)
    F42Y    d_oxy 2.877 (GOOD)  dihedral  83.1 (dead zone)

That is LESSONS 1.3's orthogonality: distance and orientation are independent,
and neither LigandMPNN (which optimises sequence likelihood) nor a steric
ladder (which optimises nothing) can aim at both. A Rosetta enzdes constraint
block scores distance, angle AND torsion into one energy, which is exactly how
denovoKemp's `theozime.cst` holds its catalytic geometry.

ATOM CONVENTION, verified rather than assumed
---------------------------------------------
Rosetta enzdes names three atoms per partner: the ligand ("downstream") gets
A1 A2 A3, the protein residue ("upstream") gets B1 B2 B3, and the six
constraints are then

    distanceAB  A1-B1
    angle_A     A2-A1-B1
    angle_B     A1-B1-B2
    torsion_A   A3-A2-A1-B1
    torsion_AB  A2-A1-B1-B2
    torsion_B   A1-B1-B2-B3

Confirmed against denovoKemp's own file, whose block 1 maps
`angle_A = C1-C7-OD`, `angle_B = C7-OD-CG`, `torsion_A = N1-C1-C7-OD`,
`torsion_AB = C1-C7-OD-CG`, `torsion_B = C7-OD-CG-OD2` for ligand atoms
`C7 C1 N1` and an OOC-type protein oxygen. That is the mapping above.

WHY THIS MAPS ONTO OUR SCORED TERMS EXACTLY
-------------------------------------------
With ligand atoms `N1 N2 N3` and Gln `NE2 CD OE1`:
    distanceAB  = N1...NE2        -> our d_oxy
    angle_A     = N2-N1-NE2       -> our angle2_oxy (centred on the acceptor)
    torsion_AB  = N2-N1-NE2-CD    -> our dihedral12_oxy, same four atoms
and with ligand `N3 N2 N1` and Asp `OD2 CG OD1`:
    distanceAB  = N3...OD2        -> our d_base
    angle_A     = N2-N3-OD2       -> our angle2_base
    torsion_AB  = N2-N3-OD2-CG    -> our dihedral12_base

So satisfying this file is, atom for atom, satisfying the scored geometry.

TARGETS come from official.TARGETS (crystal-measured, 4BS0/7K4Z) and the
tolerances from the scoring windows - NOT from denovoKemp, whose numbers are
on a different ligand with a different contact atom and do not transfer
(LOOP_DESIGN 3.2).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C                                    # noqa: E402
from kemp.official import (LAUKO, TARGETS, _angle, _dihedral)   # noqa: E402
from kemp.theozyme import load                                  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "outputs/rosetta"

#: (tag, ligand A1 A2 A3, protein resname, B1 B2 B3, residue1 code, targets)
BLOCKS = [
    {"tag": "catalytic base: Asp carboxylate to 6NT N3", "key": "base",
     "lig": ("N3", "N2", "N1"), "res3": "ASP", "res1": "D",
     "prot": ("OD2", "CG", "OD1"),
     "dist_target": TARGETS["d_base"]["target"], "dist_tol": 0.30,
     "dist_k": 100.0, "ang_k": 60.0, "tor_k": 60.0},
    {"tag": "oxyanion donor: Gln amide to 6NT N1", "key": "oxy",
     "lig": ("N1", "N2", "N3"), "res3": "GLN", "res1": "Q",
     "prot": ("NE2", "CD", "OE1"),
     "dist_target": TARGETS["d_oxy"]["target"], "dist_tol": 0.30,
     "dist_k": 100.0, "ang_k": 60.0, "tor_k": 60.0},
]

ANGLE_TARGET = TARGETS["angle"]["target"]
#: window 80-160 about 120 -> 20 deg keeps us comfortably inside
ANGLE_TOL = 20.0

#: TORSION TARGETS ARE TAKEN FROM THE CRYSTALS, NOT FROM 180.
#:
#: There are TWO allowed sp2 bands, 130-230 and 310-50. An earlier version of
#: this file targeted 180 - the centre of the first band - which would have
#: pulled the geometry AWAY from where every real structure sits. Measured:
#:
#:   torsion_AB (oxy)   4BS0  2.74 deg   s=0.954   real refined Gln50
#:                      7K4Z 355.58 deg  s=0.926   real refined Gln50
#:                      theozyme 50.90 deg  s=0.000  MODELLED Gln50
#:   torsion_AB (base)  4BS0 338.96 / 7K4Z 340.78 / theozyme 328.13
#:
#: Both crystals use the SECOND band, near 0/360. The supplied theozyme's
#: modelled Gln50 - the one with B-factor 0.00 on all 17 atoms, computationally
#: reverted from the Q50M mutant - sits ~48 deg away from the crystallographic
#: rotamer and is the ONLY one of the three that scores zero.
#:
#: Consequence, and it is the important one: every design scaffolded on this
#: theozyme inherits a donor torsion that cannot score. Aiming such a design at
#: the NEAREST band edge (230) earns s_dih ~ 0.167; aiming at the crystal
#: rotamer earns ~1.0. Target the crystals.
TORSION_TARGETS = {"base": 339.9, "oxy": 359.2}    # circular means, measured
TORSION_TOL = 40.0


def main() -> int:
    tz = load(C.THEOZYME_DESIGN if C.THEOZYME_DESIGN.exists()
              else C.THEOZYME_WITH_ORI)
    OUT.mkdir(parents=True, exist_ok=True)

    lines = [
        "# Kemp eliminase theozyme: Asp general base + Gln oxyanion donor + 6NT.",
        "# Generated by scripts/make_enzdes_cst.py - do not hand-edit.",
        "#",
        "# Targets are the CRYSTAL-measured values in kemp/official.py TARGETS",
        "# (4BS0/7K4Z), with tolerances set from the scoring windows. The",
        "# denovoKemp constraint numbers are NOT used: their base contacts a",
        "# ring CARBON on the substrate, ours contacts a ring NITROGEN on the",
        "# TSA, so the values do not transfer (LOOP_DESIGN 3.2).",
        "#",
        "# Rosetta enzdes atom convention (A = ligand, B = protein):",
        "#   distanceAB A1-B1   angle_A A2-A1-B1   angle_B A1-B1-B2",
        "#   torsion_A A3-A2-A1-B1   torsion_AB A2-A1-B1-B2   torsion_B A1-B1-B2-B3",
        "#",
        "# columns: target  tolerance  penalty_k  periodicity  n_samples",
        "",
    ]

    print("measured on the supplied theozyme (the trusted reference):")
    for b in BLOCKS:
        a1, a2, _a3 = b["lig"]
        b1, b2, _b3 = b["prot"]
        try:
            xa1 = tz.atom(C.LIGAND_CCD, a1).xyz
            xa2 = tz.atom(C.LIGAND_CCD, a2).xyz
            xb1 = tz.atom(b["res3"], b1).xyz
            xb2 = tz.atom(b["res3"], b2).xyz
        except Exception as exc:                                # noqa: BLE001
            print(f"  {b['tag']}: ATOM MISSING -> {exc}")
            return 2
        d = float(np.linalg.norm(xa1 - xb1))
        ang = _angle(xa2, xa1, xb1)            # A2-A1-B1
        tor = _dihedral(xa2, xa1, xb1, xb2)    # A2-A1-B1-B2
        print(f"  {b['tag']}")
        print(f"     distance {a1}...{b1} = {d:.3f} A   (target {b['dist_target']})")
        print(f"     angle_A  {a2}-{a1}-{b1} = {ang:.2f} deg   (target {ANGLE_TARGET})")
        print(f"     torsion_AB {a2}-{a1}-{b1}-{b2} = {tor:.2f} deg "
              f"  (crystal target {TORSION_TARGETS[b['key']]})")

        lines += [
            f"# {b['tag']}",
            f"#   measured on the theozyme: d={d:.3f} A, angle_A={ang:.2f}, "
            f"torsion_AB={tor:.2f}",
            "CST::BEGIN",
            f"  TEMPLATE::   ATOM_MAP: 1 atom_name: {' '.join(b['lig'])}",
            f"  TEMPLATE::   ATOM_MAP: 1 residue3: {C.LIGAND_CCD}",
            f"  TEMPLATE::   ATOM_MAP: 2 atom_name: {' '.join(b['prot'])}",
            f"  TEMPLATE::   ATOM_MAP: 2 residue1: {b['res1']}",
            "",
            f"  CONSTRAINT:: distanceAB: {b['dist_target']:7.2f} "
            f"{b['dist_tol']:6.2f} {b['dist_k']:7.1f}  0     1",
            f"  CONSTRAINT::    angle_A: {ANGLE_TARGET:7.2f} "
            f"{ANGLE_TOL:6.2f} {b['ang_k']:7.1f}  360.  1",
            f"  CONSTRAINT:: torsion_AB: {TORSION_TARGETS[b['key']]:7.2f} "
            f"{TORSION_TOL:6.2f} {b['tor_k']:7.1f}  360.  1",
            "CST::END",
            "",
        ]

    path = OUT / "kemp_theozyme.cst"
    path.write_text("\n".join(lines))
    print(f"\nwrote {path}")
    print("\nNOTE: this constrains BOTH distance and torsion in one energy, "
          "which is the thing\nneither LigandMPNN nor the steric ladder could do.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
