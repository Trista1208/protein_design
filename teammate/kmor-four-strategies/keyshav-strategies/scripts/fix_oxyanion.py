#!/usr/bin/env python
"""Rebuild the Gln50 oxyanion contact.

The situation, which is NOT a contradiction:

  * The organisers' modelled Gln50 DOES match the crystallographic Gln50 in
    4BS0/7K4Z - confirmed. So the handed-out rotamer is faithful to the real
    enzyme.
  * But the organisers' own Rosetta constraint file
    (5TS_HG3_17_enzdes_planar.cst) targets d_oxy = 2.4 A, and the supplied
    contact is 3.33 A.

Both are true: the REAL enzyme has a ~3.3 A oxyanion contact, while their
DESIGN IDEAL is 2.4 A. HG3.17 is a directed-evolution product, not an optimum.

Consequence under a geometric mean: `oxy` carries 0.20 and 3.33 A sits near the
3.6 A edge of the Lauko window, so scaffolding the file as handed out bakes a
weak high-weight term into every design. Lifting the crystallographic rotamer
would change nothing, because it is already that rotamer.

So the fix is a ROTAMER RE-OPTIMISATION toward the organisers' own design
target, keeping the Gln backbone where it is: rotate chi1/chi2/chi3, hit
d_oxy ~ 2.4 A and angle2 ~ 120 deg, reject clashes.

Writes Theozyme_withORI_oxyopt.pdb preserving frame provenance.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C                                    # noqa: E402
from kemp.theozyme import FRAME_REMARK, load                     # noqa: E402

# The organisers' cst target is 2.40 A, but that is UNREACHABLE from the
# supplied Gln backbone: |CB-N1| = 6.988 A and the side chain's fully extended
# reach is 4.407 A, giving a hard floor of 2.581 A. So we target just above the
# geometric floor. 3.331 -> ~2.6 A is still a large gain on a 0.20-weight term,
# and 2.6 A sits comfortably inside the Lauko window (2.0-3.6 A).
CST_TARGET_D_OXY = 2.40      # organisers' cst - unreachable with a fixed backbone
TARGET_D_OXY = 2.62          # achievable optimum for THIS backbone
TARGET_ANGLE2 = 120.0        # both cst files, and SCORING.md scores 120 -> 1.000
CLASH_MIN = 2.90             # heavy-atom clash cutoff for NON-partners

# Gln side-chain torsion chain. Each entry: (rotated_atoms, axis_from, axis_to)
CHI = [
    ("chi1", ("CG", "CD", "OE1", "NE2"), "CA", "CB"),
    ("chi2", ("CD", "OE1", "NE2"), "CB", "CG"),
    ("chi3", ("OE1", "NE2"), "CG", "CD"),
]


def rotate_about(points: np.ndarray, p0: np.ndarray, p1: np.ndarray,
                 deg: float) -> np.ndarray:
    axis = p1 - p0
    axis = axis / np.linalg.norm(axis)
    th = np.radians(deg)
    c, s = np.cos(th), np.sin(th)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    R = np.eye(3) + s * K + (1 - c) * (K @ K)
    return (points - p0) @ R.T + p0


def angle_deg(a, b, c) -> float:
    v1, v2 = a - b, c - b
    cos = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return float(np.degrees(np.arccos(max(-1.0, min(1.0, cos)))))


def main() -> int:
    tz = load(C.THEOZYME_WITH_ORI)
    gln = {a.name: a for a in tz.by_resname("GLN")}
    lig = {a.name: a for a in tz.by_resname(C.LIGAND_CCD)}
    asp = {a.name: a for a in tz.by_resname("ASP")}
    n1, n2 = lig["N1"].xyz, lig["N2"].xyz

    d0 = float(np.linalg.norm(gln["NE2"].xyz - n1))
    a0 = angle_deg(gln["NE2"].xyz, n1, n2)
    print(f"  supplied   d_oxy = {d0:.3f} A   angle2(NE2-N1-N2) = {a0:.2f} deg")
    print(f"  cst target d_oxy = {CST_TARGET_D_OXY:.2f} A  (UNREACHABLE: backbone floor 2.581 A)")
    print(f"  our target d_oxy = {TARGET_D_OXY:.2f} A   angle2 = {TARGET_ANGLE2:.0f} deg")

    # Atoms NE2 must not clash into. Critically this EXCLUDES the ligand N1,
    # which is the atom we are deliberately bringing to 2.4 A - counting the
    # intended H-bond partner as a clash makes the target unreachable by
    # construction (that was a real bug on the first attempt).
    movable = {"CG", "CD", "OE1", "NE2"}
    hbond_partners = {("6NT", "N1"), ("6NT", "N2")}
    static = np.array([a.xyz for a in tz.atoms
                       if not (a.resname == "GLN" and a.name in movable)
                       and a.element != "H"
                       and (a.resname, a.name) not in hbond_partners])

    coords = {k: v.xyz.copy() for k, v in gln.items()}
    best = None
    step = 5
    grid = range(-180, 180, step)
    for c1, c2, c3 in itertools.product(grid, grid, grid):
        cur = {k: v.copy() for k, v in coords.items()}
        for deg, (_, atoms, f, t) in zip((c1, c2, c3), CHI):
            if deg == 0:
                continue
            pts = np.array([cur[a] for a in atoms if a in cur])
            rot = rotate_about(pts, cur[f], cur[t], deg)
            for a, p in zip([a for a in atoms if a in cur], rot):
                cur[a] = p
        ne2 = cur["NE2"]
        d = float(np.linalg.norm(ne2 - n1))
        if abs(d - TARGET_D_OXY) > 0.25:
            continue
        ang = angle_deg(ne2, n1, n2)
        sc = np.array([cur[a] for a in ("CG", "CD", "OE1", "NE2")])
        dmat = np.linalg.norm(sc[:, None, :] - static[None, :, :], axis=-1)
        if dmat.min() < CLASH_MIN:
            continue
        cost = abs(d - TARGET_D_OXY) / 0.25 + abs(ang - TARGET_ANGLE2) / 40.0
        if best is None or cost < best[0]:
            best = (cost, (c1, c2, c3), d, ang, {k: v.copy() for k, v in cur.items()},
                    float(dmat.min()))

    if best is None:
        print("  NO clash-free rotamer reaches the 2.4 A target.")
        print("  Reporting failure rather than shipping a clashing motif.")
        return 1

    cost, chis, d, ang, newc, minclash = best
    print(f"  optimised  d_oxy = {d:.3f} A   angle2 = {ang:.2f} deg   "
          f"chi=({chis[0]},{chis[1]},{chis[2]})  min heavy-atom sep {minclash:.2f} A")
    print(f"  improvement: d_oxy {d0:.3f} -> {d:.3f} "
          f"(target {TARGET_D_OXY}), angle2 {a0:.2f} -> {ang:.2f} (target {TARGET_ANGLE2})")

    out = C.THEOZYMES / "Theozyme_withORI_oxyopt.pdb"
    lines = [f"{FRAME_REMARK} withORI",
             "REMARK 100 KEMP_OXYOPT Gln side chain re-optimised toward the "
             "organisers' cst target",
             f"REMARK 100 KEMP_OXYOPT d_oxy {d0:.3f} -> {d:.3f} A "
             f"(target {TARGET_D_OXY})",
             f"REMARK 100 KEMP_OXYOPT angle2 {a0:.2f} -> {ang:.2f} deg "
             f"(target {TARGET_ANGLE2})",
             f"REMARK 100 KEMP_OXYOPT chi1,chi2,chi3 = {chis[0]},{chis[1]},{chis[2]} deg",
             "REMARK 100 KEMP_OXYOPT backbone unchanged; sidechain only"]
    if tz.cryst1:
        lines.append(tz.cryst1)
    for a in tz.atoms:
        if a.resname == "ORI" or a.name == "ORI":
            continue
        raw = a.raw
        if a.resname == "GLN" and a.name in newc:
            x, y, z = newc[a.name]
            raw = f"{raw[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{raw[54:]}"
        # B-factor 0.00 on Gln is THE trap: APNovo infers motif residues from
        # atom_b_factor > 0, so a zero-B Gln is silently dropped, taking oxy
        # (0.20) and orient_oxy (0.08) with it. Stamp a non-zero B.
        if a.resname in ("GLN", "ASP", C.LIGAND_CCD):
            b = max(float(raw[60:66] or 0.0), 1.00)
            raw = f"{raw[:60]}{b:6.2f}{raw[66:]}"
        lines.append(raw[:16] + " " + raw[17:])
    lines.append("END")
    out.write_text("\n".join(lines) + "\n")
    print(f"  wrote {out}")

    tz2 = load(out)
    from kemp import official
    raw2 = official.measure_raw(tz2)
    print(f"  verify: d_base={raw2.d_base:.3f} d_oxy={raw2.d_oxy:.3f} "
          f"angle2_oxy={raw2.angle2_oxy:.2f}")
    bz = {a.bfactor for a in tz2.by_resname("GLN")}
    print(f"  Gln B-factors now: {sorted(bz)[:4]} (must be > 0 or APNovo drops it)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
