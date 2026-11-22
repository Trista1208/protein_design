"""Build a Rosetta .params for the 6NT transition-state analogue.

Rosetta cannot score a non-canonical ligand without a params file, so this is
the prerequisite for ANY Rosetta work on our complexes - FastDesign, relax,
enzdes constraints, FuncLib.

CONNECTIVITY is measured from the folded complex, not assumed. All 13 bonds
fall in 1.26-1.44 A and reproduce two fused rings:

    nitro     NO1-O21 1.287   NO1-O11 1.261   C5-NO1 1.440
    benzene   C4-C5 1.379  C5-C6 1.391  C6-C7 1.393
              C7-C7A 1.386  C7A-C3A 1.387  C3A-C4 1.397
    triazole  C7A-N1 1.367  N1-N2 1.330  N2-N3 1.347  N3-C3A 1.360

PROTONATION - a real chemical decision, not a default
-----------------------------------------------------
6NT is the transition-state ANALOGUE. In the real reaction the Asp abstracts a
proton from the substrate's C3. The analogue mimics that with a ring N-H
pointing at the same carboxylate, so the NH must sit on the nitrogen that
contacts the Asp.

Measured in our complexes: Asp OD2...N3 = 2.48 A (the catalytic contact) and
Gln NE2...N1 = 3.33 A (the oxyanion contact). So:

    N3-H ... OD2(Asp)    mimics the C3-H being abstracted   -> H on N3
    Gln NE2-H ... N1     Gln DONATES to N1                  -> N1 has no H

Kekule structure follows from that: N3 is the NH (single bonds to C3A and N2),
N1=N2 double, N1-C7A single, benzene alternating. Every atom closes its
valence and the formula is C6H4N4O2 - the organisers' formula.

The nitro is written charge-separated (N+ / O-) rather than hypervalent, which
is what molfile_to_params expects.

ATOM NAMES ARE PRESERVED. The enzdes .cst refers to ligand atoms N1, N2 and N3
by name; if the params renamed them the constraints would silently attach to
the wrong atoms, which is exactly the class of error that produces a confident
wrong answer.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC_PDB = ROOT / "outputs/refold_d045/seed2/folded.pdb"
OUT = ROOT / "outputs/rosetta"
MOLFILE_TO_PARAMS = Path("/home/mandalorian/opt/rosetta/source/scripts/python/"
                         "public/molfile_to_params.py")

#: PDB order is preserved so the params atom order is the molfile atom order.
ORDER = ["O21", "NO1", "O11", "C5", "C6", "C7", "C7A", "N1", "C3A", "C4",
         "N3", "N2"]

ELEMENT = {"O21": "O", "NO1": "N", "O11": "O", "C5": "C", "C6": "C", "C7": "C",
           "C7A": "C", "N1": "N", "C3A": "C", "C4": "C", "N3": "N", "N2": "N"}

#: (atom_a, atom_b, order). Ring bonds are order 4 = AROMATIC, because
#: molfile_to_params says so explicitly: "Aromatic bonds must be identified
#: explicitly -- alternating single/double bonds (Kekule structure) won't cut
#: it." A Kekule ring here produced a ligand Rosetta would not treat as planar.
BONDS = [
    # nitro - genuinely not aromatic
    ("NO1", "O21", 2), ("NO1", "O11", 1), ("C5", "NO1", 1),
    # benzene
    ("C4", "C5", 4), ("C5", "C6", 4), ("C6", "C7", 4),
    ("C7", "C7A", 4), ("C7A", "C3A", 4), ("C3A", "C4", 4),
    # triazole
    ("C7A", "N1", 4), ("N1", "N2", 4), ("N2", "N3", 4), ("N3", "C3A", 4),
]

#: charge-separated nitro
CHARGES = {"NO1": 1, "O11": -1}

#: Hydrogens, added explicitly. molfile_to_params does NOT add them, and
#: without the N3-H there is no donor for the Asp contact - the whole point of
#: the TSA. (heavy atom, its two ring neighbours, name, bond length)
HYDROGENS = [
    ("C4", ("C3A", "C5"), "H4", 1.08),
    ("C6", ("C5", "C7"), "H6", 1.08),
    ("C7", ("C6", "C7A"), "H7", 1.08),
    ("N3", ("N2", "C3A"), "HN3", 1.01),      # the proton the Asp abstracts
]


def read_coords() -> dict[str, tuple[float, float, float]]:
    xyz = {}
    for ln in SRC_PDB.read_text().splitlines():
        if ln.startswith("HETATM") and ln[17:20].strip() == "6NT":
            name = ln[12:16].strip()
            xyz[name] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
    return xyz


def place_hydrogens(xyz: dict) -> tuple[list[str], dict, list[tuple]]:
    """In-plane H positions, bisecting the two ring neighbours and pointing out."""
    import numpy as np
    names, extra_bonds = [], []
    for heavy, (n1, n2), hname, blen in HYDROGENS:
        p = np.array(xyz[heavy], float)
        a = np.array(xyz[n1], float)
        b = np.array(xyz[n2], float)
        # outward direction: away from both neighbours, in the ring plane
        v = (p - a) / np.linalg.norm(p - a) + (p - b) / np.linalg.norm(p - b)
        v /= np.linalg.norm(v)
        xyz[hname] = tuple(p + blen * v)
        ELEMENT[hname] = "H"
        names.append(hname)
        extra_bonds.append((heavy, hname, 1))
    return names, xyz, extra_bonds


def write_molfile(xyz: dict, path: Path) -> tuple[list[str], list[tuple]]:
    hnames, xyz, hbonds = place_hydrogens(xyz)
    order = ORDER + hnames
    bonds = BONDS + hbonds
    idx = {n: i + 1 for i, n in enumerate(order)}
    L = ["6NT", "  kemp pipeline - scripts/make_6nt_params.py", ""]
    L.append(f"{len(order):3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000")
    for n in order:
        x, y, z = xyz[n]
        L.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {ELEMENT[n]:<3} 0  0  0  0  0"
                 "  0  0  0  0  0  0")
    for a, b, o in bonds:
        L.append(f"{idx[a]:3d}{idx[b]:3d}{o:3d}  0  0  0  0")
    for n, c in CHARGES.items():
        L.append(f"M  CHG  1 {idx[n]:3d} {c:3d}")
    L += ["M  END", "$$$$", ""]
    path.write_text("\n".join(L))
    return order, bonds


def rename_to_our_atoms(params: Path, lig_pdb: Path, xyz: dict) -> dict:
    """Rewrite params/PDB atom names to OURS, matched by coordinate.

    molfile_to_params reports "Atom names contain duplications -- renaming all
    atoms" (a .mol file carries no names at all) and emits its own scheme. That
    scheme COLLIDES with ours in the worst possible way: its `N1` is our nitro
    nitrogen NO1, and its `N2` is our triazole N1. An enzdes .cst naming `N1`
    would then restrain the NITRO group instead of the catalytic ring nitrogen
    and report a perfectly plausible energy for the wrong geometry.

    Matching on coordinates is exact here (molfile_to_params preserves input
    coordinates in its companion PDB), so the mapping is unambiguous.
    """
    import numpy as np
    theirs = {}
    for ln in lig_pdb.read_text().splitlines():
        if ln.startswith(("HETATM", "ATOM  ")):
            theirs[ln[12:16].strip()] = np.array(
                [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    mapping, used = {}, set()
    for tn, tx in theirs.items():
        best = min((k for k in xyz if k not in used),
                   key=lambda k: np.linalg.norm(np.array(xyz[k], float) - tx))
        d = float(np.linalg.norm(np.array(xyz[best], float) - tx))
        if d > 0.05:
            raise SystemExit(f"atom {tn} matched {best} at {d:.3f} A - too far; "
                             "coordinate mapping is unreliable, refusing to rename")
        mapping[tn] = best
        used.add(best)

    def swap(text: str) -> str:
        """Rename WITHOUT disturbing column widths.

        Rosetta's .params reader takes atom names from fixed-width fields. A
        split()/join() rename collapses `ATOM  C6  aroC  X   0.22` into
        `ATOM C6 aroC X 0.22`, after which Rosetta reports
        "add_bond: atom C5 and/or NO1 don't exist!" - the ATOM lines no longer
        register. So substitute each name PADDED TO 4 CHARACTERS for the new
        name padded to 4, which is length-preserving by construction.
        """
        # SIMULTANEOUS, via a single regex pass. Substituting the pairs one
        # after another renames atoms that an earlier rule already renamed:
        # their C6 -> our C4, then the rule for their C4 -> our C7A fires on
        # the C4 just written. That produced a params with `C4` and `HN3`
        # appearing twice and two atoms missing entirely. re.sub does not
        # re-examine replaced text, so one pass is order-independent.
        import re
        pad = {f"{o:<4}": f"{n:<4}" for o, n in mapping.items()}
        pat = re.compile("|".join(re.escape(k) for k in
                                  sorted(pad, key=len, reverse=True)))
        out = []
        for ln in text.splitlines():
            if ln.startswith(("ATOM ", "BOND", "CHI", "NBR_ATOM",
                              "ICOOR_INTERNAL", "FIRST_SIDECHAIN_ATOM")):
                ln = pat.sub(lambda m: pad[m.group(0)], ln)
            out.append(ln)
        return "\n".join(out) + "\n"

    params.write_text(swap(params.read_text()))
    return mapping


def main() -> int:
    if not SRC_PDB.exists():
        print(f"missing {SRC_PDB}")
        return 2
    xyz = read_coords()
    missing = [n for n in ORDER if n not in xyz]
    if missing:
        print(f"atoms missing from the PDB: {missing}")
        return 2
    OUT.mkdir(parents=True, exist_ok=True)
    mol = OUT / "6NT.mol"
    order, bonds = write_molfile(xyz, mol)
    print(f"wrote {mol}  ({len(order)} atoms incl. "
          f"{len(HYDROGENS)} H, {len(bonds)} bonds)")

    # molfile_to_params is python2-era but runs under python3 in recent Rosetta.
    cmd = [sys.executable, str(MOLFILE_TO_PARAMS), "-n", "6NT", "-p", "6NT",
           "--keep-names", "--clobber", str(mol)]
    print("  " + " ".join(cmd))
    p = subprocess.run(cmd, cwd=str(OUT), capture_output=True, text=True,
                       timeout=600)
    print(p.stdout[-2500:] or "(no stdout)")
    if p.returncode != 0:
        print("STDERR:\n" + p.stderr[-2500:])
        return 3

    params = OUT / "6NT.params"
    if not params.exists():
        print("params not produced")
        return 3
    lig_pdb = OUT / "6NT_0001.pdb"
    mapping = rename_to_our_atoms(params, lig_pdb, xyz)
    print("\nrenamed to our atom names (matched by coordinate):")
    for tn, on in sorted(mapping.items()):
        print(f"   {tn:5} -> {on}" + ("   [collided]" if tn != on and
                                      tn in ("N1", "N2", "N3") else ""))

    text = params.read_text()
    names = [ln.split()[1] for ln in text.splitlines() if ln.startswith("ATOM ")]
    print(f"\n{params}")
    print(f"  atoms: {names}")

    # VERIFY IDENTITY, not mere presence. The previous version of this check
    # asked "is a thing called N1 in the file?" and passed while N1 was the
    # nitro nitrogen.
    ok = True
    bonded = {}
    for ln in text.splitlines():
        if ln.startswith("BOND"):
            t = ln.split()
            bonded.setdefault(t[1], set()).add(t[2])
            bonded.setdefault(t[2], set()).add(t[1])
    expect = {"N1": {"C7A", "N2"}, "N2": {"N1", "N3"},
              "N3": {"N2", "C3A", "HN3"}}
    for n, want in expect.items():
        got = bonded.get(n, set())
        good = got == want
        ok &= good
        print(f"  {n}: bonded to {sorted(got)}  expected {sorted(want)}"
              f"   {'OK' if good else '<-- WRONG ATOM'}")
    print("\n  catalytic contacts the .cst will make:")
    print("    Asp OD2 ... N3   (N3 carries HN3, the proton the Asp abstracts)")
    print("    Gln NE2 ... N1   (N1 accepts from the Gln amide)")
    return 0 if ok else 3


if __name__ == "__main__":
    raise SystemExit(main())
