"""Tests for kemp.strain - the post-RFD3 structural validation.

Written against the failure modes actually observed in the 2026-11-21 outputs,
not against invented ones:

* RFD3 reported ``n_chainbreaks: 0`` for a design containing a 2.68 A C-N bond.
  ``test_cut_detected_where_rfd3_reported_none`` is that case.
* Every RFD3 design had distorted peptide bonds flanking BOTH motif residues
  while the rest of the chain was within 0.141 A of ideal. The junction split
  has to survive refactoring, hence ``test_junction_split``.
* d_base/d_oxy are unmeasurable without the ligand chain. That must not read as
  a pass - ``test_missing_ligand_is_not_a_pass``.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from kemp import strain

IDEAL_CN = 1.336


# --------------------------------------------------------------------------- helpers
def _pdb(lines: list[str]) -> str:
    return "\n".join(lines) + "\nEND\n"


def _atom(serial: int, name: str, resname: str, resseq: int, xyz, chain="A",
          element=" C"):
    x, y, z = xyz
    return (f"ATOM  {serial:5d} {name:<4s}{resname:>4s} {chain}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element}")


def _glycine_chain(n: int, *, cn: float = IDEAL_CN, break_after: int | None = None,
                   break_cn: float = 4.0) -> str:
    """A straight poly-Gly chain with controllable C-N spacing.

    Geometry is not physical beyond the C-N distances - these tests are about
    the bond-length bookkeeping, and every other term is left unchecked by
    using GLY (no CB) and asserting only on integrity output.
    """
    lines, serial, x = [], 1, 0.0
    for i in range(1, n + 1):
        gap = break_cn if (break_after is not None and i - 1 == break_after) else cn
        if i > 1:
            x += gap
        for name, dx in (("N", 0.0), ("CA", 1.458), ("C", 1.458 + 1.525)):
            lines.append(_atom(serial, name, "GLY", i, (x + dx, 0.0, 0.0)))
            serial += 1
        x += 1.458 + 1.525
    return _pdb(lines)


# --------------------------------------------------------------------------- parsing
def test_parses_pdb_and_counts_bonds(tmp_path: Path):
    p = tmp_path / "chain.pdb"
    p.write_text(_glycine_chain(5))
    atoms = strain.load_atoms(p)
    assert len({(a.chain, a.resseq) for a in atoms}) == 5
    integ = strain.check_chain_integrity(atoms)
    assert integ["n_peptide_bonds"] == 4
    assert integ["n_cuts"] == 0
    assert integ["cn_mean_A"] == pytest.approx(IDEAL_CN, abs=1e-3)


def test_cif_parsed_by_column_name_not_position(tmp_path: Path):
    """RFD3 writes 24 atom_site columns in a non-PDB order.

    Reading by index would pick up B_iso as a coordinate. The column order here
    is deliberately scrambled relative to the PDB convention.
    """
    cif = """data_test
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_seq_id
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM N N . GLY A 1 99.900 1.000 2.000 3.000
ATOM C CA . GLY A 1 99.900 4.000 5.000 6.000
"""
    p = tmp_path / "t.cif"
    p.write_text(cif)
    atoms = strain.load_atoms(p)
    assert len(atoms) == 2
    assert atoms[0].xyz == (1.0, 2.0, 3.0)      # not 99.9
    assert atoms[1].name == "CA"
    assert atoms[1].resname == "GLY"


# --------------------------------------------------------------------------- 1. integrity
def test_cut_detected_where_rfd3_reported_none(tmp_path: Path):
    """The observed case: a 2.68 A C-N bond with n_chainbreaks reported as 0."""
    p = tmp_path / "cut.pdb"
    p.write_text(_glycine_chain(6, break_after=3, break_cn=2.681))
    integ = strain.check_chain_integrity(strain.load_atoms(p))
    assert integ["n_cuts"] == 1
    assert integ["cuts"][0]["C_N_A"] == pytest.approx(2.681, abs=1e-2)
    assert integ["n_nonphysical"] >= 1


def test_compressed_bond_is_nonphysical(tmp_path: Path):
    """0.911 A was observed in diag_8seq_fixed. Too short is as wrong as too long."""
    p = tmp_path / "squashed.pdb"
    p.write_text(_glycine_chain(6, break_after=2, break_cn=0.911))
    integ = strain.check_chain_integrity(strain.load_atoms(p))
    assert integ["n_cuts"] == 0           # not long, so not a "cut"
    assert integ["n_nonphysical"] == 1    # but still impossible


def test_junction_split(tmp_path: Path):
    """Distortion at the motif must not be averaged into the whole-chain stats."""
    p = tmp_path / "j.pdb"
    p.write_text(_glycine_chain(8, break_after=4, break_cn=1.80 - 0.01))
    atoms = strain.load_atoms(p)
    # residue 5 is the motif: bonds 4-5 and 5-6 are its junctions
    integ = strain.check_chain_integrity(atoms, motif_resseq={5})
    assert integ["junction_max_dev_A"] > strain.CN_DEV_DISTORTED
    assert integ["nonjunction_max_dev_A"] == pytest.approx(0.0, abs=1e-3)
    assert [b["junction"] for b in integ["junction_bonds"]] == [True, True]


def test_numbering_gap_is_not_a_bond(tmp_path: Path):
    """A gap in residue numbering is not a measurable bond, and not a break."""
    lines = []
    serial = 1
    for i, resseq in enumerate((1, 2, 10, 11)):
        x = i * 3.8
        for name, dx in (("N", 0.0), ("CA", 1.4), ("C", 2.9)):
            lines.append(_atom(serial, name, "GLY", resseq, (x + dx, 0.0, 0.0)))
            serial += 1
    p = tmp_path / "gap.pdb"
    p.write_text(_pdb(lines))
    integ = strain.check_chain_integrity(strain.load_atoms(p))
    assert integ["n_peptide_bonds"] == 2   # 1-2 and 10-11, never 2-10
    assert integ["n_cuts"] == 0


# --------------------------------------------------------------------------- geometry primitives
def test_angle_against_known_value():
    o = strain.Atom("A", 1, "GLY", "O", (0.0, 0.0, 0.0))
    x = strain.Atom("A", 1, "GLY", "X", (1.0, 0.0, 0.0))
    y = strain.Atom("A", 1, "GLY", "Y", (0.0, 1.0, 0.0))
    assert strain.angle(x, o, y) == pytest.approx(90.0, abs=1e-6)


@pytest.mark.parametrize("d_xyz,expected", [
    ((1.0, 1.0, 0.0), 0.0),       # d on the same side as a   -> cis / syn
    ((1.0, -1.0, 0.0), 180.0),    # d on the opposite side    -> trans / anti
    ((1.0, 0.0, 1.0), -90.0),     # +z is anticlockwise from +y looking down bc
    ((1.0, 0.0, -1.0), 90.0),
])
def test_dihedral_against_known_values(d_xyz, expected):
    """a at +y, b at the origin, c along +x: the dihedral is then read directly
    off where d sits relative to the bc axis, with no trigonometry to get wrong.
    """
    a = strain.Atom("A", 1, "G", "A", (0.0, 1.0, 0.0))
    b = strain.Atom("A", 1, "G", "B", (0.0, 0.0, 0.0))
    c = strain.Atom("A", 1, "G", "C", (1.0, 0.0, 0.0))
    d = strain.Atom("A", 1, "G", "D", d_xyz)
    got = strain.dihedral(a, b, c, d)
    if expected == 180.0:
        assert abs(got) == pytest.approx(180.0, abs=1e-6)
    else:
        assert got == pytest.approx(expected, abs=1e-6)


# --------------------------------------------------------------------------- 2. strain
# --- ideal-geometry construction -------------------------------------------
# Built from internal coordinates rather than typed-in Cartesians, because
# hand-placed "ideal" atoms are not ideal: the first version of this fixture put
# the carboxylate oxygens at 1.345 A instead of 1.249 and the test failed on
# the fixture, not the code.
import numpy as np


def _unit(v):
    return np.asarray(v, dtype=float) / np.linalg.norm(v)


def _two_angle_solve(n_hat, c_hat, ang_n, ang_c):
    """Unit vector at ``ang_n`` from ``n_hat`` and ``ang_c`` from ``c_hat``.

    Used for CB, which is pinned by both N-CA-CB and CB-CA-C.
    """
    A = np.vstack([n_hat, c_hat])
    b = np.array([math.cos(math.radians(ang_n)), math.cos(math.radians(ang_c))])
    v0 = np.linalg.pinv(A) @ b                 # in-plane part
    perp = _unit(np.cross(n_hat, c_hat))
    rem = 1.0 - float(v0 @ v0)
    assert rem > 0, "angles are geometrically impossible"
    return v0 + math.sqrt(rem) * perp


def _perp_to(a_hat):
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(float(a_hat @ tmp)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    return _unit(np.cross(a_hat, tmp))


def _place_one(c_prev, c, bond, angle_deg):
    """Atom at ``bond`` from ``c``, making ``angle_deg`` with c->c_prev.

    Works for any orientation, unlike a rotation about a fixed axis. The
    dihedral is left arbitrary because no test asserts on one.
    """
    a_hat = _unit(np.asarray(c_prev) - np.asarray(c))
    p = _perp_to(a_hat)
    t = math.radians(angle_deg)
    return np.asarray(c, dtype=float) + bond * (math.cos(t) * a_hat + math.sin(t) * p)


def _planar_pair(c_prev, c, b1, a1, b2, a2):
    """Two substituents on ``c``, coplanar, on opposite sides of c->c_prev.

    With a1 + a2 + angle(x1, c, x2) = 360 the carboxylate and amide groups are
    planar by construction, which is what the Engh & Huber angle triples encode.
    """
    a_hat = _unit(np.asarray(c_prev) - np.asarray(c))
    p = _perp_to(a_hat)
    x1 = c + b1 * (math.cos(math.radians(a1)) * a_hat + math.sin(math.radians(a1)) * p)
    x2 = c + b2 * (math.cos(math.radians(a2)) * a_hat - math.sin(math.radians(a2)) * p)
    return x1, x2


def _ideal_residue(resname: str, *, ca_cb: float = 1.530) -> dict[str, np.ndarray]:
    """Asp or Gln with every term in IDEAL_BOND / IDEAL_ANGLE at its ideal value.

    ``ca_cb`` is exposed so a test can perturb exactly one bond.
    """
    I = strain.IDEAL_BOND
    A = strain.IDEAL_ANGLE
    n = np.zeros(3)
    ca = np.array([I[("*", "N", "CA")], 0.0, 0.0])
    n_hat = _unit(n - ca)
    c = _place_one(n, ca, I[("*", "CA", "C")], A[("*", "N", "CA", "C")])
    o = _place_one(ca, c, I[("*", "C", "O")], A[("*", "CA", "C", "O")])
    cb_dir = _two_angle_solve(n_hat, _unit(c - ca),
                              A[("*", "N", "CA", "CB")], A[("*", "CB", "CA", "C")])
    cb = ca + ca_cb * _unit(cb_dir)
    out = {"N": n, "CA": ca, "C": c, "O": o, "CB": cb}

    if resname == "ASP":
        cg = _place_one(ca, cb, I[("ASP", "CB", "CG")],
                        A[("ASP", "CA", "CB", "CG")])
        od1, od2 = _planar_pair(cb, cg,
                                I[("ASP", "CG", "OD1")], A[("ASP", "CB", "CG", "OD1")],
                                I[("ASP", "CG", "OD2")], A[("ASP", "CB", "CG", "OD2")])
        out.update({"CG": cg, "OD1": od1, "OD2": od2})
    elif resname == "GLN":
        cg = _place_one(ca, cb, I[("GLN", "CB", "CG")],
                        A[("GLN", "CA", "CB", "CG")])
        cd = _place_one(cb, cg, I[("GLN", "CG", "CD")],
                        A[("GLN", "CB", "CG", "CD")])
        oe1, ne2 = _planar_pair(cg, cd,
                                I[("GLN", "CD", "OE1")], A[("GLN", "CG", "CD", "OE1")],
                                I[("GLN", "CD", "NE2")], A[("GLN", "CG", "CD", "NE2")])
        out.update({"CG": cg, "CD": cd, "OE1": oe1, "NE2": ne2})
    else:
        raise ValueError(resname)
    return out


def _asp_gln_pdb(*, asp_ca_cb: float = 1.530) -> str:
    """Ideal Asp + Gln, placed 20 A apart so neither perturbs the other."""
    lines, serial = [], 1
    for resname, resseq, shift, ca_cb in (("ASP", 10, 0.0, asp_ca_cb),
                                          ("GLN", 20, 20.0, 1.530)):
        for name, xyz in _ideal_residue(resname, ca_cb=ca_cb).items():
            x, y, z = (float(v) for v in xyz)
            lines.append(_atom(serial, name, resname, resseq, (x, y + shift, z)))
            serial += 1
    return _pdb(lines)


def test_fixture_is_actually_ideal():
    """Guard the guard: if this fixture drifts, every strain test below lies."""
    for resname in ("ASP", "GLN"):
        r = _ideal_residue(resname)
        atoms = {k: strain.Atom("A", 1, resname, k, tuple(v)) for k, v in r.items()}
        for (tag, a, b), ideal in strain.IDEAL_BOND.items():
            if tag in ("*", resname) and a in atoms and b in atoms:
                assert strain.dist(atoms[a], atoms[b]) == pytest.approx(ideal, abs=1e-6), \
                    f"{resname} bond {a}-{b}"
        for (tag, a, b, c), ideal in strain.IDEAL_ANGLE.items():
            if tag in ("*", resname) and a in atoms and b in atoms and c in atoms:
                assert strain.angle(atoms[a], atoms[b], atoms[c]) == pytest.approx(ideal, abs=1e-4), \
                    f"{resname} angle {a}-{b}-{c}"


def test_ideal_sidechain_is_not_strained(tmp_path: Path):
    p = tmp_path / "ideal.pdb"
    p.write_text(_asp_gln_pdb())
    s = strain.check_catalytic_strain(strain.load_atoms(p), base_resseq=10,
                                      oxy_resseq=20)
    base = s["residues"]["base"]
    assert base["present"] and base["resname"] == "ASP"
    bonds = [t for t in base["terms"] if t["kind"] == "bond"]
    assert bonds, "no bond terms measured"
    assert all(not t["strained"] for t in bonds), [t for t in bonds if t["strained"]]


def test_stretched_ca_cb_is_strained(tmp_path: Path):
    """GLN13 in d00002 had CA-CB = 1.746 A against an ideal of 1.530."""
    p = tmp_path / "strained.pdb"
    p.write_text(_asp_gln_pdb(asp_ca_cb=1.746))
    s = strain.check_catalytic_strain(strain.load_atoms(p), base_resseq=10,
                                      oxy_resseq=20)
    base = s["residues"]["base"]
    hit = [t for t in base["terms"] if t["atoms"] == "CA-CB"]
    assert hit and hit[0]["strained"]
    assert hit[0]["value"] == pytest.approx(1.746, abs=1e-3)
    assert s["n_strained_total"] >= 1


def test_missing_residue_reported_not_crashed(tmp_path: Path):
    p = tmp_path / "ideal.pdb"
    p.write_text(_asp_gln_pdb())
    s = strain.check_catalytic_strain(strain.load_atoms(p), base_resseq=10,
                                      oxy_resseq=999)
    assert s["residues"]["oxy"]["present"] is False


# --------------------------------------------------------------------------- 3. distances
def test_missing_ligand_is_not_a_pass(tmp_path: Path):
    """The campaign's RFD3 outputs had no ligand. Unmeasurable != fine."""
    p = tmp_path / "noligand.pdb"
    p.write_text(_asp_gln_pdb())
    d = strain.check_catalytic_distances(strain.load_atoms(p), base_resseq=10,
                                         oxy_resseq=20)
    assert d["ligand_present"] is False
    assert d["d_base"]["measurable"] is False
    assert d["d_base"]["within_tol"] is False      # NOT silently True
    assert d["d_oxy"]["within_tol"] is False

    v = strain.validate(p, base_resseq=10, oxy_resseq=20)
    assert not v.passed
    assert "d_base_unmeasurable" in v.failed
    assert "d_oxy_unmeasurable" in v.failed


def test_catalytic_distances_measured_against_reference(tmp_path: Path):
    """Place a ligand at exactly the reference distances and expect a match."""
    from kemp import config as C
    ref = C.REF_GEOMETRY
    lines, serial = [], 1
    # Asp OD2 at origin, Gln NE2 on the x axis far away
    lines.append(_atom(serial, "CA", "ASP", 10, (0.0, 0.0, 5.0))); serial += 1
    lines.append(_atom(serial, "OD2", "ASP", 10, (0.0, 0.0, 0.0))); serial += 1
    lines.append(_atom(serial, "CA", "GLN", 20, (ref["asp127_ca__gln50_ca"], 0.0, 5.0))); serial += 1
    # 6NT N3 at d_base from OD2; N1 at d_oxy from NE2
    n3 = (ref["asp127_od2__6nt_n3"], 0.0, 0.0)
    lines.append(_atom(serial, "N3", "6NT", 1, n3, chain="L", element=" N")); serial += 1
    n1 = (ref["asp127_od2__6nt_n3"] + 1.4, 0.0, 0.0)
    lines.append(_atom(serial, "N1", "6NT", 1, n1, chain="L", element=" N")); serial += 1
    ne2 = (n1[0] + ref["gln50_ne2__6nt_n1"], 0.0, 0.0)
    lines.append(_atom(serial, "NE2", "GLN", 20, ne2)); serial += 1
    p = tmp_path / "withlig.pdb"
    p.write_text(_pdb(lines))

    d = strain.check_catalytic_distances(strain.load_atoms(p), base_resseq=10,
                                         oxy_resseq=20)
    assert d["ligand_present"] is True
    assert d["d_base"]["value"] == pytest.approx(ref["asp127_od2__6nt_n3"], abs=1e-3)
    assert d["d_base"]["within_tol"]
    assert d["d_oxy"]["value"] == pytest.approx(ref["gln50_ne2__6nt_n1"], abs=1e-3)
    assert d["d_oxy"]["within_tol"]
    assert d["ca_ca"]["within_tol"]


# --------------------------------------------------------------------------- gate wiring
def test_validation_as_dict_is_ledger_safe(tmp_path: Path):
    import json
    p = tmp_path / "chain.pdb"
    p.write_text(_asp_gln_pdb())
    v = strain.validate(p, base_resseq=10, oxy_resseq=20)
    json.dumps(v.as_dict())      # must not raise
    assert set(v.as_dict()) >= {"passed", "failed", "junction_max_dev_A",
                                "n_cuts", "d_base_A", "d_oxy_A"}


def test_thresholds_are_overridable(tmp_path: Path):
    """Accepting a design despite an unmeasurable ligand must be explicit."""
    p = tmp_path / "noligand.pdb"
    p.write_text(_asp_gln_pdb())
    v = strain.validate(p, base_resseq=10, oxy_resseq=20,
                        thresholds={"require_catalytic_measurable": False,
                                    "require_catalytic_within_tol": False})
    assert v.passed, v.failed


# --------------------------------------------------------------------------- regression
#: phi/psi of the catalytic residues in the REFERENCE crystal structures, read
#: off the deposited coordinates. An earlier version of ramachandran_region()
#: classified every one of these as "disallowed" - i.e. it condemned the very
#: enzyme the challenge is benchmarked against. These values are the ground
#: truth for what this site's backbone geometry looks like, so the classifier
#: must never call them impossible again.
CRYSTAL_CATALYTIC_DIHEDRALS = [
    ("4BS0", "ASP", 87.94, -99.78),
    ("4BS0", "GLN", 70.63, -166.26),
    ("7K4Z", "ASP", 82.74, -101.05),
    ("7K4Z", "GLN", 71.65, -167.14),
    ("6TU6", "ASP", 87.07, -95.25),
    ("6TU6", "MET", 66.97, -164.58),
]


@pytest.mark.parametrize("pdb,resname,phi,psi", CRYSTAL_CATALYTIC_DIHEDRALS)
def test_reference_enzyme_conformations_are_allowed(pdb, resname, phi, psi):
    r = strain.ramachandran_region(resname, phi, psi)
    assert r["allowed"] is True, (
        f"{pdb} {resname} phi={phi} psi={psi} is a refined crystallographic "
        f"conformation but was classified {r['region']!r}")


def test_ramachandran_never_gates(tmp_path: Path):
    """phi/psi is reported, never a pass/fail.

    Behavioural rather than source-inspecting: a catalytic residue in the
    reference enzymes' own unusual conformation must not produce any
    Ramachandran failure or warning.
    """
    p = tmp_path / "ideal.pdb"
    p.write_text(_asp_gln_pdb())
    v = strain.validate(p, base_resseq=10, oxy_resseq=20,
                        thresholds={"require_catalytic_measurable": False,
                                    "require_catalytic_within_tol": False})
    blame = " ".join(v.failed + v.warnings).lower()
    assert "ramachandran" not in blame, blame
    # but the measurement is still recorded
    base = v.strain["residues"]["base"]
    assert "ramachandran" in base and "region" in base["ramachandran"]
