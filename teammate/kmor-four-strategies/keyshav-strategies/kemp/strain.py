"""Post-RFD3 structural validation: is the backbone sound before we pay for MPNN?

Three questions, asked directly of the coordinates rather than of RFD3's own
self-reported metrics:

1. **Is the chain unbroken?**  RFD3 reports ``n_chainbreaks`` and
   :func:`kemp.official.gate_backbone` trusts it. That number is a model output,
   not a measurement; this module measures every peptide C-N distance itself.

2. **Are the catalytic residues strain-free?**  Bond lengths, bond angles and
   backbone dihedrals of the Asp base and Gln oxyanion stabiliser, against
   Engh & Huber ideals. Nothing checked this before. It matters because a
   strained motif is a motif LigandMPNN and Boltz will quietly relax away: the
   geometry that scores is the geometry AFTER folding, so strain built in here
   shows up downstream as a ligand that drifts, not as an error.

3. **Are the catalytic distances intact?**  ``d_base`` (Asp OD2 - 6NT N3) and
   ``d_oxy`` (Gln NE2 - 6NT N1) measured on the design itself, plus the
   CA-CA scaffolding distance and the OD2-N3-N1 angle, against
   :data:`kemp.config.REF_GEOMETRY`.

WHY THIS IS A GATE AND NOT A REPORT: the official score is a geometric mean, so
a motif that arrives at Boltz already strained or displaced scores a hard 0.000
however good the fold is. Measuring here costs milliseconds; discovering it
after LigandMPNN + 8 x Boltz costs ~20 minutes of GPU.

NOTE ON THE LIGAND: ``d_base`` and ``d_oxy`` are only measurable when RFD3 was
run with the ligand on its own chain (chain ``L``, resname ``6NT``). Runs made
before that fix contain protein only, and this module reports those two checks
as ``None`` / unmeasurable rather than passing them. An unmeasurable catalytic
distance is NOT a pass - see :func:`validate`.
"""

from __future__ import annotations

import gzip
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import config as C

AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

LIGAND_RESNAMES = {"6NT", "LIG"}


# --------------------------------------------------------------------------- parsing
@dataclass(frozen=True)
class Atom:
    chain: str
    resseq: int
    resname: str
    name: str
    xyz: tuple[float, float, float]

    @property
    def v(self) -> np.ndarray:
        return np.asarray(self.xyz, dtype=float)


def load_atoms(path: Path | str) -> list[Atom]:
    """Atoms from a PDB or (optionally gzipped) mmCIF.

    Parsed here rather than via gemmi/biotite so that validation never depends
    on which conda env happens to be importable - this runs inline in the
    pipeline, in whatever interpreter the pipeline is already using.
    """
    path = Path(path)
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as fh:  # type: ignore[operator]
        text = fh.read()
    name = path.name[:-3] if path.name.endswith(".gz") else path.name
    if name.endswith((".cif", ".mmcif")):
        return _parse_cif(text)
    return _parse_pdb(text)


def _parse_pdb(text: str) -> list[Atom]:
    out = []
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        altloc = line[16]
        if altloc not in (" ", "", "A"):
            continue
        try:
            out.append(Atom(
                chain=line[21],
                resseq=int(line[22:26]),
                resname=line[17:20].strip(),
                name=line[12:16].strip(),
                xyz=(float(line[30:38]), float(line[38:46]), float(line[46:54])),
            ))
        except ValueError:
            continue
    return out


def _parse_cif(text: str) -> list[Atom]:
    """Parse the ``_atom_site`` loop by HEADER NAME, never by fixed position.

    RFD3 writes 24 atom_site columns in an order that is not the PDB's; reading
    by index would silently pick up ``B_iso`` as a coordinate.
    """
    lines = text.splitlines()
    cols: list[str] = []
    start = None
    for i, line in enumerate(lines):
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".", 1)[1].split()[0])
            start = i + 1
        elif cols and start is not None and i >= start:
            break
    if not cols:
        return []

    def idx(*names: str) -> int | None:
        for n in names:
            if n in cols:
                return cols.index(n)
        return None

    i_chain = idx("auth_asym_id", "label_asym_id")
    i_seq = idx("auth_seq_id", "label_seq_id")
    i_comp = idx("label_comp_id", "auth_comp_id")
    i_name = idx("label_atom_id", "auth_atom_id")
    i_alt = idx("label_alt_id")
    i_x, i_y, i_z = idx("Cartn_x"), idx("Cartn_y"), idx("Cartn_z")
    if None in (i_chain, i_seq, i_comp, i_name, i_x, i_y, i_z):
        return []

    out = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        f = line.split()
        if len(f) <= max(i_x, i_y, i_z):  # type: ignore[type-var]
            continue
        if i_alt is not None and f[i_alt] not in (".", "?", "A"):
            continue
        try:
            resseq = int(f[i_seq])  # type: ignore[index]
        except ValueError:
            continue
        try:
            xyz = (float(f[i_x]), float(f[i_y]), float(f[i_z]))  # type: ignore[index]
        except ValueError:
            continue
        out.append(Atom(chain=f[i_chain], resseq=resseq,  # type: ignore[index]
                        resname=f[i_comp], name=f[i_name], xyz=xyz))  # type: ignore[index]
    return out


def residues(atoms: list[Atom]) -> dict[tuple[str, int], dict[str, Atom]]:
    """{(chain, resseq): {atom_name: Atom}}, first occurrence wins."""
    out: dict[tuple[str, int], dict[str, Atom]] = {}
    for a in atoms:
        out.setdefault((a.chain, a.resseq), {}).setdefault(a.name, a)
    return out


def _resname(atoms: list[Atom], key: tuple[str, int]) -> str:
    for a in atoms:
        if (a.chain, a.resseq) == key:
            return a.resname
    return ""


# --------------------------------------------------------------------------- geometry primitives
def dist(a: Atom, b: Atom) -> float:
    return float(np.linalg.norm(a.v - b.v))


def angle(a: Atom, b: Atom, c: Atom) -> float:
    """Angle a-b-c in degrees."""
    u, v = a.v - b.v, c.v - b.v
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu == 0 or nv == 0:
        return float("nan")
    cos = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    return math.degrees(math.acos(cos))


def dihedral(a: Atom, b: Atom, c: Atom, d: Atom) -> float:
    """Dihedral a-b-c-d in degrees, in (-180, 180], IUPAC sign convention.

    NOTE the direction of ``b0``. Using ``a - b`` here instead of ``b - a``
    flips ``n1``, which offsets every result by exactly 180 deg - turning trans
    peptides into cis ones. The first version of this function had that bug and
    reported omega ~9 deg for ordinary trans peptides;
    ``test_dihedral_against_known_values`` pins the convention down.
    """
    b0, b1, b2 = b.v - a.v, c.v - b.v, d.v - c.v
    n1 = np.cross(b0, b1)
    n2 = np.cross(b1, b2)
    m = np.cross(n1, b1 / max(np.linalg.norm(b1), 1e-9))
    x = float(np.dot(n1, n2))
    y = float(np.dot(m, n2))
    return math.degrees(math.atan2(y, x))


# --------------------------------------------------------------------------- ideals
#: Engh & Huber (1991) protein-geometry ideals, as used by CCP4/REFMAC. Bond
#: lengths in Angstrom, angles in degrees. Only the atoms we actually check.
IDEAL_BOND = {
    ("*", "N", "CA"): 1.459,
    ("*", "CA", "C"): 1.525,
    ("*", "C", "O"): 1.229,
    ("*", "CA", "CB"): 1.530,
    ("peptide", "C", "N"): 1.336,
    ("ASP", "CB", "CG"): 1.516,
    ("ASP", "CG", "OD1"): 1.249,
    ("ASP", "CG", "OD2"): 1.249,
    ("GLN", "CB", "CG"): 1.520,
    ("GLN", "CG", "CD"): 1.516,
    ("GLN", "CD", "OE1"): 1.235,
    ("GLN", "CD", "NE2"): 1.333,
}

IDEAL_ANGLE = {
    ("*", "N", "CA", "C"): 111.0,
    ("*", "CA", "C", "O"): 120.8,
    ("*", "N", "CA", "CB"): 110.5,
    ("*", "CB", "CA", "C"): 110.1,
    ("ASP", "CA", "CB", "CG"): 112.6,
    ("ASP", "CB", "CG", "OD1"): 118.4,
    ("ASP", "CB", "CG", "OD2"): 118.4,
    ("ASP", "OD1", "CG", "OD2"): 123.2,
    ("GLN", "CA", "CB", "CG"): 114.1,
    ("GLN", "CB", "CG", "CD"): 112.6,
    ("GLN", "CG", "CD", "OE1"): 121.6,
    ("GLN", "CG", "CD", "NE2"): 116.7,
    ("GLN", "OE1", "CD", "NE2"): 121.7,
}

#: Peptide bond length above which the chain is cut, not merely stretched.
#: Ideal 1.336; a genuine break in a diffusion output is several Angstrom.
CHAINBREAK_CN_MAX = 1.80
#: Trans omega is 180; cis (Pro) is 0. Anything else is a twisted peptide.
OMEGA_TRANS_TOL = 30.0

#: MEASURED REFERENCE (scripts/audit_rfd3_strain.py, 2026-11-21).
#: Five Boltz-folded structures, 1000 peptide bonds: mean 1.327, sd 0.007,
#: min 1.302, max 1.370 - i.e. ZERO bonds more than 0.05 A off ideal. That is
#: what a physically valid chain measures with this code, and it is why the
#: numbers below can be trusted as artifacts rather than measurement error.
#:
#: Ten RFD3 outputs, same code: sd 0.03-0.25, with a clean separation
#:   away from the motif   |dev| <= 0.141 A   (worst case over all 10)
#:   flanking the motif    |dev|  0.15-2.94 A (ALL 10, BOTH junctions)
#: So RFD3's local geometry is sound everywhere except the four peptide bonds
#: that weld the rigid motif into the diffused backbone. Thresholds are set to
#: that measured separation, with the fatal tier placed well clear of it.
CN_DEV_DISTORTED = 0.15   # above the 0.141 A non-junction worst case
CN_DEV_NONPHYSICAL = 0.30  # a bond this far off cannot be relaxed, only rebuilt


# --------------------------------------------------------------------------- 1. chain integrity
def check_chain_integrity(atoms: list[Atom], *,
                          motif_resseq: set[int] | None = None) -> dict:
    """Measure every peptide C-N bond. Does not trust RFD3's own count.

    ``motif_resseq`` splits the result into junction bonds (those flanking a
    motif residue) and the rest. That split is the whole point: RFD3's backbone
    is sound away from the motif and distorted at it, so one summary number
    over all bonds hides the only signal that matters.
    """
    motif = motif_resseq or set()
    ideal = IDEAL_BOND[("peptide", "C", "N")]
    res = residues(atoms)
    keys = sorted((k for k in res if _resname(atoms, k) in AA3),
                  key=lambda k: (k[0], k[1]))

    bonds = []
    for (ch_a, seq_a), (ch_b, seq_b) in zip(keys, keys[1:]):
        if ch_a != ch_b or seq_b != seq_a + 1:
            continue  # a numbering gap is not a measurable bond
        c = res[(ch_a, seq_a)].get("C")
        n = res[(ch_b, seq_b)].get("N")
        if c is None or n is None:
            continue
        d = dist(c, n)
        bonds.append({
            "between": [f"{ch_a}{seq_a}", f"{ch_b}{seq_b}"],
            "C_N_A": round(d, 3),
            "dev": round(d - ideal, 3),
            "junction": bool(seq_a in motif or seq_b in motif),
            "cut": bool(d > CHAINBREAK_CN_MAX),
            "nonphysical": bool(abs(d - ideal) > CN_DEV_NONPHYSICAL),
            "distorted": bool(abs(d - ideal) > CN_DEV_DISTORTED),
        })

    def devs(sel):
        return [abs(b["dev"]) for b in bonds if sel(b)]

    jd = devs(lambda b: b["junction"])
    nd = devs(lambda b: not b["junction"])
    return {
        "n_peptide_bonds": len(bonds),
        "n_cuts": sum(1 for b in bonds if b["cut"]),
        "n_nonphysical": sum(1 for b in bonds if b["nonphysical"]),
        "n_distorted": sum(1 for b in bonds if b["distorted"]),
        "cuts": [b for b in bonds if b["cut"]][:10],
        "nonphysical": [b for b in bonds if b["nonphysical"]][:10],
        "junction_bonds": [b for b in bonds if b["junction"]],
        "junction_max_dev_A": round(max(jd), 3) if jd else None,
        "nonjunction_max_dev_A": round(max(nd), 3) if nd else None,
        "cn_mean_A": round(float(np.mean([b["C_N_A"] for b in bonds])), 4) if bonds else None,
        "cn_max_A": round(max((b["C_N_A"] for b in bonds), default=float("nan")), 4) if bonds else None,
        "cn_min_A": round(min((b["C_N_A"] for b in bonds), default=float("nan")), 4) if bonds else None,
    }


# --------------------------------------------------------------------------- 2. strain
def _bond_checks(resname: str, r: dict[str, Atom]) -> list[dict]:
    out = []
    for (tag, a, b), ideal in IDEAL_BOND.items():
        if tag == "peptide":
            continue
        if tag not in ("*", resname):
            continue
        if a not in r or b not in r:
            continue
        d = dist(r[a], r[b])
        out.append({"kind": "bond", "atoms": f"{a}-{b}", "value": round(d, 3),
                    "ideal": ideal, "dev": round(d - ideal, 3)})
    return out


def _angle_checks(resname: str, r: dict[str, Atom]) -> list[dict]:
    out = []
    for (tag, a, b, c), ideal in IDEAL_ANGLE.items():
        if tag not in ("*", resname):
            continue
        if a not in r or b not in r or c not in r:
            continue
        th = angle(r[a], r[b], r[c])
        out.append({"kind": "angle", "atoms": f"{a}-{b}-{c}",
                    "value": round(th, 2), "ideal": ideal,
                    "dev": round(th - ideal, 2)})
    return out


def backbone_dihedrals(atoms: list[Atom], key: tuple[str, int]) -> dict:
    """phi, psi, omega for one residue. None where the neighbour is absent."""
    res = residues(atoms)
    ch, seq = key
    prev, cur, nxt = res.get((ch, seq - 1)), res.get((ch, seq)), res.get((ch, seq + 1))
    out: dict[str, float | None] = {"phi": None, "psi": None, "omega": None}
    if cur is None:
        return out
    if prev and all(k in prev for k in ("C",)) and all(k in cur for k in ("N", "CA", "C")):
        out["phi"] = round(dihedral(prev["C"], cur["N"], cur["CA"], cur["C"]), 2)
    if nxt and "N" in nxt and all(k in cur for k in ("N", "CA", "C")):
        out["psi"] = round(dihedral(cur["N"], cur["CA"], cur["C"], nxt["N"]), 2)
    if nxt and all(k in nxt for k in ("N", "CA")) and all(k in cur for k in ("CA", "C")):
        out["omega"] = round(dihedral(cur["CA"], cur["C"], nxt["N"], nxt["CA"]), 2)
    return out


def ramachandran_region(resname: str, phi: float | None,
                        psi: float | None) -> dict:
    """Which broad Ramachandran region (phi, psi) falls in.

    Deliberately GENEROUS - the question is "is this conformation physically
    available to a real residue", not "is it the most probable rotamer". Only
    three regions are accepted for a general residue, covering essentially all
    of the populated map:

      beta / extended   phi < 0, psi > 90 or psi < -150
      alpha_R           phi < 0, -90 < psi < 30
      alpha_L           phi > 0 and roughly 30..100, psi roughly -10..90
                        (populated mainly by Gly, Asn and Asp)

    A POSITIVE phi outside alpha_L is the classic disallowed conformation for
    anything but glycine, which has no CB to clash.
    """
    out = {"phi": phi, "psi": psi, "region": None, "allowed": None,
           "note": ""}
    if phi is None or psi is None:
        out["note"] = "phi/psi not measurable (chain terminus or missing atom)"
        return out
    if resname == "GLY":
        out.update(region="gly_any", allowed=True,
                   note="glycine has no CB; the whole map is available")
        return out
    if phi < 0:
        if -90.0 <= psi <= 30.0:
            out.update(region="alpha_R", allowed=True)
        elif psi > 90.0 or psi < -150.0:
            out.update(region="beta", allowed=True)
        else:
            out.update(region="bridge", allowed=True,
                       note="sparsely populated but not disallowed")
        return out
    if 30.0 <= phi <= 100.0 and -10.0 <= psi <= 90.0:
        out.update(region="alpha_L", allowed=True,
                   note="left-handed; common for Asn/Asp, strained for others")
        return out
    # MEASURED, NOT ASSUMED. The catalytic residues of the reference Kemp
    # eliminases sit at POSITIVE phi with strongly negative psi:
    #   4BS0  Asp127 phi +87.9 psi  -99.8    Gln50 phi +70.6 psi -166.3
    #   7K4Z  Asp127 phi +82.7 psi -101.1    Gln50 phi +71.7 psi -167.1
    #   6TU6  Asp127 phi +87.1 psi  -95.3    Met50 phi +67.0 psi -164.6
    # Those are refined, deposited coordinates, so this region is unambiguously
    # physically available - an earlier version of this function called it
    # "disallowed" and would have condemned the reference enzyme itself.
    # Treated as unusual-but-real, and never gated on.
    if 40.0 <= phi <= 120.0:
        out.update(region="epsilon_unusual", allowed=True,
                   note="positive phi, negative psi: sparsely populated, but "
                        "exactly where the reference Kemp eliminases put their "
                        "catalytic Asp/Gln (4BS0, 7K4Z, 6TU6)")
        return out
    out.update(region="unusual", allowed=True,
               note=f"phi {phi:+.1f} is outside the common regions; reported, "
                    "not judged - the reference structures show this family's "
                    "catalytic residues occupy uncommon conformations")
    return out


def check_catalytic_strain(atoms: list[Atom], *, base_resseq: int,
                           oxy_resseq: int, chain: str = "A",
                           bond_tol: float = 0.05,
                           angle_tol: float = 8.0) -> dict:
    """Bond lengths, bond angles and dihedrals of the two catalytic residues.

    ``bond_tol`` / ``angle_tol`` are the per-term tolerances used to label a
    term strained. They are reported alongside every measurement so a caller
    can retune without re-measuring.
    """
    res = residues(atoms)
    out: dict[str, dict] = {}
    for role, seq in (("base", base_resseq), ("oxy", oxy_resseq)):
        key = (chain, seq)
        r = res.get(key)
        if not r:
            out[role] = {"resseq": seq, "present": False, "terms": [],
                         "n_strained": None}
            continue
        resname = _resname(atoms, key)
        terms = _bond_checks(resname, r) + _angle_checks(resname, r)
        for t in terms:
            tol = bond_tol if t["kind"] == "bond" else angle_tol
            t["strained"] = bool(abs(t["dev"]) > tol)
        dih = backbone_dihedrals(atoms, key)
        om = dih.get("omega")
        omega_ok = (om is None or min(abs(abs(om) - 180.0), abs(om)) <= OMEGA_TRANS_TOL)
        rama = ramachandran_region(resname, dih.get("phi"), dih.get("psi"))
        out[role] = {
            "resseq": seq,
            "resname": resname,
            "present": True,
            "n_atoms": len(r),
            "dihedrals": dih,
            "omega_planar": bool(omega_ok),
            "ramachandran": rama,
            "terms": terms,
            "n_strained": sum(1 for t in terms if t["strained"]),
            "worst_bond": max((t for t in terms if t["kind"] == "bond"),
                              key=lambda t: abs(t["dev"]), default=None),
            "worst_angle": max((t for t in terms if t["kind"] == "angle"),
                               key=lambda t: abs(t["dev"]), default=None),
        }
    strained = sum(v["n_strained"] or 0 for v in out.values())
    rama_bad = [r for r, v in out.items()
                if v.get("present") and (v.get("ramachandran") or {}).get("allowed") is False]
    return {"residues": out, "n_strained_total": strained,
            "ramachandran_outliers": rama_bad,
            "tolerances": {"bond_A": bond_tol, "angle_deg": angle_tol,
                           "omega_deg": OMEGA_TRANS_TOL}}


# --------------------------------------------------------------------------- 3. catalytic distances
def check_catalytic_distances(atoms: list[Atom], *, base_resseq: int,
                              oxy_resseq: int, chain: str = "A",
                              tol_A: float = 0.5,
                              tol_deg: float = 20.0) -> dict:
    """d_base, d_oxy, the CA-CA scaffolding distance and the OD2-N3-N1 angle.

    d_base and d_oxy need the ligand. When RFD3 was run without it they come
    back ``None`` with ``measurable: False`` - never as a pass.
    """
    res = residues(atoms)
    base = res.get((chain, base_resseq)) or {}
    oxy = res.get((chain, oxy_resseq)) or {}
    lig = {a.name: a for a in atoms if a.resname in LIGAND_RESNAMES}

    def term(name, value, ref, tol):
        ok = value is not None and abs(value - ref) <= tol
        return {"value": round(value, 3) if value is not None else None,
                "reference": ref, "tolerance": tol,
                "deviation": round(value - ref, 3) if value is not None else None,
                "measurable": value is not None, "within_tol": bool(ok)}

    d_base = dist(base["OD2"], lig["N3"]) if "OD2" in base and "N3" in lig else None
    d_oxy = dist(oxy["NE2"], lig["N1"]) if "NE2" in oxy and "N1" in lig else None
    ca_ca = (dist(base["CA"], oxy["CA"])
             if "CA" in base and "CA" in oxy else None)
    ang = (angle(base["OD2"], lig["N3"], lig["N1"])
           if "OD2" in base and "N3" in lig and "N1" in lig else None)

    ref = C.REF_GEOMETRY
    return {
        "ligand_present": bool(lig),
        "ligand_atoms": len(lig),
        "d_base": term("d_base", d_base, ref["asp127_od2__6nt_n3"], tol_A),
        "d_oxy": term("d_oxy", d_oxy, ref["gln50_ne2__6nt_n1"], tol_A),
        "ca_ca": term("ca_ca", ca_ca, ref["asp127_ca__gln50_ca"], 1.0),
        "od2_n3_n1_angle_deg": term("angle", ang, ref["od2_n3_n1_angle_deg"],
                                    tol_deg),
    }


# --------------------------------------------------------------------------- the gate
#: TWO TIERS, because the measured data demands two.
#:
#: FATAL rejects the design. A cut chain or a bond >0.30 A off ideal is not
#: something LigandMPNN repacks or Boltz relaxes - the chain has to be rebuilt,
#: and anything downstream is measuring a different molecule.
#:
#: WARN records and ranks. Junction distortion of 0.15-0.30 A is present in
#: EVERY design measured so far (see the calibration note above), so rejecting
#: on it would reject the whole sweep and produce nothing. It is reported as a
#: magnitude (``junction_max_dev_A``) so designs can be ordered by it and the
#: least-distorted ones preferred, rather than silently treated as equal.
#:
#: This split is deliberate: a gate that rejects everything is indistinguishable
#: from a broken pipeline, and a gate that warns about everything is noise. The
#: fatal tier is placed where the physics is unambiguous.
STRAIN_GATE = {
    "n_cuts_max": 0,              # measured C-N > 1.80 A: fatal, no tolerance
    "n_nonphysical_max": 0,       # |dev| > 0.30 A: fatal
    "n_strained_max": 4,          # strained bond/angle terms across both residues
    "require_omega_planar": True,
    "require_catalytic_measurable": True,   # no ligand in the output == cannot verify
    "require_catalytic_within_tol": True,
}


@dataclass
class Validation:
    passed: bool
    failed: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    integrity: dict = field(default_factory=dict)
    strain: dict = field(default_factory=dict)
    distances: dict = field(default_factory=dict)

    @property
    def junction_max_dev_A(self) -> float | None:
        return self.integrity.get("junction_max_dev_A")

    def summary(self) -> str:
        i, s, d = self.integrity, self.strain, self.distances
        db = (d.get("d_base") or {}).get("value")
        do = (d.get("d_oxy") or {}).get("value")
        verdict = "PASS" if self.passed else "FAIL:" + ",".join(self.failed)
        warn = (" WARN:" + ",".join(self.warnings)) if self.warnings else ""
        return (f"cuts={i.get('n_cuts')} nonphys={i.get('n_nonphysical')} "
                f"junction_dev={i.get('junction_max_dev_A')} "
                f"elsewhere={i.get('nonjunction_max_dev_A')} "
                f"strained={s.get('n_strained_total')} "
                f"d_base={db} d_oxy={do} {verdict}{warn}")

    def as_dict(self) -> dict:
        """Flat enough for the ledger, complete enough to re-derive the verdict."""
        return {
            "passed": self.passed,
            "failed": self.failed,
            "warnings": self.warnings,
            "junction_max_dev_A": self.junction_max_dev_A,
            "nonjunction_max_dev_A": self.integrity.get("nonjunction_max_dev_A"),
            "n_cuts": self.integrity.get("n_cuts"),
            "n_nonphysical": self.integrity.get("n_nonphysical"),
            "n_distorted": self.integrity.get("n_distorted"),
            "n_strained_total": self.strain.get("n_strained_total"),
            "d_base_A": (self.distances.get("d_base") or {}).get("value"),
            "d_oxy_A": (self.distances.get("d_oxy") or {}).get("value"),
            "ligand_present": self.distances.get("ligand_present"),
        }


def validate(structure: Path | str, *, base_resseq: int, oxy_resseq: int,
             chain: str = "A", thresholds: dict | None = None) -> Validation:
    """Run all three checks on an RFD3 output and return a gate verdict.

    Fails closed: an unmeasurable catalytic distance counts as a failure when
    ``require_catalytic_measurable`` is set, because "we could not check" and
    "it is fine" are not the same claim.
    """
    t = {**STRAIN_GATE, **(thresholds or {})}
    atoms = load_atoms(structure)
    motif = {base_resseq, oxy_resseq}
    integrity = check_chain_integrity(atoms, motif_resseq=motif)
    strain = check_catalytic_strain(atoms, base_resseq=base_resseq,
                                    oxy_resseq=oxy_resseq, chain=chain)
    distances = check_catalytic_distances(atoms, base_resseq=base_resseq,
                                          oxy_resseq=oxy_resseq, chain=chain)

    failed: list[str] = []
    warnings: list[str] = []

    # --- 1. chain integrity, measured not reported
    if (integrity["n_cuts"] or 0) > t["n_cuts_max"]:
        failed.append("chain_cut")
    if (integrity["n_nonphysical"] or 0) > t["n_nonphysical_max"]:
        failed.append("nonphysical_peptide_bond")
    jdev = integrity.get("junction_max_dev_A")
    if jdev is not None and jdev > CN_DEV_DISTORTED:
        warnings.append(f"motif_junction_distorted:{jdev:.3f}A")

    # --- 2. strain in the catalytic residues
    for role, v in strain["residues"].items():
        if not v.get("present"):
            failed.append(f"{role}_residue_missing")
        elif t["require_omega_planar"] and not v.get("omega_planar", True):
            failed.append(f"{role}_omega_nonplanar")
    if (strain["n_strained_total"] or 0) > t["n_strained_max"]:
        failed.append("catalytic_strain")
    # NO Ramachandran warning. The region boundaries were calibrated against
    # the reference crystal structures (4BS0/7K4Z/6TU6) and those put the real
    # catalytic Asp/Gln at positive phi, so "unusual phi" carries no signal for
    # this site. phi/psi stay in the record as measurements.

    # --- 3. catalytic distances
    for name in ("d_base", "d_oxy", "ca_ca", "od2_n3_n1_angle_deg"):
        term = distances[name]
        if not term["measurable"]:
            if t["require_catalytic_measurable"] and name in ("d_base", "d_oxy"):
                failed.append(f"{name}_unmeasurable")
            continue
        if t["require_catalytic_within_tol"] and not term["within_tol"]:
            failed.append(f"{name}_off")

    return Validation(passed=not failed, failed=failed, warnings=warnings,
                      integrity=integrity, strain=strain, distances=distances)
