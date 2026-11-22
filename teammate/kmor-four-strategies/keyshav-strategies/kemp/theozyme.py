"""M0 - theozyme handling.

This module exists because four of the handover's traps are silent: they do not
raise, they just quietly waste a whole batch of designs.

  1. Mixed coordinate frames (handover 1.2). Using Theozyme_noORI.pdb as the RFD3
     input while pasting ORI from the withORI file puts the origin 70.5 A from the
     ligand. RFD3 does not error.
  2. The ORI "atom" is a visualisation, not an input (1.3). Its element field is
     3 characters, which is invalid PDB and breaks RDKit and OpenBabel outright.
  3. Altlocs are mixed - Asp 'A', ligand 'B', Gln blank (4.3.4). A parser that
     keeps only altloc 'A' or blank silently DROPS THE LIGAND.
  4. Gln50 is modelled, not observed (1.1), so its 3.33 A contact is an opinion.

Deliberately numpy-only so it runs in any of the project envs without pulling
biotite/RDKit into the import path.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

from . import config as C


class TheozymeError(RuntimeError):
    """Raised instead of degrading silently. Handover 4.3.9: fail closed."""


# --------------------------------------------------------------------------- model
@dataclass
class Atom:
    record: str
    serial: int
    name: str
    altloc: str
    resname: str
    chain: str
    resseq: int
    xyz: np.ndarray
    occupancy: float
    bfactor: float
    element: str
    raw: str

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.chain, self.resseq, self.name)


#: REMARK tag used to carry frame provenance into derived files. Without this a
#: cleaned/derived PDB loses the one fact that makes ``ori_for`` safe, and the
#: guard then (correctly but uselessly) refuses to run.
FRAME_REMARK = "REMARK 100 KEMP_FRAME"
FRAME_WITH_ORI = "withORI"
FRAME_NO_ORI = "noORI"


@dataclass
class Theozyme:
    path: Path
    atoms: list[Atom] = field(default_factory=list)
    ori: np.ndarray | None = None
    cryst1: str | None = None
    #: which supplied file's coordinate system these atoms live in
    frame: str | None = None

    # -- selectors ---------------------------------------------------------
    def by_resname(self, resname: str) -> list[Atom]:
        return [a for a in self.atoms if a.resname == resname]

    def atom(self, resname: str, name: str) -> Atom:
        hits = [a for a in self.atoms if a.resname == resname and a.name == name]
        if not hits:
            raise TheozymeError(f"atom {resname}.{name} not found in {self.path.name}")
        if len(hits) > 1:
            raise TheozymeError(
                f"atom {resname}.{name} is ambiguous ({len(hits)} copies) - "
                "normalise altlocs first"
            )
        return hits[0]

    def atom_at(self, resseq: int, name: str, *, chain: str | None = None) -> Atom:
        """Atom by residue NUMBER, not residue type.

        Required for real designs: a 240-residue protein has many ASPs, so
        ``atom("ASP", "OD2")`` is ambiguous and raises. The catalytic positions
        come from RFD3's ``diffused_index_map``.
        """
        hits = [a for a in self.atoms if a.resseq == resseq and a.name == name
                and (chain is None or a.chain == chain)]
        if not hits:
            raise TheozymeError(
                f"no atom {name} in residue {resseq}"
                + (f" chain {chain}" if chain else "") + f" of {self.path.name}")
        if len(hits) > 1:
            raise TheozymeError(f"atom {name} in residue {resseq} is ambiguous "
                                f"({len(hits)} copies) - normalise altlocs first")
        return hits[0]

    def resname_at(self, resseq: int, *, chain: str | None = None) -> str | None:
        for a in self.atoms:
            if a.resseq == resseq and (chain is None or a.chain == chain):
                return a.resname
        return None

    @property
    def ligand(self) -> list[Atom]:
        return self.by_resname(C.LIGAND_CCD)

    @property
    def ligand_centroid(self) -> np.ndarray:
        lig = self.ligand
        if not lig:
            raise TheozymeError(
                f"no {C.LIGAND_CCD} atoms survived parsing - this is trap #3 "
                "(altloc filtering dropped the ligand)"
            )
        return np.mean([a.xyz for a in lig], axis=0)

    @property
    def altlocs(self) -> set[str]:
        return {a.altloc for a in self.atoms}


# --------------------------------------------------------------------------- parsing
def _parse_line(line: str) -> Atom | None:
    if not line.startswith(("ATOM", "HETATM")):
        return None
    try:
        return Atom(
            record=line[0:6].strip(),
            serial=int(line[6:11]),
            name=line[12:16].strip(),
            altloc=line[16:17].strip(),
            resname=line[17:20].strip(),
            chain=line[21:22].strip(),
            resseq=int(line[22:26]),
            xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
            occupancy=float(line[54:60] or 0.0),
            bfactor=float(line[60:66] or 0.0),
            element=line[76:78].strip(),
            raw=line.rstrip("\n"),
        )
    except (ValueError, IndexError) as exc:  # malformed -> fail closed
        raise TheozymeError(f"unparseable PDB line: {line!r}") from exc


def load(path: Path | str, *, normalise_altlocs: bool = True,
         strip_ori: bool = True, require_ligand: bool = True) -> Theozyme:
    """Load a theozyme PDB safely.

    normalise_altlocs collapses to one conformer per (chain, resseq, name) while
    PRESERVING the ligand - the handover's trap #3. strip_ori removes the ORI
    pseudo-atom and stashes its coordinate on ``.ori``.
    """
    path = Path(path)
    if not path.exists():
        raise TheozymeError(f"theozyme not found: {path}")

    tz = Theozyme(path=path)
    # Frame provenance: filename for the supplied files, REMARK for derived ones.
    if path.name == C.THEOZYME_WITH_ORI.name:
        tz.frame = FRAME_WITH_ORI
    elif path.name == C.THEOZYME_NO_ORI.name:
        tz.frame = FRAME_NO_ORI

    raw_atoms: list[Atom] = []
    for line in path.read_text().splitlines():
        if line.startswith(FRAME_REMARK):
            tz.frame = line.split()[-1].strip() or tz.frame
            continue
        if line.startswith("CRYST1"):
            tz.cryst1 = line.rstrip()
            continue
        atom = _parse_line(line)
        if atom is None:
            continue
        if atom.resname == "ORI" or atom.name == "ORI":
            # trap #2: 3-char element field, invalid PDB. Capture then discard.
            tz.ori = atom.xyz.copy()
            if strip_ori:
                continue
        raw_atoms.append(atom)

    tz.atoms = _normalise_altlocs(raw_atoms) if normalise_altlocs else raw_atoms

    # trap #3 assertion: for a theozyme the ligand MUST have survived.
    # require_ligand=False for RFD3 backbones and relaxed intermediates, which
    # legitimately have no ligand yet.
    if require_ligand and not tz.by_resname(C.LIGAND_CCD):
        raise TheozymeError(
            f"{C.LIGAND_CCD} vanished during parsing of {path.name}. "
            "Altloc normalisation dropped the ligand (handover trap #3). "
            "If this structure genuinely has no ligand, pass require_ligand=False."
        )
    return tz


def _normalise_altlocs(atoms: Iterable[Atom]) -> list[Atom]:
    """Keep the highest-occupancy conformer per atom, whatever its altloc letter.

    Crucially this does NOT prefer altloc 'A'/blank, because the ligand here
    carries altloc 'B' and would be deleted.
    """
    best: dict[tuple[str, int, str], Atom] = {}
    for a in atoms:
        prev = best.get(a.key)
        if prev is None or a.occupancy > prev.occupancy:
            best[a.key] = a
    return sorted(best.values(), key=lambda a: a.serial)


# --------------------------------------------------------------------------- frames
def ori_for(source: "Theozyme | Path | str") -> np.ndarray:
    """Return the ORI coordinate expressed in the frame of ``source``.

    Handover 1.2: the two supplied files are in different frames (identical motif,
    RMSD 0.0005 A after superposition). There is exactly one correct answer per file.

    Accepts a Theozyme (preferred - carries ``.frame`` provenance even for derived
    files) or a path. Derived files keep provenance via the FRAME_REMARK line that
    ``write_clean_pdb`` emits, so cleaning a theozyme does not lose the one fact
    that makes this check possible.
    """
    frame = None
    if isinstance(source, Theozyme):
        frame, name = source.frame, source.path.name
    else:
        name = Path(source).name
        if name == C.THEOZYME_WITH_ORI.name:
            frame = FRAME_WITH_ORI
        elif name == C.THEOZYME_NO_ORI.name:
            frame = FRAME_NO_ORI
        elif Path(source).exists():
            for line in Path(source).read_text().splitlines():
                if line.startswith(FRAME_REMARK):
                    frame = line.split()[-1].strip()
                    break
                if line.startswith(("ATOM", "HETATM")):
                    break

    if frame == FRAME_WITH_ORI:
        return np.array(C.ORI_IN_WITHORI_FRAME)
    if frame == FRAME_NO_ORI:
        return np.array(C.ORI_IN_NOORI_FRAME)
    raise TheozymeError(
        f"unknown theozyme frame for {name!r}; refusing to guess an ORI coordinate. "
        "Mixing frames puts the origin ~70 A from the ligand and silently wastes "
        f"the entire batch. Derived files must carry a '{FRAME_REMARK} <frame>' line - "
        "write them with write_clean_pdb()."
    )


def assert_ori_sane(tz: Theozyme, ori: np.ndarray) -> float:
    """Handover 4.3.2: assert ||ori - ligand_centroid|| < 12 A before EVERY run."""
    d = float(np.linalg.norm(np.asarray(ori, dtype=float) - tz.ligand_centroid))
    if d >= C.ORI_LIGAND_MAX_DIST:
        raise TheozymeError(
            f"ORI is {d:.2f} A from the ligand centroid (limit "
            f"{C.ORI_LIGAND_MAX_DIST} A). This is the mixed-frame bug: RFD3 will "
            "not error, it will build proteins centred away from the active site."
        )
    return d


def kabsch_rmsd(p: np.ndarray, q: np.ndarray) -> float:
    """RMSD after optimal superposition. Used to prove two files share a motif."""
    p = p - p.mean(axis=0)
    q = q - q.mean(axis=0)
    v, _, wt = np.linalg.svd(p.T @ q)
    if (np.linalg.det(v) * np.linalg.det(wt)) < 0.0:
        v[:, -1] = -v[:, -1]
    rot = v @ wt
    diff = (p @ rot) - q
    return float(np.sqrt((diff * diff).sum() / len(p)))


def compare_frames(a: Theozyme, b: Theozyme) -> dict:
    """Confirm handover 1.2: same motif, different frame."""
    shared = sorted(set(x.key for x in a.atoms) & set(x.key for x in b.atoms))
    if not shared:
        raise TheozymeError("no shared atoms between the two theozyme files")
    pa = np.array([next(x.xyz for x in a.atoms if x.key == k) for k in shared])
    pb = np.array([next(x.xyz for x in b.atoms if x.key == k) for k in shared])
    return {
        "n_shared_atoms": len(shared),
        "superposed_rmsd_A": kabsch_rmsd(pa, pb),
        "raw_rmsd_A": float(np.sqrt(((pa - pb) ** 2).sum(axis=1).mean())),
    }


# --------------------------------------------------------------------------- geometry
def _angle_deg(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    v1, v2 = a - b, c - b
    cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def measure_geometry(tz: Theozyme) -> dict[str, float]:
    """Reproduce the handover's reference table from a theozyme file.

    Returns the same keys as config.REF_GEOMETRY so the two can be diffed.
    """
    asp_od2 = tz.atom("ASP", "OD2").xyz
    asp_od1 = tz.atom("ASP", "OD1").xyz
    asp_ca = tz.atom("ASP", "CA").xyz
    gln_ne2 = tz.atom("GLN", "NE2").xyz
    gln_ca = tz.atom("GLN", "CA").xyz
    n3 = tz.atom(C.LIGAND_CCD, "N3").xyz
    n1 = tz.atom(C.LIGAND_CCD, "N1").xyz

    out = {
        "asp127_od2__6nt_n3": float(np.linalg.norm(asp_od2 - n3)),
        "gln50_ne2__6nt_n1": float(np.linalg.norm(gln_ne2 - n1)),
        "asp127_ca__gln50_ca": float(np.linalg.norm(asp_ca - gln_ca)),
        "od2_n3_n1_angle_deg": _angle_deg(asp_od2, n3, n1),
        "od1__n3_nonengaging": float(np.linalg.norm(asp_od1 - n3)),
    }
    try:  # hydrogens are present in the supplied files but not guaranteed
        out["asp127_od2__h3"] = float(np.linalg.norm(asp_od2 - tz.atom(C.LIGAND_CCD, "H3").xyz))
    except TheozymeError:
        pass
    return out


def diff_reference(measured: dict[str, float], tol: float = 0.05) -> dict[str, dict]:
    """Compare measured geometry against the verified reference table."""
    report = {}
    for key, ref in C.REF_GEOMETRY.items():
        if key not in measured:
            continue
        dev = measured[key] - ref
        report[key] = {
            "measured": round(measured[key], 3),
            "reference": ref,
            "deviation": round(dev, 3),
            "within_tol": abs(dev) <= (tol if "angle" not in key else tol * 20),
            "modelled_not_observed": key in C.REF_GEOMETRY_MODELLED,
        }
    return report


def write_clean_pdb(tz: Theozyme, out: Path | str, *, sort_residues: bool = True,
                    ligand_chain: str | None = "L") -> Path:
    """Write a single-conformer, ORI-free PDB safe for RDKit/Boltz/RFD3/RF3.

    Handover 4.3.3: no ORI atom downstream, ever.

    Two things beyond stripping ORI:

    * ``FRAME_REMARK`` records which coordinate system these atoms live in, so
      ``ori_for`` still works on this derived file. Without it the frame guard
      fires on our own output.
    * ``sort_residues`` renumbers into non-decreasing residue order. The supplied
      theozyme is ASP 127, GLN 50, 6NT 401 - genuinely out of order - and RF3
      rejects that outright with
      ``ValueError: res_id values are not in non-decreasing order for chain 'G'``.
      Original numbering is preserved in a REMARK so the mapping is recoverable.
    """
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)

    atoms = [a for a in tz.atoms if a.resname != "ORI" and a.name != "ORI"]
    # Put the ligand on its OWN chain. Sharing a chain with the protein makes
    # RFD3 refuse to append it ("Ligand chain(s) overlap with existing
    # chain(s)"), and dropping the `ligand` field then lets RFD3 build backbone
    # straight THROUGH the ligand volume - measured at 0.46 A closest approach,
    # which is why every folded design docked 6NT 3.5-15 A away.
    if ligand_chain:
        for a in atoms:
            if a.resname == C.LIGAND_CCD:
                a.chain = ligand_chain
                a.raw = a.raw[:21] + ligand_chain + a.raw[22:]
    lines: list[str] = []
    if tz.frame:
        lines.append(f"{FRAME_REMARK} {tz.frame}")
    if tz.cryst1:
        lines.append(tz.cryst1)

    if sort_residues:
        # stable sort by (chain, original resseq); then renumber 1..N per chain
        atoms = sorted(atoms, key=lambda a: (a.chain, a.resseq, a.serial))
        mapping: dict[tuple[str, int], int] = {}
        counter: dict[str, int] = {}
        for a in atoms:
            key = (a.chain, a.resseq)
            if key not in mapping:
                counter[a.chain] = counter.get(a.chain, 0) + 1
                mapping[key] = counter[a.chain]
        for (chain, old), new in sorted(mapping.items()):
            lines.append(f"REMARK 100 KEMP_RESMAP {chain} {old} -> {new}")
        for i, a in enumerate(atoms, start=1):
            new_res = mapping[(a.chain, a.resseq)]
            raw = a.raw
            # serial (7-11), altloc (17), resseq (23-26)
            raw = f"{raw[:6]}{i:>5}{raw[11:16]} {raw[17:22]}{new_res:>4}{raw[26:]}"
            lines.append(raw)
    else:
        for a in atoms:
            lines.append(a.raw[:16] + " " + a.raw[17:])

    lines.append("END")
    out.write_text("\n".join(lines) + "\n")
    return out
