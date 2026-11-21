#!/usr/bin/env python3
"""Kemp filter waterfall for one designed complex.

Hard fails only on organizer rules and on a broken active site:
length, alphabet, identity, more than one 6NT, missing Asp/Gln/ligand atoms,
or a non-catalytic heavy atom inside 1.5 Å of the ligand.

Geometry and burial are reported against 4BS0. They are not pass/fail:
the slides did not publish a numeric cutoff. iPDE, pLDDT, and 5-seed
consistency stay pending until a Boltz summary exists.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from pdb_io import dist, find, parse_atoms  # noqa: E402
from validate_design import DEFAULT_REFS, check_one, load_fasta  # noqa: E402

CRYSTAL_BASE = 2.53
CRYSTAL_DONOR = 3.10
THEOZYME_BASE = 2.482
THEOZYME_DONOR = 3.331
CLASH_FAIL_A = 1.5
HYDROPHOBIC = set("AILMFWYV")


def ligand_copies(atoms: list[dict], resname: str) -> list[tuple[str, int]]:
    seen: list[tuple[str, int]] = []
    for atom in find(atoms, rn=resname, het=True):
        key = (atom["ch"], atom["ri"])
        if key not in seen:
            seen.append(key)
    return seen


def closest(cands: list[dict], target: dict) -> dict | None:
    if not cands:
        return None
    return min(cands, key=lambda atom: dist(atom, target))


def sequence_from_protein(atoms: list[dict], chain: str) -> str:
    order = {"ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G",
             "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N",
             "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S", "THR": "T", "VAL": "V",
             "TRP": "W", "TYR": "Y"}
    residues: dict[int, str] = {}
    for atom in find(atoms, ch=chain, an="CA", het=False):
        residues[atom["ri"]] = order.get(atom["rn"], "X")
    return "".join(residues[i] for i in sorted(residues))


def sequons(seq: str) -> list[str]:
    hits = []
    for i in range(len(seq) - 2):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in "ST":
            hits.append(f"{i + 1}{seq[i:i + 3]}")
    return hits


def filter_complex(pdb: Path, chain: str | None, refs: list[tuple[str, str]]) -> dict:
    atoms = parse_atoms(pdb, chain=chain, alts=(" ", "A", "B"))
    errors: list[str] = []
    warnings: list[str] = []
    copies = ligand_copies(atoms, "6NT")
    if len(copies) != 1:
        errors.append(f"expected one 6NT, found {len(copies)}")
    ligand_chain, ligand_res = copies[0] if len(copies) == 1 else ("", 0)
    try:
        n3 = next(a for a in atoms if a["rn"] == "6NT" and a["an"] == "N3" and a["het"])
        n1 = next(a for a in atoms if a["rn"] == "6NT" and a["an"] == "N1" and a["het"])
    except StopIteration:
        errors.append("6NT is missing N1 or N3")
        n3 = n1 = None

    base = donor = None
    if n3 is not None:
        base = closest(find(atoms, rn="ASP", an="OD2", het=False), n3)
    if n1 is not None:
        donor = closest(find(atoms, rn="GLN", an="NE2", het=False), n1)
    if base is None:
        errors.append("no ASP OD2 near 6NT")
    if donor is None:
        errors.append("no GLN NE2 near 6NT")

    geometry = None
    if base is not None and donor is not None and n3 is not None and n1 is not None:
        base_d = dist(base, n3)
        donor_d = dist(donor, n1)
        geometry = {
            "base": f"{base['ch']}:ASP{base['ri']} OD2",
            "donor": f"{donor['ch']}:GLN{donor['ri']} NE2",
            "ASP_OD2_to_6NT_N3_A": round(base_d, 3),
            "GLN_NE2_to_6NT_N1_A": round(donor_d, 3),
            "delta_vs_4BS0_A": {
                "base": round(base_d - CRYSTAL_BASE, 3),
                "donor": round(donor_d - CRYSTAL_DONOR, 3),
            },
            "delta_vs_Theozyme_withORI_A": {
                "base": round(base_d - THEOZYME_BASE, 3),
                "donor": round(donor_d - THEOZYME_DONOR, 3),
            },
        }

    burial = None
    prot_chain = base["ch"] if base is not None else chain
    if base is not None and prot_chain:
        base_ca = next((a for a in find(atoms, ch=prot_chain, ri=base["ri"], an="CA", het=False)), None)
        if base_ca is None:
            errors.append("catalytic Asp has no CA")
        else:
            cas = find(atoms, ch=prot_chain, an="CA", het=False)
            n_incl = sum(1 for ca in cas if dist(ca, base_ca) <= 12.0)
            burial = {
                "n_CA_within_12A_excluding_self": n_incl - 1,
                "crystal_4BS0_excluding_self": 47,
            }
            if n_incl - 1 < 8:
                warnings.append("base burial is motif-scale; not a scaffold yet")

    clash = None
    if n3 is not None and base is not None and donor is not None:
        ligand_atoms = [a for a in atoms if a["het"] and a["rn"] == "6NT" and a["heavy"]]
        protected = {(base["ch"], base["ri"]), (donor["ch"], donor["ri"])}
        others = [
            a for a in atoms
            if (not a["het"]) and a["heavy"] and (a["ch"], a["ri"]) not in protected
        ]
        if others and ligand_atoms:
            minimum = min(dist(a, b) for a in others for b in ligand_atoms)
            clash = {"min_noncatalytic_heavy_to_6NT_A": round(minimum, 3)}
            if minimum < CLASH_FAIL_A:
                errors.append(f"non-catalytic clash with 6NT at {minimum:.2f} Å")

    seq = sequence_from_protein(atoms, prot_chain) if prot_chain else ""
    sequence = check_one(pdb.stem, seq, refs) if seq else {
        "length": 0, "max_identity": None, "closest_published": "",
        "errors": ["no protein sequence in the PDB"], "warnings": [], "pass": False,
    }
    errors.extend(sequence["errors"])
    warnings.extend(sequence["warnings"])
    if seq:
        hits = sequons(seq)
        if hits:
            warnings.append("N-X-S/T sequons: " + ",".join(hits))
        hydro = sum(aa in HYDROPHOBIC for aa in seq) / len(seq)
        if hydro > 0.45 and len(seq) >= 40:
            warnings.append(f"hydrophobic fraction {hydro:.0%}")

    pending = ["iPDE", "pLDDT", "pocket_PAE", "consistency_5_seeds", "TM_vs_published"]
    hard_pass = not errors
    return {
        "pdb": str(pdb),
        "hard_pass": hard_pass,
        "lab_ready": False,
        "lab_ready_reason": "iPDE, pLDDT, and 5-seed consistency are still pending",
        "errors": errors,
        "warnings": warnings,
        "geometry": geometry,
        "burial": burial,
        "clash": clash,
        "sequence": {
            "length": sequence["length"],
            "max_identity": sequence["max_identity"],
            "closest_published": sequence["closest_published"],
            "sequence": seq,
        },
        "pending": pending,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Score one Kemp complex with the local filter waterfall.")
    parser.add_argument("pdb", nargs="?", type=Path)
    parser.add_argument("--chain", default=None, help="protein/ligand chain; default keeps every chain")
    parser.add_argument("--refs", type=Path, default=DEFAULT_REFS)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    refs = load_fasta(args.refs) if args.refs.exists() else []

    if args.self_test:
        motif = HERE.parent / "theozyme" / "Theozyme_withORI.pdb"
        report = filter_complex(motif, chain="G", refs=refs)
        geometry = report["geometry"] or {}
        ok = (
            report["hard_pass"]
            and geometry.get("base", "").endswith("ASP127 OD2")
            and geometry.get("donor", "").endswith("GLN50 NE2")
            and abs(geometry.get("ASP_OD2_to_6NT_N3_A", 99) - THEOZYME_BASE) < 0.02
            and abs(geometry.get("GLN_NE2_to_6NT_N1_A", 99) - THEOZYME_DONOR) < 0.02
            and report["lab_ready"] is False
        )
        print(json.dumps(report, indent=2))
        print("self-test", "OK" if ok else "BROKEN")
        return 0 if ok else 1

    if args.pdb is None:
        parser.error("provide a PDB or --self-test")
    report = filter_complex(args.pdb, args.chain, refs)
    json.dump(report, sys.stdout, indent=2)
    print()
    return 0 if report["hard_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
