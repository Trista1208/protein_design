"""Can this pocket host the ligand at the IDEAL catalytic geometry?

This is the question Boltz cannot answer for us. Boltz folds at pLDDT 0.41 on
our designs and swings d_oxy by 4 A between runs, so "Boltz put the ligand
somewhere else" tells us nothing about whether the pocket COULD hold it.

So we stop asking the predictor and place the pose ourselves:

  1. Superpose the theozyme's FUNCTIONAL ATOMS - Asp carboxylate (CG, OD1, OD2)
     and Gln amide (CD, OE1, NE2) - onto the design's own catalytic residues.
     The resulting RMSD is itself a headline number: it measures whether the
     design presents the two functional groups in the theozyme's mutual
     arrangement at all. If it is large, NO rigid ligand placement can satisfy
     both contacts, whatever the predictor does.
  2. Apply that transform to the theozyme's 6NT. The ligand now sits at exactly
     the reference catalytic geometry (d_base 2.48, d_oxy 3.33) relative to
     this design's residues - by construction, which the script verifies.
  3. Ask whether that pose clashes with the rest of the protein.

!!! THE CLASH VERDICT FROM RIGID TRANSPLANT IS INVALID - DO NOT USE IT !!!

Measured against the controls: 4BS0 and 7K4Z are refined crystal structures
with 6NT actually bound, and this method reports 21 clashes for each, with
contacts down to 0.393 A. In the crystal its own ligand sits 2.99-3.48 A from
Met237 - ordinary vdW contact. The transplanted ligand matches the crystal's
own ligand to 0.32 A at N2/N3 but diverges by up to 3.93 A at the nitro group.

WHY: two H-bond contacts DO NOT DETERMINE A POSE. They pin the triazole
nitrogens, leaving free rotation about the axis through the contact atoms, plus
a possible ring flip. The theozyme's choice of that free rotation is arbitrary
and does not transfer to another pocket. So a rigid transplant places the
contact end correctly and the far end anywhere, and the resulting "clashes" are
an artefact of the arbitrary rotation.

CONSEQUENCE: single-point scoring of a hand-placed pose cannot work for this
problem. The remaining degrees of freedom have to be OPTIMISED subject to the
two distance restraints - restrained minimisation, not --score_only. See
scripts/restrained_score.py.

WHAT IS STILL VALID HERE: the functional-group superposition RMSD. It involves
only the protein's own Asp carboxylate and Gln amide, no ligand, so the
rotational ambiguity does not touch it. Validated against the controls:
  4BS0  0.348 A      7K4Z  0.341 A      (crystal-to-crystal scale)
That is the reference against which a design's fg-RMSD should be read.

AMBIGUOUS ATOM NAMING: the Asp carboxylate oxygens are chemically equivalent,
and the Gln amide OE1/NE2 assignment is a well-known modelling ambiguity (the
"amide flip"). Both assignments are tried and the better RMSD is reported,
because penalising a design for an arbitrary naming choice would be measuring
the file, not the structure.

NO force field is used here - this is pure geometry, so it needs nothing
installed and cannot be wrong about parameters. The energetic score
(the `--score_only` equivalent) is a separate step.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from kemp import config as C
from kemp import strain

ROOT = Path(__file__).resolve().parents[1]
THEOZYME = ROOT / "theozymes-discovered" / "Theozyme_withORI_bfix.pdb"

#: Atoms defining each contact geometry. Order matters: paired element-wise.
ASP_ATOMS = ["CG", "OD1", "OD2"]
GLN_ATOMS = ["CD", "OE1", "NE2"]
#: The H-bond pairs themselves - excluded from the clash check, since being
#: close is the entire point of them.
CONTACT_PAIRS = {("OD1", "N3"), ("OD2", "N3"), ("NE2", "N1"), ("OE1", "N1")}

#: Clash thresholds, Angstrom. Two N/O atoms may legitimately sit at H-bond
#: range; anything involving carbon or sulphur may not.
POLAR = set("NO")
MIN_POLAR_POLAR = 2.50
MIN_OTHER = 2.80
SEVERE = 2.20


def kabsch(mobile: np.ndarray, target: np.ndarray):
    """Rotation+translation taking `mobile` onto `target`, and the RMSD."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    p, q = mobile - mc, target - tc
    v, s, wt = np.linalg.svd(p.T @ q)
    d = np.sign(np.linalg.det(v @ wt))
    dm = np.diag([1.0, 1.0, d])
    r = (v @ dm @ wt).T
    aligned = (r @ p.T).T
    rmsd = float(np.sqrt(((aligned - q) ** 2).sum(axis=1).mean()))
    return r, mc, tc, rmsd


def element_of(name: str) -> str:
    for ch in name:
        if ch.isalpha():
            return ch.upper()
    return "X"


def load_theozyme():
    atoms = strain.load_atoms(THEOZYME)
    res = strain.residues(atoms)
    asp = gln = lig = None
    for key, group in res.items():
        rn = strain._resname(atoms, key)
        if rn == "ASP":
            asp = group
        elif rn == "GLN":
            gln = group
        elif rn in strain.LIGAND_RESNAMES:
            lig = group
    if not (asp and gln and lig):
        raise SystemExit(f"theozyme incomplete in {THEOZYME}")
    return asp, gln, lig


def feasibility(design: Path, base_resseq: int, oxy_resseq: int,
                chain: str = "A") -> dict:
    t_asp, t_gln, t_lig = load_theozyme()
    atoms = strain.load_atoms(design)
    res = strain.residues(atoms)
    d_asp = res.get((chain, base_resseq))
    d_gln = res.get((chain, oxy_resseq))
    out: dict = {"design": design.name, "base": f"{chain}{base_resseq}",
                 "oxy": f"{chain}{oxy_resseq}"}
    if not d_asp or not d_gln:
        out["error"] = "catalytic residues not found"
        return out

    # try both ambiguous assignments, keep the better fit
    best = None
    for asp_names in (ASP_ATOMS, ["CG", "OD2", "OD1"]):
        for gln_names in (GLN_ATOMS, ["CD", "NE2", "OE1"]):
            try:
                mob = np.array([t_asp[n].xyz for n in asp_names]
                               + [t_gln[n].xyz for n in gln_names])
                tgt = np.array([d_asp[n].xyz for n in asp_names]
                               + [d_gln[n].xyz for n in gln_names])
            except KeyError as exc:
                out["error"] = f"missing atom {exc}"
                return out
            r, mc, tc, rmsd = kabsch(mob, tgt)
            if best is None or rmsd < best[-1]:
                best = (asp_names, gln_names, r, mc, tc, rmsd)
    asp_names, gln_names, r, mc, tc, rmsd = best
    out["superposition_rmsd_A"] = round(rmsd, 3)
    out["assignment"] = {"asp": asp_names, "gln": gln_names}

    # transplant the ligand through the same transform
    placed = {}
    for name, a in t_lig.items():
        if element_of(name) == "H":
            continue
        placed[name] = (r @ (a.v - mc)) + tc

    # verify the construction: the contacts should now BE the reference values
    def dist_to(group, aname, lname):
        if aname not in group or lname not in placed:
            return None
        return float(np.linalg.norm(group[aname].v - placed[lname]))

    out["placed_d_base"] = (lambda v: None if v is None else round(v, 3))(
        dist_to(d_asp, "OD2", "N3"))
    out["placed_d_oxy"] = (lambda v: None if v is None else round(v, 3))(
        dist_to(d_gln, "NE2", "N1"))
    out["reference_d_base"] = C.REF_GEOMETRY["asp127_od2__6nt_n3"]
    out["reference_d_oxy"] = C.REF_GEOMETRY["gln50_ne2__6nt_n1"]

    # clash check against every protein heavy atom
    contacts = []
    for lname, lxyz in placed.items():
        le = element_of(lname)
        for a in atoms:
            if a.resname in strain.LIGAND_RESNAMES:
                continue              # the design's own ligand is replaced
            if element_of(a.name) == "H":
                continue
            if (a.chain, a.resseq) in ((chain, base_resseq), (chain, oxy_resseq)):
                if (a.name, lname) in CONTACT_PAIRS:
                    continue          # the H-bonds themselves
            d = float(np.linalg.norm(a.v - lxyz))
            if d > 4.0:
                continue
            ae = element_of(a.name)
            limit = (MIN_POLAR_POLAR if (ae in POLAR and le in POLAR)
                     else MIN_OTHER)
            if d < limit:
                contacts.append({
                    "ligand_atom": lname,
                    "protein": f"{strain._resname(atoms, (a.chain, a.resseq))}"
                               f"{a.resseq}:{a.name}",
                    "dist_A": round(d, 3), "limit_A": limit,
                    "severe": bool(d < SEVERE),
                })
    contacts.sort(key=lambda c: c["dist_A"])
    out["n_clashes"] = len(contacts)
    out["n_severe"] = sum(1 for c in contacts if c["severe"])
    out["worst"] = contacts[:6]
    out["min_contact_A"] = contacts[0]["dist_A"] if contacts else None
    # Reported for diagnosis only. NOT a feasibility verdict - see the module
    # docstring: the controls fail this check, so it cannot discriminate.
    out["clash_check_valid"] = False
    out["clash_check_note"] = ("artefact of the arbitrary rotation about the "
                               "contact axis; crystals fail this too")
    return out


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        print("usage: pose_feasibility.py <structure.cif> [<structure.cif> ...]")
        return 2

    results = []
    for path in argv[1:]:
        p = Path(path)
        atoms = strain.load_atoms(p)
        res = strain.residues(atoms)
        lig = {a.name: a for a in atoms if a.resname in strain.LIGAND_RESNAMES}
        # locate the catalytic pair the same way the scoring probe does
        base = oxy = None
        if "N3" in lig:
            cands = [(strain.dist(res[k][an], lig["N3"]), k) for k in res
                     for an in ("OD1", "OD2")
                     if strain._resname(atoms, k) == "ASP" and an in res[k]]
            if cands:
                base = min(cands)[1]
        if "N1" in lig:
            cands = [(strain.dist(res[k]["NE2"], lig["N1"]), k) for k in res
                     if strain._resname(atoms, k) == "GLN" and "NE2" in res[k]]
            if cands:
                oxy = min(cands)[1]
        if base is None or oxy is None:
            print(f"{p.name}: could not locate the catalytic pair")
            continue
        r = feasibility(p, base[1], oxy[1], chain=base[0])
        results.append(r)

        print(f"\n=== {r['design']}   base={r['base']} oxy={r['oxy']}")
        if r.get("error"):
            print(f"    ERROR {r['error']}")
            continue
        print(f"    functional-group superposition RMSD : {r['superposition_rmsd_A']} A")
        print(f"    transplanted pose d_base/d_oxy      : "
              f"{r['placed_d_base']} / {r['placed_d_oxy']} "
              f"(reference {r['reference_d_base']} / {r['reference_d_oxy']})")
        print(f"    clashes <limit: {r['n_clashes']}  (severe <{SEVERE} A: {r['n_severe']})"
              f"   closest contact {r['min_contact_A']} A")
        for c in r["worst"]:
            flag = "  SEVERE" if c["severe"] else ""
            print(f"       6NT {c['ligand_atom']:5} -- {c['protein']:18} "
                  f"{c['dist_A']:.3f} A  (limit {c['limit_A']}){flag}")
        print("    => clash counts above are DIAGNOSTIC ONLY (the crystal "
              "controls fail this check too)")

    print("\n\n================ SUMMARY ================")
    print(f"{'design':36} {'fg-RMSD':>9}   vs crystal reference 0.34-0.35 A")
    for r in results:
        if r.get("error"):
            continue
        rm = r["superposition_rmsd_A"]
        ratio = rm / 0.345
        print(f"{r['design']:36} {rm:>9}   {ratio:4.1f}x the crystal scale")
    print("""
HOW TO READ fg-RMSD: how well the design's Asp carboxylate and Gln amide
reproduce the theozyme's mutual arrangement. Small means a rigid ligand
placement can satisfy BOTH contacts at once; large means the two groups are
mispositioned relative to each other and no single ligand pose can.""")
    out = ROOT / "outputs" / "posefeas"
    out.mkdir(parents=True, exist_ok=True)
    (out / "results.json").write_text(json.dumps(results, indent=2, default=str))
    print(f"\nwrote {out / 'results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
