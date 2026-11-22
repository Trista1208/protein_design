"""Define the redesign shell for the d045 insurance track.

GOAL: bank >=5 compliant, foldable, submittable designs. The organisers take
five regardless of score (slide 7), and we currently have ZERO banked, with
exactly one design that folds above the pLDDT gate.

WHY A SHELL, AND WHY THIS ONE
-----------------------------
LESSONS 7.6 and the HER2/Catalyst method (LESSONS 4.1): do not let a model
write a whole sequence - perturb a parent that already folds. d045 folds at
pLDDT 0.963-0.967, so it is the only parent available.

denovoKemp (Methods) defines its design shell as 8 A around the ligand plus
the catalytic residues, with catalytic identities held fixed for the first
rounds. We use the same radius, measured here rather than assumed.

The shell is centred on BOTH the ligand and the Gln46 side chain, because
LOOP_DESIGN 6.2 showed `oxy` and `orient_oxy` are the same Gln46 tip: d_oxy is
NE2...N1 and orient_oxy is the dihedral CD-NE2...N1-N2. One chi2/chi3 change
moves both, worth 0.28 of the score. So the positions that can reposition that
rotamer are exactly the positions worth opening.

WHAT IS FROZEN
--------------
  * Asp102 and Gln46 - the catalytic pair. Never designable.
  * Everything outside the shell - so edits stay local and the fold survives.

This script only REPORTS the shell; it writes a JSON the design step consumes.
Nothing is mutated here.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C                                   # noqa: E402
from kemp.theozyme import load                                 # noqa: E402

PARENT_PDB = Path("outputs/refold_d045/seed2/folded.pdb")      # best seed
BASE_IDX, OXY_IDX = 102, 46
CHAIN = "A"
SHELL_A = 8.0                                                  # denovoKemp
OUT = Path("outputs/insurance")

#: Never designable. The catalytic pair IS the enzyme; LigandMPNN will happily
#: mutate it away without this (runners.ligandmpnn docstring, handover 4.3.1).
CATALYTIC = {BASE_IDX, OXY_IDX}


def main() -> int:
    if not PARENT_PDB.exists():
        print(f"parent structure missing: {PARENT_PDB}")
        return 2
    tz = load(str(PARENT_PDB), require_ligand=False)

    prot = [a for a in tz.atoms if a.chain == CHAIN]
    lig = [a for a in tz.atoms if a.resname == C.LIGAND_CCD]
    if not lig:
        print("no 6NT in the parent structure")
        return 2

    # centres: every ligand heavy atom, plus the Gln46 side-chain tip atoms
    centres = [a.xyz for a in lig]
    gln_tip = [a.xyz for a in prot
               if a.resseq == OXY_IDX and a.name in ("CD", "OE1", "NE2")]
    centres += gln_tip
    centres = np.array(centres)

    # residue -> min distance to any centre
    per_res: dict[int, tuple[float, str]] = {}
    for a in prot:
        d = float(np.min(np.linalg.norm(centres - a.xyz, axis=1)))
        cur = per_res.get(a.resseq)
        if cur is None or d < cur[0]:
            per_res[a.resseq] = (d, a.resname)

    shell = sorted(r for r, (d, _) in per_res.items()
                   if d <= SHELL_A and r not in CATALYTIC)
    designable = shell
    all_res = sorted(per_res)
    fixed = [r for r in all_res if r not in designable]

    print(f"parent: {PARENT_PDB}  ({len(all_res)} residues, chain {CHAIN})")
    print(f"centres: {len(lig)} ligand atoms + {len(gln_tip)} Gln{OXY_IDX} tip atoms")
    print(f"shell: <= {SHELL_A} A, catalytic Asp{BASE_IDX}/Gln{OXY_IDX} excluded\n")
    print(f"{len(designable)} designable positions:")
    for r in designable:
        d, rn = per_res[r]
        print(f"   {rn}{r:<4} {d:5.2f} A")
    print(f"\n{len(fixed)} positions frozen "
          f"({100*len(fixed)/len(all_res):.0f}% of the chain)")

    # LigandMPNN's fixed_residues form is space-separated "A10 A11"
    fixed_spec = " ".join(f"{CHAIN}{r}" for r in fixed)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "shell.json").write_text(json.dumps({
        "parent_pdb": str(PARENT_PDB),
        "chain": CHAIN,
        "shell_radius_A": SHELL_A,
        "shell_rationale": ("8 A around ligand + Gln46 side-chain tip; "
                            "denovoKemp Methods radius, measured not assumed"),
        "catalytic_frozen": sorted(CATALYTIC),
        "designable": designable,
        "designable_detail": {str(r): {"resname": per_res[r][1],
                                       "min_dist_A": round(per_res[r][0], 2)}
                              for r in designable},
        "n_fixed": len(fixed),
        "fixed_residues_spec": fixed_spec,
    }, indent=2))
    print(f"\nwrote {OUT / 'shell.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
