"""Refold Rosetta enzdes designs with Boltz-2. The only test that counts.

    python scripts/refold_rosetta.py <rosetta_run_dir> [seeds]

WHY THIS STEP IS NOT OPTIONAL
-----------------------------
Rosetta reports our constraints essentially satisfied on its own models -
allcst 0.01, atom_pair 0.00, angle 0.00, dihedral 0.01 - which no Boltz-folded
variant has achieved. But Rosetta achieved it by MOVING COORDINATES: repacking
the Gln rotamer and settling the ligand into the constrained pose.

LESSONS 6 already tested exactly that reasoning and it failed. A coordinate
scan found 360 clash-free Gln rotamers satisfying both the distance window and
the dihedral band, and 0 of 5 Boltz seeds sampled any of them. A geometry
Rosetta likes is a hypothesis about a SEQUENCE; the competition metric is
computed on what Boltz predicts from that sequence.

So the question here is narrow and empirical: does the designed SEQUENCE make
Boltz reproduce the geometry Rosetta built?

The bound-ligand filter applies as it does everywhere else. A dihedral
measured on a design whose ligand has drifted is noise, not an angle
(LESSONS 2.3) - that filter has already overturned two apparent breakthroughs
tonight, so it runs before any geometry is quoted.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C, objective, official, pipeline     # noqa: E402
from kemp.runners import boltz_predict                          # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BASE_IDX, OXY_IDX = 102, 46
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
         "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
         "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
         "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
IN_WINDOW = (2.0, 3.6)


def chain_sequence(pdb: Path, chain: str = "A") -> str:
    s: dict[int, str] = {}
    for ln in pdb.read_text().splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA" and ln[21] == chain:
            s[int(ln[22:26])] = THREE.get(ln[17:20].strip(), "X")
    return "".join(s[k] for k in sorted(s))


def main() -> int:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    run_dir = Path(sys.argv[1])
    n_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    pdbs = sorted(p for p in run_dir.glob("input_0*.pdb"))
    if not pdbs:
        raise SystemExit(f"no Rosetta output PDBs in {run_dir}")

    out = run_dir / "refold"
    out.mkdir(parents=True, exist_ok=True)
    site = official.CatalyticSite(base_resseq=BASE_IDX, oxy_resseq=OXY_IDX,
                                  chain="A")
    print(f"{len(pdbs)} Rosetta design(s), {n_seeds} seed(s) each\n")
    print("  Rosetta satisfied the constraints on its own models "
          "(allcst ~0.01).\n  The question is whether Boltz reproduces that "
          "from the sequence alone.\n")

    rows = []
    for i, pdb in enumerate(pdbs):
        seq = chain_sequence(pdb)
        tag = pdb.stem
        if len(seq) < 50:
            print(f"  {tag}: sequence too short ({len(seq)}) - skipped")
            continue
        if seq[BASE_IDX - 1] != "D" or seq[OXY_IDX - 1] != "Q":
            print(f"  {tag}: CATALYTIC PAIR LOST "
                  f"({seq[BASE_IDX-1]}{BASE_IDX}/{seq[OXY_IDX-1]}{OXY_IDX}) - skipped")
            continue
        for seed in C.BOLTZ_SEEDS[:n_seeds]:
            w = out / f"{tag}_seed{seed}"
            w.mkdir(parents=True, exist_ok=True)
            t0 = time.time()
            try:
                r = boltz_predict(seq, w, seed=seed, ligand_ccd=C.LIGAND_CCD)
                folded = next((p for p in r.structures), None)
                if folded is None:
                    print(f"  {tag} seed {seed}: FOLD FAILED")
                    continue
                p = (folded if folded.suffix == ".pdb"
                     else pipeline.cif_to_pdb(folded, w / "folded.pdb"))
                raw = official.measure_raw(str(p), site=site)
                conf = r.confidence or {}
                gap = objective.dihedral_band_gap(raw.dihedral12_oxy)
                bound = (raw.d_base is not None and raw.d_oxy is not None
                         and IN_WINDOW[0] <= raw.d_base <= IN_WINDOW[1]
                         and IN_WINDOW[0] <= raw.d_oxy <= IN_WINDOW[1])
                sc = official.score(raw, strict=False)
                rows.append({"design": tag, "seed": seed, "sequence": seq,
                             "d_base": raw.d_base, "d_oxy": raw.d_oxy,
                             "dih_oxy": raw.dihedral12_oxy, "gap": gap,
                             "bound": bound,
                             "plddt": conf.get("complex_plddt"),
                             "ipde": conf.get("complex_ipde"),
                             "subscores": sc["subscores"]})
                flag = ""
                if not bound:
                    flag = "   LIGAND NOT BOUND - geometry is noise (LESSONS 2.3)"
                elif gap == 0.0:
                    flag = "   <<< BOUND *AND* DIHEDRAL IN BAND"
                print(f"  {tag} seed {seed}: d_base {raw.d_base:6.3f}  "
                      f"d_oxy {raw.d_oxy:6.3f}  dih {raw.dihedral12_oxy:6.1f}  "
                      f"gap {gap:5.1f}  pLDDT {conf.get('complex_plddt', 0):.3f}"
                      f"  ({time.time()-t0:.0f}s){flag}", flush=True)
            except Exception as exc:                            # noqa: BLE001
                print(f"  {tag} seed {seed}: ERROR {type(exc).__name__}: "
                      f"{str(exc)[:120]}")
        (out / "refold.json").write_text(json.dumps(rows, indent=2, default=str))

    print("\n" + "=" * 72)
    bound = [r for r in rows if r["bound"]]
    print(f"folded {len(rows)};  ligand bound: {len(bound)}")
    if not bound:
        print("  Rosetta's constrained geometry did NOT survive refolding.")
        print("  That is LESSONS 6 again: a pose Rosetta can build is not a")
        print("  pose Boltz predicts from the sequence. Rosetta's value here")
        print("  is the POSITIONS it chose, not the coordinates it produced.")
    else:
        inband = [r for r in bound if r["gap"] == 0.0]
        print(f"  of those, dihedral in band: {len(inband)}")
        for r in sorted(bound, key=lambda x: x["gap"]):
            print(f"   {r['design']} seed {r['seed']}: gap {r['gap']:.1f} deg, "
                  f"pLDDT {r['plddt']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
