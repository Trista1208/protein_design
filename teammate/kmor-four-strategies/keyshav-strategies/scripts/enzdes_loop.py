"""The local enzdes optimisation loop: design -> refold -> filter -> score.

    python scripts/enzdes_loop.py [n_parents] [nstruct] [--untargeted]

WHY LOCALLY, INSTEAD OF FuncLib
-------------------------------
FuncLib's defining step is a PSSM from a multiple sequence alignment: it
restricts mutations to ones evolution tolerates at each position, then applies
Rosetta ddG, then enumerates combinations. Our designs are de novo and the MSA
test (LOOP_DESIGN 6.1) measured their alignment depth at exactly 1 - the
ColabFold server returns the query and nothing else. FuncLib's evolutionary
filter therefore has nothing to work with here.

Rosetta enzdes design has no such dependency. It scores catalytic geometry and
energy directly, which is what we actually need.

THE LOOP
--------
  1. PARENTS  - designs from the screen with a BOUND ligand (both distances
                inside 2.0-3.6 A) and pLDDT >= C.PLDDT_MIN. Anything else is
                unmeasurable, and its geometry is noise (LESSONS 2.3, 7.1).
  2. DESIGN   - Rosetta FastDesign under the enzdes .cst, which scores
                distance AND torsion in one energy. This is the thing neither
                LigandMPNN (sequence likelihood) nor a steric ladder could do,
                and the reason every single mutant so far fixed one term while
                breaking the other.
  3. REFOLD   - Boltz-2, 1 seed. NOT optional. Rosetta satisfies the
                constraints by moving coordinates; LESSONS 6 showed 360
                clash-free in-band rotamers existed that 0/5 Boltz seeds ever
                sampled. The metric is computed on a Boltz prediction, so a
                Rosetta pose is a hypothesis about a sequence.
  4. FILTER   - bound ligand, then pLDDT, then geometry. In that order.
  5. RANK     - dense percentile objective, so a batch where every official
                score is 0.000 still orders (kemp/objective.py).

Survivors go to scripts/score_finalist.py for the 5-seed official score.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C, objective, official, pipeline     # noqa: E402
from kemp.runners import boltz_predict                          # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RO = ROOT / "outputs/rosetta"
SCREEN = ROOT / "outputs/insurance/screen.json"
OUT = ROOT / "outputs/enzdes_loop"
BASE_IDX, OXY_IDX = 102, 46
IN_WINDOW = (2.0, 3.6)

#: independently supported by the bound-ligand screen AND by the untargeted
#: Rosetta run, which chose F26L/M, F42M, V84L/W, F103L unprompted
TARGET_POSITIONS = [23, 26, 27, 42, 45, 49, 84, 99, 103]

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
         "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
         "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
         "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def bound(d_base, d_oxy) -> bool:
    return (d_base is not None and d_oxy is not None
            and IN_WINDOW[0] <= d_base <= IN_WINDOW[1]
            and IN_WINDOW[0] <= d_oxy <= IN_WINDOW[1])


def pick_parents(n: int) -> list[dict]:
    """Bound ligand + measurable fold, best distances first."""
    rows = [r for r in json.loads(SCREEN.read_text())
            if r.get("status") == "ok"]
    ok = [r for r in rows
          if bound(r.get("d_base"), r.get("d_oxy"))
          and (r.get("plddt") or 0) >= C.PLDDT_MIN]
    tb, to = official.TARGETS["d_base"], official.TARGETS["d_oxy"]
    ok.sort(key=lambda r: -(official._window_score(r["d_base"], tb["target"], *tb["window"])
                            * official._window_score(r["d_oxy"], to["target"], *to["window"])))
    return ok[:n]


def chain_sequence(pdb: Path, chain: str = "A") -> str:
    s: dict[int, str] = {}
    for ln in pdb.read_text().splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA" and ln[21] == chain:
            s[int(ln[22:26])] = THREE.get(ln[17:20].strip(), "X")
    return "".join(s[k] for k in sorted(s))


def run_design(parent_pdb: Path, label: str, nstruct: int,
               targeted: bool) -> list[Path]:
    """One Rosetta enzdes FastDesign job. Returns the design PDBs."""
    cmd = [sys.executable, str(ROOT / "scripts/run_rosetta_design.py"),
           str(parent_pdb), str(nstruct), label]
    env_xml = "enzdes_targeted.xml" if targeted else "enzdes_design.xml"
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=10800,
                       env={**__import__("os").environ,
                            "KEMP_ROOT": str(ROOT),
                            "ENZDES_XML": env_xml,
                            "ENZDES_DES_POS": ",".join(
                                str(x) for x in TARGET_POSITIONS)})
    if p.returncode != 0:
        print(f"    rosetta failed: {p.stdout[-300:]}")
        return []
    return sorted((RO / "runs" / label).glob("input_0*.pdb"))


def main() -> int:
    n_parents = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    nstruct = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    targeted = "--untargeted" not in sys.argv

    OUT.mkdir(parents=True, exist_ok=True)
    parents = pick_parents(n_parents)
    if not parents:
        print("No parent has a bound ligand AND passes the pLDDT gate.")
        print("The loop has nothing to start from - that is a result, not a bug.")
        return 3

    print(f"{len(parents)} parent(s), nstruct={nstruct}, "
          f"{'targeted at ' + str(TARGET_POSITIONS) if targeted else 'full 8 A shell'}\n")
    for p in parents:
        print(f"  {p['design_id']:34} d_base {p['d_base']:.3f}  "
              f"d_oxy {p['d_oxy']:.3f}  dih {p['dih_oxy']:6.1f}  "
              f"pLDDT {p['plddt']:.3f}")

    site = official.CatalyticSite(base_resseq=BASE_IDX, oxy_resseq=OXY_IDX,
                                  chain="A")
    results: list[dict] = []
    t_start = time.time()

    for p in parents:
        pid = p["design_id"]
        parent_pdb = ROOT / "outputs/insurance/folds" / pid / "folded.pdb"
        if not parent_pdb.exists():
            print(f"\n  {pid}: folded.pdb missing - skipped")
            continue
        print(f"\n=== designing on {pid} ===")
        designs = run_design(parent_pdb, f"loop_{pid}", nstruct, targeted)
        print(f"  {len(designs)} design(s) produced")

        parent_seq = chain_sequence(parent_pdb)
        for d in designs:
            seq = chain_sequence(d)
            muts = [f"{parent_seq[i]}{i+1}{seq[i]}"
                    for i in range(min(len(parent_seq), len(seq)))
                    if parent_seq[i] != seq[i]]
            if seq[BASE_IDX - 1] != "D" or seq[OXY_IDX - 1] != "Q":
                print(f"    {d.stem}: catalytic pair lost - skipped")
                continue
            w = OUT / f"{pid}__{d.stem}"
            w.mkdir(parents=True, exist_ok=True)
            try:
                r = boltz_predict(seq, w, seed=C.BOLTZ_SEEDS[0],
                                  ligand_ccd=C.LIGAND_CCD)
                folded = next((x for x in r.structures), None)
                if folded is None:
                    print(f"    {d.stem}: FOLD FAILED")
                    continue
                fp = (folded if folded.suffix == ".pdb"
                      else pipeline.cif_to_pdb(folded, w / "folded.pdb"))
                raw = official.measure_raw(str(fp), site=site)
                conf = r.confidence or {}
                gap = objective.dihedral_band_gap(raw.dihedral12_oxy)
                is_bound = bound(raw.d_base, raw.d_oxy)
                plddt = conf.get("complex_plddt") or 0.0
                sc = official.score(raw, strict=False)
                rec = {"parent": pid, "design": d.stem, "sequence": seq,
                       "n_mutations": len(muts), "mutations": muts,
                       "d_base": raw.d_base, "d_oxy": raw.d_oxy,
                       "dih_oxy": raw.dihedral12_oxy, "gap": gap,
                       "bound": is_bound, "plddt": plddt,
                       "ipde": conf.get("complex_ipde"),
                       "measurable": plddt >= C.PLDDT_MIN,
                       "subscores": sc["subscores"],
                       "badness": objective.badness(raw)}
                results.append(rec)
                flag = ""
                if not is_bound:
                    flag = "  ligand unbound -> geometry is noise"
                elif plddt < C.PLDDT_MIN:
                    flag = f"  pLDDT {plddt:.3f} < {C.PLDDT_MIN}, reserve"
                elif gap == 0.0:
                    flag = "  <<< BOUND, MEASURABLE, DIHEDRAL IN BAND"
                print(f"    {d.stem}: {len(muts)} mut  d_base {raw.d_base:6.3f}  "
                      f"d_oxy {raw.d_oxy:6.3f}  dih {raw.dihedral12_oxy:6.1f}  "
                      f"gap {gap:5.1f}  pLDDT {plddt:.3f}{flag}", flush=True)
            except Exception as exc:                            # noqa: BLE001
                print(f"    {d.stem}: ERROR {type(exc).__name__}: {str(exc)[:110]}")
            (OUT / "loop.json").write_text(json.dumps(results, indent=2,
                                                      default=str))

    print("\n" + "=" * 78)
    good = [r for r in results if r["bound"] and r["measurable"]]
    print(f"{len(results)} refolded;  bound + measurable: {len(good)}")
    if good:
        ranked = objective.rank_batch([r["badness"] for r in good])
        order = sorted(range(len(good)), key=lambda i: -ranked[i]["dense_rank"])
        print(f"\n  {'design':40}{'d_oxy':>8}{'dih':>8}{'gap':>7}{'rank':>7}")
        for i in order[:10]:
            r = good[i]
            print(f"  {r['parent'][:20]+'/'+r['design']:40}{r['d_oxy']:8.3f}"
                  f"{r['dih_oxy']:8.1f}{r['gap']:7.1f}"
                  f"{ranked[i]['dense_rank']:7.3f}")
        inband = [r for r in good if r["gap"] == 0.0]
        if inband:
            print(f"\n  *** {len(inband)} design(s) with a bound ligand AND the "
                  "dihedral in band - run score_finalist.py on these ***")
        else:
            print("\n  None reached the dihedral band. Rosetta can BUILD the "
                  "geometry\n  (allcst ~0.01) but the sequence does not make "
                  "Boltz reproduce it.")
    else:
        print("  Nothing survived the bound-ligand and pLDDT filters.")
    print(f"\nwall {(time.time()-t_start)/60:.1f} min -> {OUT/'loop.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
