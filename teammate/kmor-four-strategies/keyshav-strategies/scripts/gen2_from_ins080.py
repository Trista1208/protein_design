"""Generation 2: third mutations on ins080, to stabilise the oxyanion rotamer.

    python scripts/gen2_from_ins080.py            # generate only
    python scripts/gen2_from_ins080.py --fold     # generate and fold

WHY ins080 IS THE NEW PARENT
----------------------------
Five-seed official score of ins080 (F42Y + H103Y):

    base 0.998   oxy 0.377   ipde 0.939   orient_base 0.615
    orient_oxy 0.000   bury 1.000   cons 0.853
    sd_base 0.0734 (HG3.17 reference 0.068)   mean pLDDT 0.9788

Six of seven terms are non-zero and `cons` is essentially at the reference
enzyme's level, so replicate consistency - one of the three blockers that
zeroed the original parent d045 - is solved. Exactly one term remains.

And the remaining gap is small. Holding the other six fixed:

    orient_oxy 0.05  -> official 0.596
    orient_oxy 0.197 -> official 0.665     (seed 1 ALREADY reached 0.197)
    orient_oxy 0.50  -> official 0.717

THE ACTUAL PROBLEM IS ROTAMER STABILITY, NOT GEOMETRY
-----------------------------------------------------
Per-seed oxyanion dihedral on ins080:

    48.2   94.0   95.9   249.8   245.9      in band on 1 of 5 seeds

The rotamer samples three distinct states: ~48 (in the 310-50 crystal band),
~95 (the 50-130 dead zone) and ~248 (the 230-310 dead zone). It is not that
the pocket cannot reach the band - it reaches it one seed in five. It is that
nothing holds it there.

So generation 2 does NOT look for a new geometry. It looks for a third
mutation that makes the ~48 state the one Boltz keeps picking. That is a
packing question: add or remove bulk next to the Gln46 tip so the other two
rotamer wells become unfavourable.

Positions are the rotamer shell minus the two already used (42, 103) and
minus the catalytic pair.
"""
from __future__ import annotations

import itertools
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C, objective, official, pipeline     # noqa: E402
from kemp.runners import RunnerError, boltz_predict, ligandmpnn  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PARENT_PDB = ROOT / "outputs/finalists/ins080_F42Y_H103Y/seed1/folded.pdb"
OUT = ROOT / "outputs/gen2"
BASE_IDX, OXY_IDX, CHAIN = 102, 46, "A"
IN_WINDOW = (2.0, 3.6)

#: ins080 = d045 + F42Y + H103Y
PARENT_SEQ = None          # read from the manifest, not hard-coded

#: rotamer shell minus the positions ins080 already uses
GEN2_POSITIONS = [26, 23, 27, 45, 49, 84, 99]
STERIC_LADDER = ["A", "V", "L", "I", "F", "W", "M", "T"]
OMIT_EVERYWHERE = "C"
TEMPERATURE = 0.3
N_SEQS, N_BATCHES = 8, 2
MAX_EDITS = 1              # ONE further edit; ins080 is already a double


def parent_sequence() -> str:
    man = json.loads((ROOT / "outputs/insurance/variants.json").read_text())
    for v in man:
        if v["design_id"] == "ins080_r103_42_F42Y_H103Y":
            return v["sequence"]
    raise SystemExit("ins080 not found in variants.json")


def main() -> int:
    seq0 = parent_sequence()
    assert seq0[BASE_IDX - 1] == "D" and seq0[OXY_IDX - 1] == "Q"
    OUT.mkdir(parents=True, exist_ok=True)
    n_res = len(seq0)
    print(f"parent ins080 (F42Y+H103Y), {n_res} aa, Asp{BASE_IDX}/Gln{OXY_IDX}")
    print(f"target: make the ~48 deg rotamer the one Boltz keeps "
          f"(currently 1/5 seeds)\n")

    seen: dict[str, dict] = {}

    # --- model proposals, native forbidden so a substitution is forced
    if PARENT_PDB.exists():
        for pos in GEN2_POSITIONS:
            w = OUT / "mpnn" / f"p{pos}"
            w.mkdir(parents=True, exist_ok=True)
            try:
                res = ligandmpnn(
                    PARENT_PDB, w,
                    redesigned_residues=f"{CHAIN}{pos}",
                    omit_per_residue={f"{CHAIN}{pos}":
                                      seq0[pos - 1] + OMIT_EVERYWHERE},
                    n_seqs=N_SEQS, n_batches=N_BATCHES,
                    temperature=TEMPERATURE, seed=37)
            except RunnerError as exc:
                print(f"  mpnn p{pos}: {str(exc)[:90]}")
                continue
            if res.returncode != 0:
                continue
            for _h, s in res.sequences:
                if len(s) != n_res or s[BASE_IDX - 1] != "D" or s[OXY_IDX - 1] != "Q":
                    continue
                d = [i for i in range(n_res) if s[i] != seq0[i]]
                if len(d) != MAX_EDITS or s in seen:
                    continue
                i = d[0]
                seen[s] = {"tag": f"mpnn{pos}", "edit": f"{seq0[i]}{i+1}{s[i]}",
                           "source": "ligandmpnn"}
    else:
        print(f"  (no parent PDB at {PARENT_PDB} - ladder only)")

    n_mpnn = len(seen)
    # --- mechanistic proposals: bulk up / down beside the Gln46 tip
    for pos in GEN2_POSITIONS:
        native = seq0[pos - 1]
        for aa in STERIC_LADDER:
            if aa == native:
                continue
            s = seq0[:pos - 1] + aa + seq0[pos:]
            if s not in seen:
                seen[s] = {"tag": f"ladder{pos}", "edit": f"{native}{pos}{aa}",
                           "source": "steric_ladder"}

    manifest = [{"design_id": f"g2_{m['tag']}_{m['edit']}", "sequence": s,
                 "parent": "ins080_r103_42_F42Y_H103Y", **m}
                for s, m in seen.items()]
    (OUT / "variants.json").write_text(json.dumps(manifest, indent=2))
    print(f"{len(manifest)} gen-2 variants "
          f"({n_mpnn} LigandMPNN, {len(manifest)-n_mpnn} ladder)")
    print(f"wrote {OUT/'variants.json'}")

    if "--fold" not in sys.argv:
        print("\n(generate-only; pass --fold to screen them)")
        return 0

    site = official.CatalyticSite(base_resseq=BASE_IDX, oxy_resseq=OXY_IDX,
                                  chain=CHAIN)
    rows, t0 = [], time.time()
    print(f"\nfolding {len(manifest)} at 1 seed (~41 s each, "
          f"~{len(manifest)*41/60:.0f} min)\n")
    for i, m in enumerate(manifest):
        w = OUT / "folds" / m["design_id"]
        w.mkdir(parents=True, exist_ok=True)
        try:
            r = boltz_predict(m["sequence"], w, seed=C.BOLTZ_SEEDS[0],
                              ligand_ccd=C.LIGAND_CCD)
            f = next((x for x in r.structures), None)
            if f is None:
                continue
            p = f if f.suffix == ".pdb" else pipeline.cif_to_pdb(f, w / "folded.pdb")
            raw = official.measure_raw(str(p), site=site)
            raw.ipde = (r.confidence or {}).get("complex_ipde")
            plddt = (r.confidence or {}).get("complex_plddt") or 0.0
            gap = objective.dihedral_band_gap(raw.dihedral12_oxy)
            bound = (IN_WINDOW[0] <= (raw.d_base or 99) <= IN_WINDOW[1]
                     and IN_WINDOW[0] <= (raw.d_oxy or 99) <= IN_WINDOW[1])
            rows.append({**m, "d_base": raw.d_base, "d_oxy": raw.d_oxy,
                         "dih_oxy": raw.dihedral12_oxy, "gap": gap,
                         "bound": bound, "plddt": plddt, "ipde": raw.ipde,
                         "badness": objective.badness(raw)})
            flag = ""
            if bound and gap == 0.0 and plddt >= C.PLDDT_MIN:
                flag = "  <<< BOUND + IN BAND + MEASURABLE"
            elif not bound:
                flag = "  unbound"
            print(f"  [{i+1}/{len(manifest)}] {m['design_id']:28} "
                  f"d_base {raw.d_base:6.3f} d_oxy {raw.d_oxy:6.3f} "
                  f"dih {raw.dihedral12_oxy:6.1f} gap {gap:5.1f} "
                  f"pLDDT {plddt:.3f}{flag}", flush=True)
        except Exception as exc:                                # noqa: BLE001
            print(f"  [{i+1}] {m['design_id']}: {type(exc).__name__}")
        (OUT / "screen.json").write_text(json.dumps(rows, indent=2, default=str))

    good = [r for r in rows if r["bound"] and r["plddt"] >= C.PLDDT_MIN]
    hits = [r for r in good if r["gap"] == 0.0]
    print(f"\n{len(rows)} folded; bound+measurable {len(good)}; in band {len(hits)}")
    for r in sorted(good, key=lambda x: x["gap"])[:10]:
        print(f"   {r['design_id']:28} gap {r['gap']:5.1f}  "
              f"d_oxy {r['d_oxy']:6.3f}  pLDDT {r['plddt']:.3f}")
    print(f"\nwall {(time.time()-t0)/60:.1f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
