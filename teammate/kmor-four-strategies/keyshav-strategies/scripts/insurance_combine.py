"""Combine the single mutants whose defects are complementary.

THE OBSERVATION THIS IS BUILT ON
--------------------------------
The 84-variant screen separated the two oxyanion defects onto DIFFERENT
positions. Scored with official.py:

    variant   s_oxy    s_dih   torsion band
    parent    0.000    0.000   dead zone (234.4 deg)
    F26Y      0.000    0.727   310-50, the CRYSTAL band (16.4 deg)
    H103F     0.817    0.000   dead zone
    F42Y      0.797    0.000   dead zone
    L49R      0.606    0.000   dead zone
    V27I      0.713    0.360   130-230 (the other band)

F26Y is the only variant that reaches the band both crystals use, and it
reaches it well (0.727 against 4BS0's 0.954). Its failure is the distance -
the ligand drifts to 9.4 A. H103F and F42Y fix exactly that and nothing else.

So the single mutants are not competing solutions, they are complementary
halves. If their effects compose at all, F26Y+H103F has s_oxy ~0.8 AND
s_dih ~0.73, which would be the first design with no zero in the oxyanion
terms - worth 0.28 of the official score that is currently 0.

NO CLAIM IS MADE THAT THEY COMPOSE. Second-shell effects routinely do not, and
a double mutant can break the fold outright. This is a 4-minute experiment that
answers it, which is cheaper than reasoning about it.

Each double is also exactly what FuncLib would enumerate at these positions, so
this is a same-night preview of that run rather than a replacement for it.
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
PARENT_SEQ = ("SSFLELMHKAVSAAKQGDFTSAISDFVEACNLLLPDDKKDDFEAMQKALLSTLLAHHDQDRAKDALTTLLTVTLILFE"
              "KHHDDVVELMQKHYDHFHDMVSDDHKIVTEYLKGELDYEDYLTSLVHNRGDALNLSEAEITEIYEEAKKMLS")
BASE_IDX, OXY_IDX = 102, 46
OUT = ROOT / "outputs/insurance/combos"

#: (position, native, mutant, what it fixed on its own)
TORSION_FIXERS = [(26, "F", "Y", "s_dih 0.727, crystal band")]
DISTANCE_FIXERS = [(103, "H", "F", "s_oxy 0.817"),
                   (42, "F", "Y", "s_oxy 0.797"),
                   (49, "L", "R", "s_oxy 0.606"),
                   (27, "V", "I", "s_oxy 0.713 + s_dih 0.360")]

#: V27I reached both terms alone but folds at pLDDT 0.863, under the gate.
#: Pairing it with a distance fixer is a long shot; pairing the two
#: both-term-movers (26 + 27) is the more interesting one.
EXTRA_PAIRS = [((27, "V", "I"), (103, "H", "F")),
               ((27, "V", "I"), (42, "F", "Y"))]


def mutate(seq: str, muts: list[tuple[int, str, str]]) -> str:
    s = list(seq)
    for pos, native, mut in muts:
        assert s[pos - 1] == native, f"position {pos} is {s[pos-1]}, not {native}"
        s[pos - 1] = mut
    return "".join(s)


def main() -> int:
    combos: list[tuple[str, list[tuple[int, str, str]]]] = []
    for t in TORSION_FIXERS:
        for d in DISTANCE_FIXERS:
            combos.append((f"{t[1]}{t[0]}{t[2]}_{d[1]}{d[0]}{d[2]}", [t, d]))
    for a, b in EXTRA_PAIRS:
        combos.append((f"{a[1]}{a[0]}{a[2]}_{b[1]}{b[0]}{b[2]}", [a, b]))

    OUT.mkdir(parents=True, exist_ok=True)
    site = official.CatalyticSite(base_resseq=BASE_IDX, oxy_resseq=OXY_IDX,
                                  chain="A")
    print(f"{len(combos)} double mutants, 1 seed each (~41 s)\n")
    print("  reference: parent s_oxy 0.000 s_dih 0.000 | "
          "F26Y s_dih 0.727 | H103F s_oxy 0.817\n")

    rows = []
    for i, (tag, muts) in enumerate(combos):
        did = f"combo_{tag}"
        seq = mutate(PARENT_SEQ, muts)
        w = OUT / did
        w.mkdir(parents=True, exist_ok=True)
        row = {"design_id": did, "edits": [f"{n}{p}{m}" for p, n, m in muts],
               "sequence": seq}
        try:
            r = boltz_predict(seq, w, seed=C.BOLTZ_SEEDS[0],
                              ligand_ccd=C.LIGAND_CCD)
            folded = next((p for p in r.structures), None)
            if folded is None:
                row["status"] = "fold_failed"
                print(f"  [{i+1}/{len(combos)}] {did:26} FOLD FAILED")
                rows.append(row)
                continue
            pdb = (folded if folded.suffix == ".pdb"
                   else pipeline.cif_to_pdb(folded, w / "folded.pdb"))
            raw = official.measure_raw(str(pdb), site=site)
            conf = r.confidence or {}
            sc = official.score(raw, strict=False)
            subs = sc["subscores"]
            row.update({
                "status": "ok", "d_oxy": raw.d_oxy, "d_base": raw.d_base,
                "dih_oxy": raw.dihedral12_oxy,
                "gap": objective.dihedral_band_gap(raw.dihedral12_oxy),
                "s_oxy": subs.get("oxy"), "s_orient_oxy": subs.get("orient_oxy"),
                "s_base": subs.get("base"), "plddt": conf.get("complex_plddt"),
                "ipde": conf.get("complex_ipde"),
            })
            both = (subs.get("oxy") or 0) > 0 and (subs.get("orient_oxy") or 0) > 0
            print(f"  [{i+1}/{len(combos)}] {did:26} "
                  f"d_oxy {raw.d_oxy:6.3f} s_oxy {subs.get('oxy') or 0:.3f}  "
                  f"dih {raw.dihedral12_oxy:6.1f} s_or {subs.get('orient_oxy') or 0:.3f}  "
                  f"pLDDT {conf.get('complex_plddt', 0):.3f}"
                  + ("   <<< BOTH NON-ZERO" if both else ""), flush=True)
        except Exception as exc:                                # noqa: BLE001
            row["status"] = f"error: {type(exc).__name__}: {str(exc)[:140]}"
            print(f"  [{i+1}/{len(combos)}] {did:26} ERROR")
        rows.append(row)
        (OUT / "combos.json").write_text(json.dumps(rows, indent=2, default=str))

    ok = [r for r in rows if r.get("status") == "ok"]
    both = [r for r in ok if (r.get("s_oxy") or 0) > 0
            and (r.get("s_orient_oxy") or 0) > 0]
    print("\n" + "=" * 70)
    print(f"folded {len(ok)}/{len(rows)};  both oxyanion terms non-zero: {len(both)}")
    for r in sorted(both, key=lambda x: -(x["s_oxy"] * x["s_orient_oxy"])):
        print(f"   {r['design_id']:26} s_oxy {r['s_oxy']:.3f}  "
              f"s_orient_oxy {r['s_orient_oxy']:.3f}  pLDDT {r['plddt']:.3f}")
    if not both:
        print("   none - the two effects do not compose at these positions.")
        print("   That is a real answer: FuncLib's combinatorial search over")
        print("   26/27/42/103 is then the right next move, not more pairs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
