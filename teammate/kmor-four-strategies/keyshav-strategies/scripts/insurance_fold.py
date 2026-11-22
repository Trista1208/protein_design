"""Stage 2: fold and score the d045 insurance family.

One Boltz seed per variant as a SCREEN (~41 s each, measured). Five seeds are
for finalists only - the Consistency term needs them, but spending 5x on a
variant that has not yet shown anything would burn the budget.

Protocol is pinned to the parent's: `msa: empty`, diffusion_samples 1,
recycling 3, no pocket constraint. The MSA question is closed - the server
returns depth 1 for these de novo sequences and the folds are bit-identical
(LOOP_DESIGN 6.1) - so single-sequence is both correct and cheaper.

WHAT IS BEING LOOKED FOR, in priority order:
  1. d_oxy inside 2.0-3.6 A. The parent's best seed is 3.654 - it misses the
     window by 0.054 A on EVERY seed, so `oxy` is a hard 0.000 (LOOP_DESIGN
     6.2). Any variant that lands inside the window has done the single most
     valuable thing available.
  2. oxyanion dihedral out of the 230-310 dead zone. Parent: 233.4-239.8,
     best gap 3.4 deg.
  3. pLDDT >= C.PLDDT_MIN. Below it the geometry is noise, not measurement
     (LESSONS 7.1) - such a variant is UNMEASURABLE, not bad, and is kept in a
     reserve tier rather than deleted, because the organisers accept five
     designs regardless of score.

Nothing here decides a winner. It produces a ranked, compliant pool.
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
VARIANTS = ROOT / "outputs/insurance/variants.json"
OUT = ROOT / "outputs/insurance/folds"
BASE_IDX, OXY_IDX = 102, 46
SEED = C.BOLTZ_SEEDS[0]

#: the parent, measured over 5 seeds (outputs/refold_d045/verdict.json)
PARENT = {"design_id": "rfd4_design_045 (parent)", "d_base": 2.864,
          "d_oxy": 3.654, "dih_o": 234.4, "gap": 4.4, "plddt": 0.967,
          "note": "best single seed; oxy=0.000 on all five"}


def main() -> int:
    variants = json.loads(VARIANTS.read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"{len(variants)} variants, 1 seed each (~41 s) "
          f"-> ~{len(variants) * 41 / 60:.0f} min")
    print(f"parent best seed: d_oxy {PARENT['d_oxy']} (window closes 3.6), "
          f"dihedral {PARENT['dih_o']} (band edge 230)\n")

    site = official.CatalyticSite(base_resseq=BASE_IDX, oxy_resseq=OXY_IDX,
                                  chain="A")
    rows: list[dict] = []
    t_start = time.time()
    for i, v in enumerate(variants):
        did, seq = v["design_id"], v["sequence"]
        w = OUT / did
        w.mkdir(parents=True, exist_ok=True)
        row = {**{k: v[k] for k in ("design_id", "edits", "n_edits", "tag")},
               "sequence": seq}
        try:
            r = boltz_predict(seq, w, seed=SEED, ligand_ccd=C.LIGAND_CCD)
            folded = next((p for p in r.structures), None)
            if folded is None:
                row["status"] = "fold_failed"
                rows.append(row)
                print(f"  [{i+1}/{len(variants)}] {did:38} FOLD FAILED")
                continue
            pdb = (folded if folded.suffix == ".pdb"
                   else pipeline.cif_to_pdb(folded, w / "folded.pdb"))
            raw = official.measure_raw(str(pdb), site=site)
            conf = r.confidence or {}
            plddt = conf.get("complex_plddt")
            ipde = conf.get("complex_ipde")
            # measure_raw() does NOT set ipde - that normally happens in
            # combine_replicates(). Without this line objective.badness()
            # returns ipde=None for every design and the dense rank silently
            # drops a 0.20-weight term. The degeneracy alarm caught exactly
            # this on the first full screen (fatal=['ipde']).
            raw.ipde = ipde
            gap = objective.dihedral_band_gap(raw.dihedral12_oxy)
            row.update({
                "status": "ok",
                "d_base": raw.d_base, "d_oxy": raw.d_oxy,
                "dih_base": raw.dihedral12_base, "dih_oxy": raw.dihedral12_oxy,
                "gap_dih_oxy": gap, "plddt": plddt, "ipde": ipde,
                "measurable": (plddt is not None and plddt >= C.PLDDT_MIN),
                "d_oxy_in_window": (raw.d_oxy is not None
                                    and 2.0 <= raw.d_oxy <= 3.6),
                "dih_in_band": (gap == 0.0 if gap is not None else False),
                "badness": objective.badness(raw),
            })
            flag = ""
            if row["d_oxy_in_window"]:
                flag += "  <<< d_oxy IN WINDOW"
            if row["dih_in_band"]:
                flag += "  <<< DIHEDRAL IN BAND"
            if not row["measurable"]:
                flag += f"  (pLDDT {plddt:.3f} < {C.PLDDT_MIN}, reserve tier)"
            print(f"  [{i+1}/{len(variants)}] {did:38} "
                  f"d_oxy {raw.d_oxy:6.3f}  dih {raw.dihedral12_oxy:6.1f}  "
                  f"gap {gap:5.1f}  pLDDT {plddt:.3f}{flag}", flush=True)
        except Exception as exc:                                # noqa: BLE001
            row["status"] = f"error: {type(exc).__name__}: {str(exc)[:160]}"
            print(f"  [{i+1}/{len(variants)}] {did:38} ERROR {row['status'][:70]}")
        rows.append(row)
        (OUT.parent / "screen.json").write_text(json.dumps(rows, indent=2,
                                                           default=str))

    ok = [r for r in rows if r.get("status") == "ok"]
    if not ok:
        print("\nNo variant folded. Nothing to rank.")
        return 3

    # Degeneracy check: if every variant reports the same geometry we are
    # measuring the parent, not the variants (LESSONS 2.5).
    cp = objective.campaign_checkpoint([r["badness"] for r in ok],
                                       after=min(10, len(ok)))
    print("\n" + "=" * 78)
    print(f"degeneracy checkpoint: {'PASS' if cp['passed'] else 'FAIL'}"
          + (f"  fatal={cp['fatal']}" if cp.get("fatal") else ""))

    measurable = [r for r in ok if r["measurable"]]
    reserve = [r for r in ok if not r["measurable"]]
    print(f"folded {len(ok)}/{len(rows)}   measurable (pLDDT >= "
          f"{C.PLDDT_MIN}): {len(measurable)}   reserve: {len(reserve)}")

    hits_w = [r for r in measurable if r["d_oxy_in_window"]]
    hits_b = [r for r in measurable if r["dih_in_band"]]
    print(f"d_oxy inside 2.0-3.6 A : {len(hits_w)}   "
          f"dihedral in band: {len(hits_b)}")

    pool = measurable or ok
    ranked = objective.rank_batch([r["badness"] for r in pool])
    order = sorted(range(len(pool)), key=lambda i: -ranked[i]["dense_rank"])
    print(f"\ntop {min(12, len(order))} by dense rank "
          "(internal surrogate, NOT a score):")
    print(f"  {'design':40} {'d_oxy':>7} {'dih':>7} {'gap':>6} "
          f"{'pLDDT':>6} {'rank':>6}")
    for i in order[:12]:
        r = pool[i]
        print(f"  {r['design_id']:40} {r['d_oxy']:7.3f} {r['dih_oxy']:7.1f} "
              f"{r['gap_dih_oxy']:6.1f} {r['plddt']:6.3f} "
              f"{ranked[i]['dense_rank']:6.3f}")
    print(f"\n  parent for comparison:                   "
          f"{PARENT['d_oxy']:7.3f} {PARENT['dih_o']:7.1f} "
          f"{PARENT['gap']:6.1f} {PARENT['plddt']:6.3f}")

    (OUT.parent / "screen_ranked.json").write_text(json.dumps({
        "protocol": objective.protocol_stamp(replicates=1),
        "parent": PARENT,
        "degeneracy": cp,
        "n_folded": len(ok), "n_measurable": len(measurable),
        "n_d_oxy_in_window": len(hits_w), "n_dihedral_in_band": len(hits_b),
        "ranked": [{**pool[i], "dense_rank": ranked[i]["dense_rank"],
                    "worst_term": ranked[i]["worst_term"]} for i in order],
    }, indent=2, default=str))
    print(f"\nwall {(time.time()-t_start)/60:.1f} min  "
          f"-> {OUT.parent / 'screen_ranked.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
