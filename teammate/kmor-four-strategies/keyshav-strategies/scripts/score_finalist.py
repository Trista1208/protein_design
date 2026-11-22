"""Full 5-seed official score for one design. The only verdict that counts.

    python scripts/score_finalist.py <design_id|sequence> [label]

Five Boltz-2 seeds is the organisers' protocol ("measured on only 5 replicates
per sequence", slide 7), and it is the only way to compute `cons` - the
standard deviation of d_base across replicates. A 1-seed screen can rank, but
it cannot score.

Reads design_ids from the insurance manifests, or takes a raw sequence.

WHAT TO WATCH, beyond the total
-------------------------------
* `cons`. The parent d045 has FOUR excellent seeds (d_base sd ~0.075, matching
  HG3.17's 0.068) and one that loses the ligand entirely (d_base 5.272,
  d_oxy 9.654). That single seed takes sd_base to 1.047 and zeroes `cons` by
  itself. Watch the per-seed table, not the aggregate.
* The spread of the oxyanion dihedral. A variant whose rotamer is reproducible
  is a design; one that scatters across seeds is noise wearing a good number.
* pLDDT. Below C.PLDDT_MIN the geometry is unmeasurable rather than bad
  (LESSONS 7.1), and this script says so rather than ranking it.
"""
from __future__ import annotations

import json
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C, objective, official, pipeline     # noqa: E402
from kemp.runners import boltz_predict                          # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BASE_IDX, OXY_IDX = 102, 46
OUT = ROOT / "outputs/finalists"
MANIFESTS = [ROOT / "outputs/insurance/variants.json",
             ROOT / "outputs/insurance/combos/combos.json"]


def find_sequence(key: str) -> tuple[str, str]:
    if set(key.upper()) <= set("ACDEFGHIKLMNPQRSTVWY") and len(key) > 30:
        return key.upper(), "raw_sequence"
    for m in MANIFESTS:
        if not m.exists():
            continue
        for row in json.loads(m.read_text()):
            if row.get("design_id") == key:
                return row["sequence"], str(m.name)
    raise SystemExit(f"design_id {key!r} not found in {[m.name for m in MANIFESTS]}")


def main() -> int:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    key = sys.argv[1]
    label = sys.argv[2] if len(sys.argv) > 2 else key
    seq, src = find_sequence(key)
    assert seq[BASE_IDX - 1] == "D" and seq[OXY_IDX - 1] == "Q", \
        "catalytic pair not preserved in this sequence"

    w_root = OUT / label
    w_root.mkdir(parents=True, exist_ok=True)
    site = official.CatalyticSite(base_resseq=BASE_IDX, oxy_resseq=OXY_IDX,
                                  chain="A")
    print(f"{label}  ({len(seq)} aa, from {src})")
    print(f"{C.REPLICATES} seeds, msa: empty, Asp{BASE_IDX}/Gln{OXY_IDX}\n")

    raws, ipdes, plddts, rows = [], [], [], []
    t0 = time.time()
    for seed in C.BOLTZ_SEEDS[:C.REPLICATES]:
        w = w_root / f"seed{seed}"
        w.mkdir(parents=True, exist_ok=True)
        r = boltz_predict(seq, w, seed=seed, ligand_ccd=C.LIGAND_CCD)
        folded = next((p for p in r.structures), None)
        if folded is None:
            print(f"  seed {seed}: FOLD FAILED")
            continue
        pdb = (folded if folded.suffix == ".pdb"
               else pipeline.cif_to_pdb(folded, w / "folded.pdb"))
        raw = official.measure_raw(str(pdb), site=site)
        raws.append(raw)
        conf = r.confidence or {}
        if conf.get("complex_ipde") is not None:
            ipdes.append(conf["complex_ipde"])
        if conf.get("complex_plddt") is not None:
            plddts.append(float(conf["complex_plddt"]))
        gap = objective.dihedral_band_gap(raw.dihedral12_oxy)
        rows.append({"seed": seed, "d_base": raw.d_base, "d_oxy": raw.d_oxy,
                     "dih_oxy": raw.dihedral12_oxy, "gap": gap,
                     "plddt": conf.get("complex_plddt"),
                     "ipde": conf.get("complex_ipde")})
        print(f"  seed {seed}: d_base {raw.d_base:6.3f}  d_oxy {raw.d_oxy:6.3f}  "
              f"dih_oxy {raw.dihedral12_oxy:6.1f}  gap {gap:5.1f}  "
              f"pLDDT {conf.get('complex_plddt', 0):.3f}", flush=True)

    if not raws:
        print("\nnothing folded")
        return 3

    comb = official.combine_replicates(raws, ipde_values=ipdes)
    sc = official.score(comb, strict=True)
    # PER-REPLICATE is the number that counts. Averaging the geometry and
    # scoring once invented 0.624 for a design where 0/5 replicates scored,
    # and erased a real 0.743 on another. See official.score_replicates().
    rep = official.score_replicates(raws, ipde_values=ipdes)
    mean_plddt = statistics.fmean(plddts) if plddts else None
    measurable = mean_plddt is not None and mean_plddt >= C.PLDDT_MIN

    print("\n" + "=" * 72)
    print(f"{label}   {len(raws)} replicates")
    print(f"  PER-REPLICATE totals : {rep['per_replicate_total']}")
    print(f"  worst / median / best: {rep['worst_total']} / "
          f"{rep['median_total']} / {rep['best_total']}")
    print(f"  replicates scoring >0: {rep['n_nonzero_replicates']}/"
          f"{rep['n_replicates']}   zeroed on any: "
          f"{rep['terms_zeroed_on_any_replicate']}")
    print(f"  (averaged-geometry total, DIAGNOSTIC ONLY: "
          f"{rep['averaged_geometry_total']})")
    print("  subscores: " + "  ".join(
        f"{k}={v:.3f}" if v is not None else f"{k}=None"
        for k, v in sc["subscores"].items()))
    print(f"  official_score = {sc.get('official_score')}   "
          f"limiting = {sc.get('limiting_term')}")
    print(f"  sd_base = {comb.sd_base if comb.sd_base is None else round(comb.sd_base,4)}"
          f"   (HG3.17 reference 0.068)")
    if comb.sd_base is not None and len(raws) >= 3:
        db = sorted(r.d_base for r in raws if r.d_base is not None)
        trimmed = statistics.stdev(db[:-1]) if len(db) > 2 else None
        if trimmed is not None:
            print(f"  sd_base excluding the worst seed = {trimmed:.4f} "
                  "(DIAGNOSTIC ONLY - the official term allows no exclusion)")
    print(f"  mean pLDDT = {mean_plddt:.4f}  "
          f"({'measurable' if measurable else f'BELOW {C.PLDDT_MIN} - unmeasurable, not bad'})")

    dihs = [r["dih_oxy"] for r in rows if r["dih_oxy"] is not None]
    if len(dihs) >= 2:
        print(f"  oxyanion dihedral spread: {min(dihs):.1f} - {max(dihs):.1f} deg"
              f"   in band on {sum(1 for r in rows if r['gap']==0)}/{len(rows)} seeds")

    (w_root / "verdict.json").write_text(json.dumps({
        "label": label, "sequence": seq, "source": src,
        "protocol": objective.protocol_stamp(replicates=len(raws)),
        "replicates": rows, "subscores": sc["subscores"],
        "official_score": rep["median_total"],
        "official_score_worst": rep["worst_total"],
        "official_score_best": rep["best_total"],
        "n_nonzero_replicates": rep["n_nonzero_replicates"],
        "per_replicate_total": rep["per_replicate_total"],
        "terms_zeroed_on_any_replicate": rep["terms_zeroed_on_any_replicate"],
        "averaged_geometry_total_DIAGNOSTIC": rep["averaged_geometry_total"],
        "limiting_term": sc.get("limiting_term"),
        "sd_base": comb.sd_base, "mean_plddt": mean_plddt,
        "measurable": measurable,
        "badness": objective.badness(comb),
    }, indent=2, default=str))
    print(f"\n  wall {(time.time()-t0)/60:.1f} min -> {w_root/'verdict.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
