"""Re-score every finalist per replicate, from the structures on disk.

    python scripts/rescore_all.py [--write]

Any verdict.json written before 05:40 on 22 Nov holds an averaged-geometry
score, which is wrong in both directions:

    ins079  averaged 0.624  but 0/5 replicates score  -> false positive
    ins080  averaged 0.000  but 1/5 replicates scores 0.743 -> erased a result

It re-measures from the folded PDBs rather than trusting the stored numbers,
because the stored replicate records do not carry `angle2_oxy`. Reconstructing
a RawMetrics from them leaves that term MISSING rather than zero, and
`score(strict=False)` then silently drops it and inflates the total — the
LESSONS §1.5 failure. Measuring from the structure keeps every term present.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import official                                      # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIN = ROOT / "outputs/finalists"
SITE = official.CatalyticSite(base_resseq=102, oxy_resseq=46, chain="A")


def main() -> int:
    write = "--write" in sys.argv
    rows = []
    for vd in sorted(FIN.glob("*/verdict.json")):
        d = vd.parent
        pdbs = sorted(d.glob("seed*/folded.pdb"))
        if not pdbs:
            continue
        v = json.loads(vd.read_text())
        raws = [official.measure_raw(str(p), site=SITE) for p in pdbs]
        ipdes = [r.get("ipde") for r in v.get("replicates", [])
                 if r.get("ipde") is not None]
        rep = official.score_replicates(raws, ipde_values=ipdes)
        old = v.get("official_score")
        rows.append({"design": d.name, "n": len(raws), "old": old,
                     "worst": rep["worst_total"], "median": rep["median_total"],
                     "best": rep["best_total"],
                     "nonzero": rep["n_nonzero_replicates"],
                     "per_rep": rep["per_replicate_total"],
                     "zeroed": rep["terms_zeroed_on_any_replicate"]})
        if write:
            v["official_score"] = rep["median_total"]
            v["official_score_worst"] = rep["worst_total"]
            v["official_score_best"] = rep["best_total"]
            v["n_nonzero_replicates"] = rep["n_nonzero_replicates"]
            v["per_replicate_total"] = rep["per_replicate_total"]
            v["terms_zeroed_on_any_replicate"] = rep["terms_zeroed_on_any_replicate"]
            v["averaged_geometry_total_DIAGNOSTIC"] = rep["averaged_geometry_total"]
            v["rescored_note"] = ("re-measured from structures and scored per "
                                  "replicate; the previous value averaged the "
                                  "geometry and scored once")
            vd.write_text(json.dumps(v, indent=2, default=str))

    rows.sort(key=lambda r: -(r["best"] or 0))
    print(f"{'design':34}{'old':>9}{'best':>9}{'median':>9}{'worst':>8}"
          f"{'nz':>5}  zeroed-on-any")
    for r in rows:
        print(f"  {r['design']:32}{str(r['old']):>9}{r['best']:9.4f}"
              f"{r['median']:9.4f}{r['worst']:8.4f}"
              f"{r['nonzero']:>3}/{r['n']}  {','.join(r['zeroed'])}")
    real = [r for r in rows if r["best"] > 0]
    print(f"\n{len(real)}/{len(rows)} design(s) score on at least one replicate")
    for r in real:
        print(f"   {r['design']}: {r['per_rep']}")
    if not write:
        print("\n(dry run; pass --write to update verdict.json)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
