"""Does the d_oxy target choice change WHICH designs we submit?

    python scripts/rank_sensitivity.py

The submitted FASTA is a set of 20 sequences; the organisers re-fold and
re-score every one themselves. So our own absolute scores never reach them.
The only thing our scoring decides is WHICH 20 of the 30 scoring designs go
in, and in what order we'd defend them.

This re-measures every scoring design from its structures under both oxyanion
targets - 3.10 A (measured from 4BS0 / 7K4Z crystals, what we use) and the
retired 2.40 A (from a constraint file for a DIFFERENT acceptor atom) - and
compares the resulting top-20 sets.

Note the two targets share the same WINDOW (2.0-3.6 A). The window decides
whether a replicate scores at all; the target only decides where inside it the
subscore peaks. So reproducibility (n/5) should be identical under both, and
only the median tiebreak can move. That prediction is checked, not assumed.
"""
from __future__ import annotations
import json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import official as O                                   # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIN = ROOT / "outputs/finalists"
SITE = O.CatalyticSite(base_resseq=102, oxy_resseq=46, chain="A")


def rank_all(target: float) -> list[dict]:
    old = dict(O.TARGETS["d_oxy"])
    O.TARGETS["d_oxy"] = {**old, "target": target}
    out = []
    try:
        for vd in FIN.glob("*/verdict.json"):
            try:
                v = json.loads(vd.read_text())
            except Exception:                                    # noqa: BLE001
                continue
            pdbs = sorted(vd.parent.glob("seed*/folded.pdb"))
            if not pdbs or not v.get("sequence"):
                continue
            raws = [O.measure_raw(str(p), site=SITE) for p in pdbs]
            ip = [x.get("ipde") for x in v.get("replicates", [])
                  if x.get("ipde") is not None]
            rep = O.score_replicates(raws, ipde_values=ip)
            if rep["n_nonzero_replicates"] == 0:
                continue
            out.append({"id": vd.parent.name,
                        "nz": rep["n_nonzero_replicates"],
                        "n": len(rep["per_replicate_total"]),
                        "med": rep["median_total"],
                        "best": rep["best_total"]})
    finally:
        O.TARGETS["d_oxy"] = old
    out.sort(key=lambda r: (-r["nz"], -r["med"], -r["best"]))
    return out


def main() -> int:
    a = rank_all(3.10)
    b = rank_all(2.40)
    print(f"{len(a)} scoring designs under 3.10 A, {len(b)} under 2.40 A\n")

    # prediction check: reproducibility must be identical (same window)
    ra = {r["id"]: r["nz"] for r in a}
    rb = {r["id"]: r["nz"] for r in b}
    diff = {k for k in ra if ra[k] != rb.get(k)}
    print(f"reproducibility (n/5) differs for {len(diff)} design(s): "
          f"{sorted(diff) if diff else 'none - as predicted'}\n")

    sa = [r["id"] for r in a[:20]]
    sb = [r["id"] for r in b[:20]]
    print(f"{'rank':>4}  {'at 3.10 A (ours)':34} {'med':>7}   "
          f"{'at 2.40 A (retired)':34} {'med':>7}")
    for i in range(20):
        x, y = a[i], b[i]
        mark = "  <-- different design" if x["id"] != y["id"] else ""
        print(f"{i+1:4d}  {x['id']:34} {x['med']:7.4f}   "
              f"{y['id']:34} {y['med']:7.4f}{mark}")

    only_a, only_b = set(sa) - set(sb), set(sb) - set(sa)
    print()
    print(f"SET of submitted 20: {'IDENTICAL' if not only_a and not only_b else 'DIFFERENT'}")
    if only_a:
        print(f"  in ours only : {sorted(only_a)}")
    if only_b:
        print(f"  in 2.40 only : {sorted(only_b)}")
    moved = sum(1 for i in range(20) if sa[i] != sb[i])
    print(f"positions that change order: {moved}/20")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
