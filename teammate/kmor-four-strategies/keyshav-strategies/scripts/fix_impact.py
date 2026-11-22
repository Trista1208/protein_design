"""Which of our scores would differ WITHOUT each bug fix?

    python scripts/fix_impact.py

Answers a direct question with a measurement instead of an assertion: for each
design in the current top 20, recompute the score the OLD way and report the
delta.

FIXES TESTED
  A  averaged geometry  - mean the geometry over 5 seeds, score once, versus
     scoring each replicate and aggregating the scores
  B  d_oxy target 2.40  - the retired target, taken from a constraint file for
     a DIFFERENT acceptor atom, versus 3.10 measured from crystals
  C  linear vs circular mean of the dihedrals, which only bites when geometry
     is averaged at all

Everything is re-measured from the folded structures; no stored numbers.
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import official as O                                   # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIN = ROOT / "outputs/finalists"
SITE = O.CatalyticSite(base_resseq=102, oxy_resseq=46, chain="A")


def linear_mean(vals):
    return sum(vals) / len(vals)


def avg_raw(raws, circular: bool):
    """Average the geometry across replicates, the way the old code did."""
    f = O._circular_mean if circular else linear_mean
    def m(attr, circ=False):
        v = [getattr(r, attr) for r in raws if getattr(r, attr) is not None]
        if not v:
            return None
        return (f(v) if circ else linear_mean(v))
    kw = dict(
        d_base=m("d_base"), d_oxy=m("d_oxy"),
        angle2_base=m("angle2_base"), angle2_oxy=m("angle2_oxy"),
        dihedral12_base=m("dihedral12_base", True),
        dihedral12_oxy=m("dihedral12_oxy", True),
        ca_neighbors_mean=m("ca_neighbors_mean"),
        n_replicates=len(raws))
    db = [r.d_base for r in raws if r.d_base is not None]
    if len(db) > 1:
        mu = linear_mean(db)
        kw["sd_base"] = math.sqrt(sum((x - mu) ** 2 for x in db) / (len(db) - 1))
    return O.RawMetrics(**kw)


def main() -> int:
    top = json.loads((ROOT / "outputs/top20.json").read_text())
    print(f"{'design':32} {'now':>7} {'avg-geom':>9} {'d_oxy=2.40':>11} "
          f"{'linear-mean':>12}")
    print("-" * 76)
    n_avg = n_oxy = n_lin = 0
    for r in top:
        did = r["id"]
        d = FIN / did
        vd = d / "verdict.json"
        if not vd.exists():
            continue
        v = json.loads(vd.read_text())
        pdbs = sorted(d.glob("seed*/folded.pdb"))
        if not pdbs:
            continue
        raws = [O.measure_raw(str(p), site=SITE) for p in pdbs]
        ip = [x.get("ipde") for x in v.get("replicates", [])
              if x.get("ipde") is not None]

        now = O.score_replicates(raws, ipde_values=ip)["median_total"]

        # A: averaged geometry, scored once (circular mean, so this isolates
        #    the averaging itself from the circular-mean bug)
        ipm = linear_mean(ip) if ip else None
        a = avg_raw(raws, circular=True)
        a.ipde = ipm
        s_avg = O.score(a, strict=False)["geometric_total"] or 0.0

        # C: same, but with the LINEAR mean of the dihedrals
        c = avg_raw(raws, circular=False)
        c.ipde = ipm
        s_lin = O.score(c, strict=False)["geometric_total"] or 0.0

        # B: per-replicate, but with the retired 2.40 A oxyanion target
        old = dict(O.TARGETS["d_oxy"])
        O.TARGETS["d_oxy"] = {**old, "target": 2.40}
        try:
            s_oxy = O.score_replicates(raws, ipde_values=ip)["median_total"]
        finally:
            O.TARGETS["d_oxy"] = old

        n_avg += abs(s_avg - now) > 0.005
        n_oxy += abs(s_oxy - now) > 0.005
        n_lin += abs(s_lin - s_avg) > 0.005
        print(f"{did:32} {now:7.4f} {s_avg:9.4f} {s_oxy:11.4f} {s_lin:12.4f}")

    print("-" * 76)
    print(f"differ from the current score:  averaged-geometry {n_avg}/{len(top)}"
          f",  d_oxy=2.40 {n_oxy}/{len(top)}")
    print(f"linear vs circular mean differ (only when geometry is averaged): "
          f"{n_lin}/{len(top)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
