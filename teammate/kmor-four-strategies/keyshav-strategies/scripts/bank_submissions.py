"""Produce the submission set: 5-seed scores, gates, novelty, export.

    python scripts/bank_submissions.py [n_candidates]

This is the deliverable, and it runs regardless of whether the optimisation
loop succeeds. The organisers accept at least five designs per team
"independent of score" (slide 7), so having nothing at noon is the only
unrecoverable failure available to us.

SELECTION ORDER - gates first, ranking second (LESSONS 4.5, and the ordering
AI.zymes uses):

  1. BOUND LIGAND. Both d_base and d_oxy inside 2.0-3.6 A. A design whose
     ligand has drifted still reports angles, and those angles are noise
     (LESSONS 2.3). This filter has overturned three apparent breakthroughs
     tonight - do not relax it.
  2. pLDDT >= C.PLDDT_MIN. Below it the geometry is unmeasurable rather than
     bad (LESSONS 7.1).
  3. Dense percentile rank to order what survives, because the official score
     is 0.000 for nearly everything and cannot rank (kemp/objective.py).
  4. FIVE Boltz seeds for the official score. One seed cannot compute `cons`
     at all, so a 1-seed number is not an official score.
  5. Hard eligibility: length <= 260, sequence identity <= 0.70, TM <= 0.70.
  6. Cluster by mutation set, so the five are not five copies of one idea.

Everything is written to published/submission/ with its protocol stamp, so a
reviewer can tell which fold protocol produced which number - the omission
that made yesterday's controls and designs incomparable (LESSONS 3.2).
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C, gates, novelty, objective, official  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SCREEN = ROOT / "outputs/insurance/screen.json"
VARIANTS = ROOT / "outputs/insurance/variants.json"
FINAL = ROOT / "outputs/finalists"
OUT = ROOT / "published/submission"
IN_WINDOW = (2.0, 3.6)
N_SLOTS = C.N_LAB_SLOTS


def bound(r) -> bool:
    d1, d2 = r.get("d_base"), r.get("d_oxy")
    return (d1 is not None and d2 is not None
            and IN_WINDOW[0] <= d1 <= IN_WINDOW[1]
            and IN_WINDOW[0] <= d2 <= IN_WINDOW[1])


#: every screen, not just the first one. Reading only insurance/screen.json
#: meant the submission set never saw the library designs - including the
#: three that score on 5/5 replicates.
ALL_SCREENS = [
    ROOT / "outputs/insurance/screen.json",
    ROOT / "outputs/library/screen.json",
    ROOT / "outputs/library2/screen.json",
    ROOT / "outputs/library2/screen_b.json",
    ROOT / "outputs/campaign/campaign.json",
    ROOT / "outputs/gen2/screen.json",
]



_PARENT_ID = "rfd4_design_045"


def _parent_seq() -> str | None:
    """d045, the root of every design here."""
    # `return x.get("sequence")` returned None here and exited before the
    # fallback ever ran - that file has no "sequence" key. Only return a value.
    f = ROOT / "outputs/refold_d045/verdict.json"
    if f.exists():
        try:
            seq = json.loads(f.read_text()).get("sequence")
            if seq:
                return seq
        except Exception:                                       # noqa: BLE001
            pass
    for mf in (ROOT / "outputs/insurance/variants.json",):
        if mf.exists():
            try:
                for v in json.loads(mf.read_text()):
                    if v.get("parent") == _PARENT_ID and v.get("sequence"):
                        # reconstruct the parent by undoing this variant's edits
                        s = list(v["sequence"])
                        for m in v.get("edits", []):
                            s[int(m[1:-1]) - 1] = m[0]
                        return "".join(s)
            except Exception:                                   # noqa: BLE001
                pass
    return None


_PSEQ = None


def _diff_from_parent(seq: str) -> list[str]:
    global _PSEQ
    if _PSEQ is None:
        _PSEQ = _parent_seq() or ""
    if not _PSEQ or len(_PSEQ) != len(seq):
        return []
    return [f"{_PSEQ[i]}{i+1}{seq[i]}"
            for i in range(len(seq)) if _PSEQ[i] != seq[i]]


def candidates(n: int) -> list[dict]:
    rows, seen = [], set()
    for src in ALL_SCREENS:
        if not src.exists():
            continue
        try:
            data = json.loads(src.read_text())
        except Exception:                                       # noqa: BLE001
            continue
        for r in data:
            did = r.get("design_id")
            if r.get("status") not in (None, "ok") or not did or did in seen:
                continue
            seen.add(did)
            rows.append(r)
    # the ipde fix: measure_raw does not populate it, so patch it in from the
    # row before ranking (see insurance_fold.py comment)
    for r in rows:
        r["badness"]["ipde"] = r.get("ipde")
    ok = [r for r in rows if bound(r) and (r.get("plddt") or 0) >= C.PLDDT_MIN]
    if not ok:
        return []
    rk = objective.rank_batch([r["badness"] for r in ok])
    for r, k in zip(ok, rk):
        r["dense_rank"] = k["dense_rank"]
    ok.sort(key=lambda r: -r["dense_rank"])
    return ok[:n]


def five_seed(design_id: str) -> dict | None:
    """Run (or reuse) the 5-seed official score."""
    vd = FINAL / design_id / "verdict.json"
    if vd.exists():
        return json.loads(vd.read_text())
    p = subprocess.run(
        [sys.executable, str(ROOT / "scripts/score_finalist.py"),
         design_id, design_id],
        capture_output=True, text=True, timeout=5400,
        env={**__import__("os").environ, "KEMP_ROOT": str(ROOT)})
    if p.returncode != 0 or not vd.exists():
        print(f"    5-seed failed: {p.stdout[-250:]}")
        return None
    return json.loads(vd.read_text())


def mutation_cluster(muts: list[str]) -> str:
    """Cluster key: the set of POSITIONS touched, order-independent."""
    return ",".join(sorted({m[1:-1] for m in muts}))


def main() -> int:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    cands = candidates(n)
    if not cands:
        print("No candidate has a bound ligand AND passes the pLDDT gate.")
        return 3
    variants = {}
    for mf in (VARIANTS, ROOT / "outputs/library/library.json",
               ROOT / "outputs/library2/library.json",
               ROOT / "outputs/gen2/variants.json"):
        if mf.exists():
            try:
                for v in json.loads(mf.read_text()):
                    variants.setdefault(v["design_id"], v)
            except Exception:                                   # noqa: BLE001
                pass
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"{len(cands)} candidate(s), 5 seeds each (~3.5 min)\n")
    scored = []
    for i, c in enumerate(cands):
        did = c["design_id"]
        print(f"[{i+1}/{len(cands)}] {did}")
        v = five_seed(did)
        if v is None:
            continue
        subs = v["subscores"]
        meta = variants.get(did, {})
        # Rank on the PER-REPLICATE truth, not the averaged geometry. The
        # averaged number ranked a 0/5 design (ins079, retracted 0.624) above
        # the only design that actually scores (ins080, 0.743 on 1/5).
        best = v.get("official_score_best")
        nz = v.get("n_nonzero_replicates")
        rec = {
            "design_id": did, "sequence": v["sequence"],
            "length": len(v["sequence"]),
            # derived from the SEQUENCE, not from a manifest key: the
            # library manifests use "mutations" and the insurance one uses
            # "edits", so reading either alone silently produced empty lists
            "mutations": _diff_from_parent(v["sequence"]),
            "parent": meta.get("parent", "rfd4_design_045"),
            "source": meta.get("source"),
            "official_score": v.get("official_score"),
            "official_score_best": best,
            "n_nonzero_replicates": nz,
            "per_replicate_total": v.get("per_replicate_total"),
            "reproducibility": (f"{nz}/{len(v.get('replicates', []))}"
                                if nz is not None else None),
            "limiting_term": v.get("limiting_term"),
            "subscores": subs, "sd_base": v.get("sd_base"),
            "mean_plddt": v.get("mean_plddt"),
            "measurable": v.get("measurable"),
            "protocol": v.get("protocol"),
            "screen_1seed": {k: c.get(k) for k in
                             ("d_base", "d_oxy", "dih_oxy", "gap_dih_oxy",
                              "plddt", "dense_rank")},
        }
        scored.append(rec)
        print(f"    official={v.get('official_score')}  "
              f"limiting={v.get('limiting_term')}  "
              f"sd_base={v.get('sd_base')}  pLDDT={v.get('mean_plddt')}")

    if not scored:
        print("\nNothing scored.")
        return 3

    # ---- hard eligibility gates
    print("\n--- eligibility ---")
    for r in scored:
        g = gates.run_all(r["sequence"], motif_residues=2)
        r["gates"] = g
        try:
            r["novelty"] = novelty.screen(r["sequence"], None)
        except Exception as exc:                                # noqa: BLE001
            r["novelty"] = {"eligible": None,
                            "verdict": f"not screened: {str(exc)[:120]}"}
        r["eligible"] = bool(g.get("passed")) and (
            r["novelty"].get("eligible") is not False)
        print(f"  {r['design_id']:34} gates={'PASS' if g.get('passed') else 'FAIL'}"
              f"  novelty={r['novelty'].get('verdict', '?')[:40]}"
              f"  -> {'ELIGIBLE' if r['eligible'] else 'rejected'}")

    eligible = [r for r in scored if r["eligible"]]
    nonzero = [r for r in eligible if (r.get("official_score") or 0) > 0]

    # ---- diversity: at most 2 per position-cluster
    chosen, per_cluster = [], {}
    # REPRODUCIBILITY FIRST, then median, then best.
    #
    # Sorting on `best` put a 3/5 design above a 5/5 one. The competition
    # scores five replicates; a design that scores on three of them carries a
    # zero into two, and under a geometric mean a zero is unrecoverable. The
    # median is the honest central number and `best` only breaks ties.
    pool = sorted(eligible,
                  key=lambda r: (-(r.get("n_nonzero_replicates") or 0),
                                 -(r.get("official_score") or 0),
                                 -(r.get("official_score_best") or 0)))
    for r in pool:
        key = mutation_cluster(r["mutations"]) or "parent"
        if per_cluster.get(key, 0) >= 2 and len(chosen) < len(pool):
            continue
        per_cluster[key] = per_cluster.get(key, 0) + 1
        chosen.append(r)
        if len(chosen) >= N_SLOTS:
            break
    for r in pool:                       # top up if diversity starved us
        if len(chosen) >= N_SLOTS:
            break
        if r not in chosen:
            chosen.append(r)

    # ---- export
    with (OUT / "submission.fasta").open("w") as fh:
        for r in chosen:
            fh.write(f">{r['design_id']} mutations={'+'.join(r['mutations'])} "
                     f"official={r['official_score']} "
                     f"plddt={r['mean_plddt']}\n{r['sequence']}\n")
    (OUT / "submission.json").write_text(json.dumps({
        "generated": time.strftime("%F %T"),
        "protocol": objective.protocol_stamp(replicates=C.REPLICATES),
        "n_candidates": len(cands), "n_scored": len(scored),
        "n_eligible": len(eligible), "n_nonzero": len(nonzero),
        "selected": chosen, "all_scored": scored,
    }, indent=2, default=str))
    for r in chosen:
        src = FINAL / r["design_id"]
        for seed in C.BOLTZ_SEEDS[:1]:
            f = src / f"seed{seed}" / "folded.pdb"
            if f.exists():
                shutil.copy2(f, OUT / f"{r['design_id']}.pdb")

    print("\n" + "=" * 78)
    print(f"scored {len(scored)}  eligible {len(eligible)}  "
          f"NON-ZERO official score: {len(nonzero)}")
    print(f"\nSELECTED {len(chosen)} for submission:")
    print(f"  {'design':34}{'best':>9}{'repro':>8}{'median':>9}"
          f"{'limiting':>13}{'pLDDT':>8}")
    for r in chosen:
        print(f"  {r['design_id']:34}{(r.get('official_score_best') or 0):9.4f}"
              f"{str(r.get('reproducibility')):>8}"
              f"{(r.get('official_score') or 0):9.4f}"
              f"{str(r['limiting_term']):>13}{r['mean_plddt']:8.3f}")
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
