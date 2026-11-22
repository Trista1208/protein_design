"""5-seed every in-band design from every screen, best first.

    python scripts/score_inband.py [--until 11:20] [--max 40]

WHY THIS AND NOT THE SUPERVISOR
-------------------------------
The supervisor crashed in phase 4 (`NameError: _pin`), so no 5-seed scoring
ran while library 2 was screening. Restarting it would redo phase 3 and
re-fold all 140 library-2 designs. This skips straight to the work that has
value: the designs whose oxyanion dihedral is already in band.

That set is now substantial. Library 1 produced 15 in-band designs out of 146
bound+measurable (10%); library 2's PAIRS are producing them at roughly 47%,
which is the first large effect anything has shown.

SCORED PER REPLICATE, re-measured from the structures. A 1-seed screen cannot
compute `cons` and cannot show reproducibility, so an in-band screening hit is
a candidate, never a result (docs/good_practices.md 5c).
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C, official                          # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIN = ROOT / "outputs/finalists"
SITE = official.CatalyticSite(base_resseq=102, oxy_resseq=46, chain="A")
SCREENS = [ROOT / "outputs/library4/screen.json",
           ROOT / "outputs/library3/screen.json",
           ROOT / "outputs/library2/screen_b.json",
           ROOT / "outputs/library2/screen.json",
           ROOT / "outputs/library/screen.json",
           ROOT / "outputs/campaign/campaign.json"]


def candidates(limit: int) -> list[dict]:
    seen, out = set(), []
    for src in SCREENS:
        if not src.exists():
            continue
        try:
            rows = json.loads(src.read_text())
        except Exception:                                        # noqa: BLE001
            continue
        for r in rows:
            did = r.get("design_id")
            if (r.get("status") not in (None, "ok") or not did or did in seen
                    or not r.get("bound") or not r.get("measurable")
                    or not r.get("sequence")):
                continue
            gap = r.get("gap")
            if gap is None:
                continue
            seen.add(did)
            out.append({**r, "_gap": gap})
    # in band first, then by how close, then by how good the distances are
    tb, to = official.TARGETS["d_base"], official.TARGETS["d_oxy"]

    def quality(r):
        return -(official._window_score(r["d_base"], tb["target"], *tb["window"])
                 * official._window_score(r["d_oxy"], to["target"], *to["window"]))
    out.sort(key=lambda r: (r["_gap"] != 0.0, r["_gap"], quality(r)))
    done = {p.parent.name for p in FIN.glob("*/verdict.json")}
    return [r for r in out if r["design_id"] not in done][:limit]


def five_seed(did: str, seq: str) -> dict | None:
    p = subprocess.run(
        [sys.executable, str(ROOT / "scripts/score_finalist.py"), seq, did],
        capture_output=True, text=True, timeout=5400,
        env={**__import__("os").environ, "KEMP_ROOT": str(ROOT)})
    d = FIN / did
    pdbs = sorted(d.glob("seed*/folded.pdb"))
    vd = d / "verdict.json"
    if not pdbs or not vd.exists():
        return None
    v = json.loads(vd.read_text())
    raws = [official.measure_raw(str(x), site=SITE) for x in pdbs]
    ip = [r.get("ipde") for r in v.get("replicates", [])
          if r.get("ipde") is not None]
    rep = official.score_replicates(raws, ipde_values=ip)
    v.update({"official_score": rep["median_total"],
              "official_score_best": rep["best_total"],
              "official_score_worst": rep["worst_total"],
              "n_nonzero_replicates": rep["n_nonzero_replicates"],
              "per_replicate_total": rep["per_replicate_total"],
              "terms_zeroed_on_any_replicate": rep["terms_zeroed_on_any_replicate"]})
    vd.write_text(json.dumps(v, indent=2, default=str))
    return rep


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--until", default="11:20")
    ap.add_argument("--max", type=int, default=40)
    a = ap.parse_args()
    hh, mm = (int(x) for x in a.until.split(":"))
    stop = datetime.now().replace(hour=hh, minute=mm, second=0).timestamp()

    hits = []
    while time.time() < stop - 300:
        cands = candidates(a.max)
        if not cands:
            print(f"[{datetime.now():%H:%M}] nothing new to score; waiting",
                  flush=True)
            time.sleep(240)
            continue
        inband = sum(1 for c in cands if c["_gap"] == 0.0)
        print(f"[{datetime.now():%H:%M}] {len(cands)} candidates "
              f"({inband} in band)", flush=True)
        for c in cands:
            if time.time() > stop - 240:
                break
            did = c["design_id"]
            try:
                rep = five_seed(did, c["sequence"])
            except Exception as exc:                             # noqa: BLE001
                print(f"   {did}: {type(exc).__name__}", flush=True)
                continue
            if rep is None:
                continue
            tag = ""
            if rep["best_total"] > 0:
                tag = "   <<< SCORES"
                hits.append({"design_id": did, "best": rep["best_total"],
                             "repro": rep["n_nonzero_replicates"],
                             "mutations": c.get("mutations")})
            print(f"   {did[:38]:38} gap {c['_gap']:5.1f}  "
                  f"best {rep['best_total']:.4f}  "
                  f"{rep['n_nonzero_replicates']}/{rep['n_replicates']}"
                  f"  zeroed {','.join(rep['terms_zeroed_on_any_replicate'])}"
                  f"{tag}", flush=True)
            subprocess.run(
                [sys.executable, str(ROOT / "scripts/bank_submissions.py"), "25"],
                capture_output=True, text=True, timeout=7200,
                env={**__import__("os").environ, "KEMP_ROOT": str(ROOT)})

    print(f"\n{len(hits)} design(s) score on >=1 replicate:")
    for h in sorted(hits, key=lambda x: -x["best"]):
        print(f"   {h['design_id']:38} {h['best']:.4f}  {h['repro']}/5  "
              f"{'+'.join(h.get('mutations') or [])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
