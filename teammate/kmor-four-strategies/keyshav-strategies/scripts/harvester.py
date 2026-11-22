"""Autonomous finalist harvester: 5-seed whatever the campaign turns up.

    python scripts/harvester.py --until 11:15

Runs beside the campaign. Every cycle it re-reads the campaign output, picks
designs that deserve the 5-seed cost, scores them PER REPLICATE, and re-banks
the submission set. No manual step anywhere.

WHO EARNS 5 SEEDS, in priority order
------------------------------------
  1. bound ligand AND dihedral in band on the screening seed - the only
     designs that can possibly score, since `orient_oxy` is the term zeroing
     everything
  2. bound ligand AND gap <= 20 deg - near misses worth confirming
  3. best dense rank among the rest

A 1-seed screen cannot compute `cons` and cannot show reproducibility, so a
screening hit is a candidate, never a result (docs/good_practices.md §2).

SCORING IS PER REPLICATE, ALWAYS
--------------------------------
`official.score_replicates()`, re-measured from the folded structures. The
averaged-geometry route invented 0.624 for a design where 0/5 replicates
scored and erased a real 0.743 on another, both inside one hour. Results are
reported as "best achievable, n/5 reproducibility".

The harvester shares one GPU with the campaign. Boltz on a 150 aa complex
uses roughly 4 GB of 24 GB, so the two coexist; the cycle sleep keeps the
contention mild rather than eliminating it.
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
from kemp import config as C, objective, official                # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / "outputs/campaign/campaign.json"
GEN2 = ROOT / "outputs/gen2/screen.json"
FIN = ROOT / "outputs/finalists"
STATE = ROOT / "outputs/harvester_state.json"
SITE = official.CatalyticSite(base_resseq=102, oxy_resseq=46, chain="A")
IN_WINDOW = (2.0, 3.6)


def load_rows() -> list[dict]:
    rows = []
    for src in (CAMPAIGN, GEN2):
        if src.exists():
            try:
                rows += json.loads(src.read_text())
            except Exception:                                   # noqa: BLE001
                pass
    return rows


def deserving(rows: list[dict], done: set[str], limit: int) -> list[dict]:
    ok = [r for r in rows
          if r.get("design_id") not in done
          and r.get("bound") and r.get("measurable")
          and r.get("sequence")]
    if not ok:
        return []
    in_band = [r for r in ok if (r.get("gap") if r.get("gap") is not None else 99) == 0.0]
    near = [r for r in ok if 0 < (r.get("gap") if r.get("gap") is not None else 99) <= 20.0]
    rest = [r for r in ok if (r.get("gap") if r.get("gap") is not None else 99) > 20.0]
    for group in (in_band, near, rest):
        group.sort(key=lambda r: (r.get("gap") or 99))
    return (in_band + near + rest)[:limit]


def five_seed(design_id: str, sequence: str) -> dict | None:
    """Fold 5 seeds and score PER REPLICATE, re-measured from structures."""
    d = FIN / design_id
    d.mkdir(parents=True, exist_ok=True)
    (d / "sequence.txt").write_text(sequence + "\n")
    p = subprocess.run(
        [sys.executable, str(ROOT / "scripts/score_finalist.py"),
         sequence, design_id],
        capture_output=True, text=True, timeout=5400,
        env={**__import__("os").environ, "KEMP_ROOT": str(ROOT)})
    vd = d / "verdict.json"
    if not vd.exists():
        print(f"    {design_id}: 5-seed produced nothing "
              f"({p.stdout[-160:].strip()})")
        return None
    pdbs = sorted(d.glob("seed*/folded.pdb"))
    if not pdbs:
        return None
    v = json.loads(vd.read_text())
    raws = [official.measure_raw(str(x), site=SITE) for x in pdbs]
    ipdes = [r.get("ipde") for r in v.get("replicates", [])
             if r.get("ipde") is not None]
    rep = official.score_replicates(raws, ipde_values=ipdes)
    v.update({"official_score": rep["median_total"],
              "official_score_best": rep["best_total"],
              "official_score_worst": rep["worst_total"],
              "n_nonzero_replicates": rep["n_nonzero_replicates"],
              "per_replicate_total": rep["per_replicate_total"],
              "terms_zeroed_on_any_replicate": rep["terms_zeroed_on_any_replicate"],
              "averaged_geometry_total_DIAGNOSTIC": rep["averaged_geometry_total"]})
    vd.write_text(json.dumps(v, indent=2, default=str))
    return rep


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--until", default="11:15")
    ap.add_argument("--cycle", type=int, default=420)
    ap.add_argument("--per-cycle", type=int, default=4)
    a = ap.parse_args()
    hh, mm = (int(x) for x in a.until.split(":"))
    stop = datetime.now().replace(hour=hh, minute=mm, second=0).timestamp()

    done = set()
    if STATE.exists():
        try:
            done = set(json.loads(STATE.read_text()).get("done", []))
        except Exception:                                       # noqa: BLE001
            pass
    done |= {p.parent.name for p in FIN.glob("*/verdict.json")}

    hits: list[dict] = []
    while time.time() < stop - 240:
        rows = load_rows()
        picks = deserving(rows, done, a.per_cycle)
        if not picks:
            time.sleep(a.cycle)
            continue
        print(f"[{datetime.now():%H:%M}] {len(rows)} screened, "
              f"5-seeding {len(picks)}", flush=True)
        for r in picks:
            did = r["design_id"]
            done.add(did)
            try:
                rep = five_seed(did, r["sequence"])
            except Exception as exc:                            # noqa: BLE001
                print(f"    {did}: {type(exc).__name__}: {str(exc)[:110]}")
                continue
            if rep is None:
                continue
            msg = (f"    {did:26} best {rep['best_total']:.4f}  "
                   f"repro {rep['n_nonzero_replicates']}/{rep['n_replicates']}"
                   f"  median {rep['median_total']:.4f}  "
                   f"zeroed {','.join(rep['terms_zeroed_on_any_replicate'])}")
            if rep["best_total"] > 0:
                msg += "   <<< SCORES"
                hits.append({"design_id": did, **rep,
                             "mutations": r.get("mutations")})
            print(msg, flush=True)
            STATE.write_text(json.dumps(
                {"done": sorted(done), "hits": hits}, indent=2, default=str))

        # re-bank after every cycle so the submission set is always current
        subprocess.run(
            [sys.executable, str(ROOT / "scripts/bank_submissions.py"), "14"],
            capture_output=True, text=True, timeout=7200,
            env={**__import__("os").environ, "KEMP_ROOT": str(ROOT)})
        time.sleep(a.cycle)

    print(f"\nharvester done. {len(hits)} design(s) scored on >=1 replicate:")
    for h in sorted(hits, key=lambda x: -x["best_total"]):
        print(f"   {h['design_id']:26} best {h['best_total']:.4f}  "
              f"repro {h['n_nonzero_replicates']}/{h['n_replicates']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
