"""Parallel Boltz-2 screening. Measured first, then used.

    python scripts/parallel_fold.py --bench            # measure the speedup
    python scripts/parallel_fold.py --manifest f.json --workers 6

WHY
---
Measured on this box while the serial campaign was running:

    13 folds in 542 s   = 41.7 s/fold
    GPU utilisation 5%, memory 1771 MiB of 24576 MiB

A 150 aa protein-ligand complex uses about 7% of this card's memory and
almost none of its compute. The wall time is model loading and CPU-side
preparation, not GPU math, so folds should overlap almost freely. Running
them one at a time wastes most of the hardware.

The worker count is NOT assumed. `--bench` folds the same small set serially
and in parallel and reports the real speedup, because a GPU that looks idle
can still serialise on kernel launches, and 96 GB of RAM can still thrash if
every worker loads its own copy of the model.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C, objective, official, pipeline      # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BASE_IDX, OXY_IDX = 102, 46
IN_WINDOW = (2.0, 3.6)


def _fold_one(task: tuple[str, str, str, int]) -> dict:
    """Runs in a separate process: fold one sequence and measure it."""
    design_id, seq, outdir, seed = task
    from kemp import official as _o, objective as _ob, pipeline as _p, config as _c
    from kemp.runners import boltz_predict as _bp
    site = _o.CatalyticSite(base_resseq=BASE_IDX, oxy_resseq=OXY_IDX, chain="A")
    w = Path(outdir)
    w.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    try:
        r = _bp(seq, w, seed=seed, ligand_ccd=_c.LIGAND_CCD)
        f = next((x for x in r.structures), None)
        if f is None:
            return {"design_id": design_id, "status": "fold_failed"}
        p = f if f.suffix == ".pdb" else _p.cif_to_pdb(f, w / "folded.pdb")
        raw = _o.measure_raw(str(p), site=site)
        conf = r.confidence or {}
        raw.ipde = conf.get("complex_ipde")
        plddt = conf.get("complex_plddt") or 0.0
        gap = _ob.dihedral_band_gap(raw.dihedral12_oxy)
        bnd = (raw.d_base is not None and raw.d_oxy is not None
               and IN_WINDOW[0] <= raw.d_base <= IN_WINDOW[1]
               and IN_WINDOW[0] <= raw.d_oxy <= IN_WINDOW[1])
        return {"design_id": design_id, "status": "ok", "sequence": seq,
                "pdb": str(p), "seed": seed,
                "d_base": raw.d_base, "d_oxy": raw.d_oxy,
                "dih_oxy": raw.dihedral12_oxy, "gap": gap,
                "plddt": plddt, "ipde": raw.ipde,
                "bound": bnd, "measurable": plddt >= _c.PLDDT_MIN,
                "subscores": _o.score(raw, strict=False)["subscores"],
                "badness": _ob.badness(raw),
                "wall_s": round(time.time() - t0, 1)}
    except Exception as exc:                                     # noqa: BLE001
        return {"design_id": design_id, "status": f"error: {type(exc).__name__}",
                "detail": str(exc)[:200], "wall_s": round(time.time() - t0, 1)}


def run(tasks: list[tuple], workers: int, out_json: Path | None = None,
        quiet: bool = False) -> list[dict]:
    rows: list[dict] = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_fold_one, t): t[0] for t in tasks}
        for i, fut in enumerate(as_completed(futs)):
            r = fut.result()
            rows.append(r)
            if not quiet and r.get("status") == "ok":
                flag = ""
                if r["bound"] and r["measurable"]:
                    flag = "  BOUND" + ("  <<< IN BAND" if r["gap"] == 0 else "")
                print(f"  [{i+1}/{len(tasks)}] {r['design_id']:26} "
                      f"db {r['d_base']:5.2f} do {r['d_oxy']:5.2f} "
                      f"dih {r['dih_oxy']:6.1f} gap {r['gap']:5.1f} "
                      f"pL {r['plddt']:.3f}{flag}", flush=True)
            elif not quiet:
                print(f"  [{i+1}/{len(tasks)}] {r['design_id']:26} "
                      f"{r.get('status')}", flush=True)
            if out_json:
                out_json.write_text(json.dumps(rows, indent=2, default=str))
    if not quiet:
        dt = time.time() - t0
        print(f"\n  {len(tasks)} folds in {dt/60:.1f} min "
              f"= {dt/max(1,len(tasks)):.1f} s/fold effective "
              f"({workers} workers)")
    return rows


def bench() -> int:
    man = json.loads((ROOT / "outputs/insurance/variants.json").read_text())
    seqs = [(v["design_id"], v["sequence"]) for v in man[:8]]
    base = ROOT / "outputs/bench_parallel"
    results = {}
    for w in (1, 4, 8):
        tasks = [(f"b{w}_{d}", s, str(base / f"w{w}" / d), 1)
                 for d, s in seqs[:w if w > 1 else 2]]
        n = len(tasks)
        t0 = time.time()
        run(tasks, w, quiet=True)
        dt = time.time() - t0
        results[w] = dt / n
        print(f"  workers={w:2}  {n} folds  {dt:6.1f} s  "
              f"=> {dt/n:5.1f} s/fold")
    base_rate = results.get(1, 42.0)
    print("\nspeedup vs serial:")
    for w, r in results.items():
        print(f"  workers={w:2}  {base_rate/r:.2f}x  "
              f"({3600/r:.0f} folds/hour)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", action="store_true")
    ap.add_argument("--manifest")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out", default="outputs/parallel/screen.json")
    ap.add_argument("--seed", type=int, default=C.BOLTZ_SEEDS[0])
    a = ap.parse_args()
    if a.bench:
        return bench()
    if not a.manifest:
        raise SystemExit("need --manifest or --bench")
    man = json.loads(Path(a.manifest).read_text())
    outdir = ROOT / "outputs/parallel/folds"
    tasks = [(v["design_id"], v["sequence"],
              str(outdir / v["design_id"]), a.seed) for v in man]
    out = ROOT / a.out
    out.parent.mkdir(parents=True, exist_ok=True)
    print(f"{len(tasks)} designs, {a.workers} workers\n")
    run(tasks, a.workers, out_json=out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
