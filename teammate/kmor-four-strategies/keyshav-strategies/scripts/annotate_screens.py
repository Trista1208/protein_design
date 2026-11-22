"""Back-fill literature filters and the Gln46 pinner onto screened designs.

    python scripts/annotate_screens.py [--loop 180 --until 11:55]

WHY A SEPARATE PASS
-------------------
`parallel_fold.py` was launched at 06:05 and never called `litfilters`; the
pinner metric was not written until 06:24. Restarting the screen to pick the
new code up would discard ~70 completed folds for a measurement that can be
recovered from the structures already on disk.

So this reads each screen's JSON, recomputes the literature filters and
`gln_pinner` from the saved PDB, and writes them back. Running it on a loop
keeps the live table current while the screen continues underneath.

The pinner is the generalised form of the Tyr103 observation and the reason
this matters: measured across 50 bound+measurable designs, the nearest
H-bond donor to the Gln46 amide sits at a median 2.34 A in the 10 designs
whose dihedral is in band and 3.34 A in the 40 that are not. It is a 1-seed
proxy for something that otherwise costs five folds to establish - but it is
necessary, not sufficient, so it ranks and explains rather than gates.
"""
from __future__ import annotations

import argparse
import json
import time
from datetime import datetime
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import litfilters                                      # noqa: E402
from kemp.theozyme import load                                   # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BASE_IDX, OXY_IDX = 102, 46
SCREENS = [
    (ROOT / "outputs/library4/screen.json", ROOT / "outputs/parallel/folds"),
    (ROOT / "outputs/library/screen.json", ROOT / "outputs/parallel/folds"),
    (ROOT / "outputs/library2/screen.json", ROOT / "outputs/parallel/folds"),
    (ROOT / "outputs/campaign/campaign.json", ROOT / "outputs/campaign"),
]


def find_pdb(design_id: str, base: Path, rec: dict) -> Path | None:
    if rec.get("pdb") and Path(rec["pdb"]).exists():
        return Path(rec["pdb"])
    direct = base / design_id / "folded.pdb"
    if direct.exists():
        return direct
    hits = list(base.rglob(f"{design_id}/folded.pdb"))
    return hits[0] if hits else None


def annotate_one(path: Path, base: Path) -> tuple[int, int]:
    if not path.exists():
        return 0, 0
    try:
        rows = json.loads(path.read_text())
    except Exception:                                            # noqa: BLE001
        return 0, 0
    added = 0
    for r in rows:
        if r.get("status") not in (None, "ok"):
            continue
        if r.get("lit") is not None and r.get("pinner") is not None:
            continue
        # only designs whose geometry means anything (good_practices 2)
        if not (r.get("bound") and r.get("measurable")):
            continue
        p = find_pdb(r.get("design_id", ""), base, r)
        if p is None:
            continue
        try:
            tz = load(str(p), require_ligand=False)
            lf = litfilters.evaluate(tz, base_resseq=BASE_IDX,
                                     oxy_resseq=OXY_IDX, chain="A")
            r["lit"] = {"values": {k: v for k, v in lf.values.items()
                                   if not isinstance(v, dict)},
                        "hard_fail": lf.hard_fail, "soft_pass": lf.soft_pass}
            r["pinner"] = lf.values.get("pinner")
            sasa = litfilters.base_relative_sasa(p, BASE_IDX, "A")
            r["lit"]["values"]["H2_base_rel_sasa"] = sasa.get("rel_sasa_pct")
            added += 1
        except Exception as exc:                                 # noqa: BLE001
            r["lit"] = {"error": f"{type(exc).__name__}: {str(exc)[:110]}"}
    if added:
        path.write_text(json.dumps(rows, indent=2, default=str))
    return added, len(rows)


def report(paths) -> None:
    allrows = []
    for path, _ in paths:
        if path.exists():
            try:
                allrows += json.loads(path.read_text())
            except Exception:                                    # noqa: BLE001
                pass
    pinned = [r for r in allrows if isinstance(r.get("pinner"), dict)
              and r["pinner"].get("distance_A") is not None]
    if not pinned:
        return
    inb = [r for r in pinned if (r.get("gap") if r.get("gap") is not None else 99) == 0.0]
    outb = [r for r in pinned if (r.get("gap") if r.get("gap") is not None else 99) != 0.0]
    import statistics as st
    f = lambda g: (st.median(r["pinner"]["distance_A"] for r in g)  # noqa: E731
                   if g else float("nan"))
    print(f"   pinner distance: in-band n={len(inb)} median {f(inb):.2f} A | "
          f"out-of-band n={len(outb)} median {f(outb):.2f} A", flush=True)
    who = {}
    for r in inb:
        who[r["pinner"]["residue"]] = who.get(r["pinner"]["residue"], 0) + 1
    if who:
        print(f"   in-band pinners: {who}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--loop", type=int, default=0)
    ap.add_argument("--until", default=None)
    a = ap.parse_args()
    stop = None
    if a.until:
        hh, mm = (int(x) for x in a.until.split(":"))
        stop = datetime.now().replace(hour=hh, minute=mm, second=0).timestamp()
    while True:
        tot = 0
        for path, base in SCREENS:
            n, _ = annotate_one(path, base)
            tot += n
        print(f"[{datetime.now():%H:%M:%S}] annotated {tot} design(s)",
              flush=True)
        report(SCREENS)
        if not a.loop or (stop and time.time() > stop):
            return 0
        time.sleep(a.loop)


if __name__ == "__main__":
    raise SystemExit(main())
