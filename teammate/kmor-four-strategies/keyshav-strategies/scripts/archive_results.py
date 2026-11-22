"""Append-only archive of every design that ever scored.

    python scripts/archive_results.py

WHY THIS EXISTS
---------------
The top-20 table is a *view*: it changes as new libraries land, and a design
that was rank 8 this morning can fall off it this afternoon. That is fine for
the view and not fine for the record - the result was real when it was
measured and it stays real.

This file is the record. It only ever ADDS. A design already in it is updated
in place (its numbers can improve when re-measured) but never removed, and
the first time it was seen is kept.

ARCHIVE.md is sorted best-first for reading, but nothing is dropped from it,
including designs that score on a single replicate. run_monitor.md 4: a
partial pass is a real result reported with its reproducibility stated.
"""
from __future__ import annotations
import json
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIN = ROOT / "outputs/finalists"
ARCH = ROOT / "published/archive.json"


def main() -> int:
    old = {}
    if ARCH.exists():
        try:
            old = {r["design_id"]: r for r in json.loads(ARCH.read_text())}
        except Exception:                                        # noqa: BLE001
            old = {}
    now = datetime.now().isoformat(timespec="seconds")
    added = updated = 0

    for vd in sorted(FIN.glob("*/verdict.json")):
        try:
            v = json.loads(vd.read_text())
        except Exception:                                        # noqa: BLE001
            continue
        per = v.get("per_replicate_total") or []
        nz = sum(1 for x in per if x and x > 0)
        if nz == 0:
            continue
        did = vd.parent.name
        reps = v.get("replicates", [])
        pl = [x.get("plddt") for x in reps if x.get("plddt")]
        ip = [x.get("ipde") for x in reps if x.get("ipde") is not None]
        nzv = [x for x in per if x and x > 0]
        rec = {
            "design_id": did,
            "sequence": v.get("sequence"),
            "length": len(v["sequence"]) if v.get("sequence") else None,
            "reproducibility": f"{nz}/{len(per)}",
            "n_nonzero_replicates": nz,
            "official_median": v.get("official_score"),
            "official_best": v.get("official_score_best"),
            "official_worst": min(nzv) if nzv else None,
            "per_replicate_total": per,
            "terms_zeroed_on_any_replicate":
                v.get("terms_zeroed_on_any_replicate"),
            "plddt_mean": sum(pl) / len(pl) if pl else None,
            "ipde_mean": sum(ip) / len(ip) if ip else None,
            "last_measured": now,
        }
        if did in old:
            rec["first_seen"] = old[did].get("first_seen", now)
            if old[did].get("official_best") != rec["official_best"]:
                updated += 1
        else:
            rec["first_seen"] = now
            added += 1
        old[did] = rec

    rows = sorted(old.values(),
                  key=lambda r: (-(r["n_nonzero_replicates"] or 0),
                                 -(r["official_median"] or 0),
                                 -(r["official_best"] or 0)))
    ARCH.parent.mkdir(parents=True, exist_ok=True)
    ARCH.write_text(json.dumps(rows, indent=2))

    L = ["# Archive — every design that ever scored", "",
         f"_{len(rows)} designs. Append-only: nothing is ever removed from "
         f"this file. A design that drops out of the current top 20 stays "
         f"here with the numbers it was measured at._", "",
         "Scored **per replicate**, re-measured from the folded structures. "
         "`repro` is how many of the five Boltz-2 seeds produced a non-zero "
         "score; a design scoring on one seed is a real result reported with "
         "its reproducibility stated, not a hidden one.", "",
         "| design | repro | median | best | worst | pLDDT | iPDE | zeroed |",
         "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        def f(x, n=4):
            return "–" if x is None else f"{x:.{n}f}"
        L.append(f"| `{r['design_id']}` | {r['reproducibility']} | "
                 f"{f(r['official_median'])} | {f(r['official_best'])} | "
                 f"{f(r['official_worst'])} | {f(r['plddt_mean'],3)} | "
                 f"{f(r['ipde_mean'],2)} | "
                 f"{','.join(r.get('terms_zeroed_on_any_replicate') or []) or '–'} |")
    L.append("")
    (ROOT / "published/ARCHIVE.md").write_text("\n".join(L))
    print(f"archive: {len(rows)} designs ({added} new, {updated} re-measured)")
    print(f"  -> published/ARCHIVE.md, published/archive.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
