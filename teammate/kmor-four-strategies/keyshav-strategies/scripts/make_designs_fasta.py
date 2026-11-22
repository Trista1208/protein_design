"""Build designs.fasta exactly as SUBMISSION.md requires.

THE RULES, taken from SUBMISSION.md rather than from memory:
  - at most 20 designs            (we were preparing 5 - this was leaving
                                   15 slots unused)
  - 50-260 residues
  - 20 standard uppercase amino acids
  - name: letters, digits, _ and -, unique, <= 64 chars
  - "base=" and "oxy=" on EVERY header, 1-BASED indices into the submitted
    sequence; base must be D or E, oxy must be one of Q N S T Y K R
  - NOTHING else in the header. "Free text after the name is rejected" - so
    the median=/best=/repro= annotations we had would have failed validation.
  - identical sequences are scored once and waste a slot, so dedupe.

Ranked by reproducibility (n/5 non-zero replicates) -> median -> best, scored
per replicate (never on averaged geometry).
"""
from __future__ import annotations
import json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIN = ROOT / "outputs/finalists"
OUT = ROOT / "published/submission"
MAX_DESIGNS = 20
BASE_1B, OXY_1B = 102, 46
BASE_OK, OXY_OK = set("DE"), set("QNSTYKR")
NAME_RE = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
AA = set("ACDEFGHIKLMNPQRSTVWY")


def collect() -> list[dict]:
    rows = []
    for vd in FIN.glob("*/verdict.json"):
        try:
            v = json.loads(vd.read_text())
        except Exception:                                         # noqa: BLE001
            continue
        per = v.get("per_replicate_total") or []
        seq = v.get("sequence")
        nz = sum(1 for x in per if x and x > 0)
        if not seq or nz == 0:
            continue
        rows.append({"id": vd.parent.name, "seq": seq, "nz": nz,
                     "n": len(per),
                     "med": v.get("official_score") or 0.0,
                     "best": v.get("official_score_best") or 0.0})
    rows.sort(key=lambda r: (-r["nz"], -r["med"], -r["best"]))
    return rows


def validate(name: str, seq: str) -> list[str]:
    bad = []
    if not NAME_RE.match(name):
        bad.append(f"name {name!r} breaks the charset/length rule")
    if not (50 <= len(seq) <= 260):
        bad.append(f"length {len(seq)} outside 50-260")
    if set(seq) - AA:
        bad.append(f"non-standard residues {sorted(set(seq)-AA)}")
    for idx, allowed, label in ((BASE_1B, BASE_OK, "base"),
                                (OXY_1B, OXY_OK, "oxy")):
        if idx > len(seq):
            bad.append(f"{label}={idx} beyond length {len(seq)}")
        elif seq[idx - 1] not in allowed:
            bad.append(f"{label}={idx} is {seq[idx-1]}, must be one of "
                       f"{''.join(sorted(allowed))}")
    if BASE_1B == OXY_1B:
        bad.append("base and oxy must be different residues")
    return bad


def main() -> int:
    rows = collect()
    print(f"{len(rows)} designs score on >=1 replicate")
    seen, chosen, dropped = {}, [], []
    for r in rows:
        if r["seq"] in seen:
            dropped.append((r["id"], seen[r["seq"]]))
            continue
        bad = validate(r["id"], r["seq"])
        if bad:
            print(f"  REJECT {r['id']}: {'; '.join(bad)}")
            continue
        seen[r["seq"]] = r["id"]
        chosen.append(r)
        if len(chosen) >= MAX_DESIGNS:
            break
    if dropped:
        print(f"\ndeduped {len(dropped)} identical sequence(s) "
              f"(would have wasted a slot):")
        for d, keep in dropped:
            print(f"    {d}  ==  {keep}")

    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    for r in chosen:
        lines.append(f">{r['id']} base={BASE_1B} oxy={OXY_1B}")
        lines.append(r["seq"])
    fa = OUT / "designs.fasta"
    fa.write_text("\n".join(lines) + "\n")

    meta = [{"rank": i, "design_id": r["id"], "reproducibility": f"{r['nz']}/{r['n']}",
             "official_median": r["med"], "official_best": r["best"],
             "length": len(r["seq"]), "base": BASE_1B, "oxy": OXY_1B}
            for i, r in enumerate(chosen, 1)]
    (OUT / "designs_manifest.json").write_text(json.dumps(meta, indent=2))

    # A second, LONGER file for the team: every design that scores, up to 30,
    # same header format. This is NOT the submission (the limit is 20) - it
    # exists so nothing good is hidden just because it fell below the cut.
    # Designs past rank 20 keep their numbers in ARCHIVE.md regardless.
    ext_seen, ext = set(), []
    for r in rows:
        if r["seq"] in ext_seen or validate(r["id"], r["seq"]):
            continue
        ext_seen.add(r["seq"])
        ext.append(r)
        if len(ext) >= 30:
            break
    ex_lines = []
    for i, r in enumerate(ext, 1):
        ex_lines.append(f">{r['id']} base={BASE_1B} oxy={OXY_1B}")
        ex_lines.append(r["seq"])
    (OUT / "designs_top30.fasta").write_text("\n".join(ex_lines) + "\n")
    (OUT / "designs_top30.json").write_text(json.dumps(
        [{"rank": i, "design_id": r["id"],
          "reproducibility": f"{r['nz']}/{r['n']}",
          "official_median": r["med"], "official_best": r["best"],
          "submitted": i <= MAX_DESIGNS}
         for i, r in enumerate(ext, 1)], indent=2))
    print(f"wrote {OUT/'designs_top30.fasta'}  ({len(ext)} designs, "
          f"top {MAX_DESIGNS} of them are the submission)")

    print(f"\nwrote {fa}  ({len(chosen)} designs, {MAX_DESIGNS} allowed)")
    for m in meta:
        print(f"  {m['rank']:2d}. {m['design_id']:34} {m['reproducibility']}  "
              f"median {m['official_median']:.4f}  best {m['official_best']:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
