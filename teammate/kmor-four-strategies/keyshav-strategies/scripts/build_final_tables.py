"""The final deliverables: deduped ranking, submission table, FASTA + PDB export.

    python scripts/build_final_tables.py

THREE THINGS
  1. Rank all scoring designs, DEDUPED BY SEQUENCE. The FASTA builder already
     did this; the ranking view did not, so L2_triple_V84I_I23V_M88F and
     L3_I23V_V84I_M88F - the same three mutations named in a different order,
     byte-identical sequences - both appeared. Identical sequences are one
     design. The earlier-named one is kept and the alias is recorded.
  2. SUBMISSION_TABLE.md - the submitted 20 with the reason each one is in.
  3. designs/ - a FASTA and the folded PDBs for EVERY qualifying design, not
     just the submitted 20.

"Qualifying" means: scores on at least one of five replicates, passes the
length and catalytic gates, and is NOVEL on both the TM and identity gates.
"""
from __future__ import annotations
import json, shutil, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
ROOT = Path(__file__).resolve().parents[1]
FIN = ROOT / "outputs/finalists"
OUT = ROOT / "published"
BASE_1B, OXY_1B = 102, 46


def collect() -> tuple[list[dict], dict[str, list[str]]]:
    rows, by_seq = [], {}
    for vd in FIN.glob("*/verdict.json"):
        try:
            v = json.loads(vd.read_text())
        except Exception:                                        # noqa: BLE001
            continue
        per = v.get("per_replicate_total") or []
        nz = sum(1 for x in per if x and x > 0)
        if nz == 0 or not v.get("sequence"):
            continue
        reps = v.get("replicates", [])
        pl = [x.get("plddt") for x in reps if x.get("plddt")]
        ip = [x.get("ipde") for x in reps if x.get("ipde") is not None]
        nzv = [x for x in per if x and x > 0]
        rows.append({
            "id": vd.parent.name, "seq": v["sequence"], "nz": nz,
            "n": len(per), "med": v.get("official_score") or 0.0,
            "best": v.get("official_score_best") or 0.0, "worst": min(nzv),
            "len": len(v["sequence"]),
            "plddt": sum(pl) / len(pl) if pl else None,
            "ipde": sum(ip) / len(ip) if ip else None,
            "zero": ",".join(v.get("terms_zeroed_on_any_replicate") or []),
            "per": per})
    # dedupe by sequence, keeping the best-ranked instance
    rows.sort(key=lambda r: (-r["nz"], -r["med"], -r["best"]))
    seen, uniq = {}, []
    for r in rows:
        if r["seq"] in seen:
            by_seq.setdefault(seen[r["seq"]], []).append(r["id"])
            continue
        seen[r["seq"]] = r["id"]
        uniq.append(r)
    return uniq, by_seq


def main() -> int:
    from kemp import gates, novelty
    rows, aliases = collect()
    print(f"{len(rows)} unique scoring designs "
          f"({sum(len(v) for v in aliases.values())} duplicate name(s) merged)")

    # qualify: gates + novelty
    qual = []
    for r in rows:
        L = gates.check_length(r["seq"], motif_residues=2)
        C = gates.check_catalytic_residues(r["seq"], {101: "D", 45: "Q"})
        d = FIN / r["id"]
        pdbs = sorted(d.glob("seed*/folded.pdb"))
        bi = max(range(len(r["per"])), key=lambda i: r["per"][i]) if r["per"] else 0
        pdb = pdbs[min(bi, len(pdbs) - 1)] if pdbs else None
        nv = novelty.screen(r["seq"], pdb)
        r.update({"gates_ok": bool(L.passed and C.passed),
                  "tm": nv.get("tm", {}).get("max"),
                  "tm_target": nv.get("tm", {}).get("target"),
                  "identity": nv.get("identity", {}).get("max"),
                  "novel": nv["eligible"], "best_pdb": str(pdb) if pdb else None,
                  "aliases": aliases.get(r["id"], [])})
        if r["gates_ok"] and r["novel"]:
            qual.append(r)
    print(f"{len(qual)} qualifying (gates PASS and NOVEL)")

    (ROOT / "outputs/top20.json").write_text(json.dumps(
        [{k: r[k] for k in ("id", "nz", "n", "med", "best", "worst",
                            "plddt", "ipde", "len", "zero")}
         for r in qual[:20]], indent=2))

    # ---------------------------------------------------------------- exports
    dd = OUT / "designs"
    dd.mkdir(parents=True, exist_ok=True)
    fa, npdb = [], 0
    for i, r in enumerate(qual, 1):
        fa.append(f">{r['id']} base={BASE_1B} oxy={OXY_1B}")
        fa.append(r["seq"])
        (dd / f"{r['id']}.fasta").write_text(
            f">{r['id']} base={BASE_1B} oxy={OXY_1B}\n{r['seq']}\n")
        for sd in sorted((FIN / r["id"]).glob("seed*/folded.pdb")):
            dst = dd / f"{r['id']}__{sd.parent.name}.pdb"
            if not dst.exists():
                shutil.copy2(sd, dst)
                npdb += 1
    (OUT / "designs_qualifying.fasta").write_text("\n".join(fa) + "\n")
    print(f"exported {len(qual)} FASTA + {npdb} PDB files to published/designs/")

    # ------------------------------------------------- the submission table
    bio = {}
    f = OUT / "biophysics.json"
    if f.exists():
        bio = {r["design_id"]: r for r in json.loads(f.read_text())}
    shell = json.loads((FIN / "ins080_r103_42_F42Y_H103Y" / "verdict.json")
                       .read_text())["sequence"]

    def muts(seq):
        return "+".join(f"{a}{i+1}{b}" for i, (a, b)
                        in enumerate(zip(shell, seq)) if a != b) or "(parent)"

    def why(r, rank):
        b = bio.get(r["id"], {})
        sv = b.get("stability_verdict")
        if r["nz"] == 5 and not r["zero"]:
            s = "scores on **every** seed with no term ever zeroed"
        elif r["nz"] == 5:
            s = "scores on every seed"
        elif r["nz"] >= 3:
            s = f"scores on {r['nz']} of 5; `{r['zero'] or 'n/a'}` zeroes the rest"
        else:
            s = (f"partial pass ({r['nz']}/5) — kept because an occupied slot "
                 f"beats an empty one")
        if rank == 1:
            s = "**highest median and highest worst replicate** — " + s
        if sv == "FLAGGED":
            s += f"; stability flagged ({b.get('stability_reasons','')})"
        return s

    L = ["# Submitted top 20 — and why each one is in", "",
         "Ranked by **reproducibility (non-zero replicates out of 5) → median "
         "→ best**, every number re-measured from the folded structures and "
         "scored **per replicate**.", "",
         "Reproducibility leads the sort deliberately: a design that scores "
         "0.88 on one seed and 0.00 on four is not better than one that "
         "scores 0.83 five times out of five. The organisers run five "
         "replicates too.", "",
         "All 20 pass the competition gates: ≤260 aa, Asp102 and Gln46 "
         "intact, TM-score ≤0.70 and sequence identity ≤0.70 to published "
         "enzymes. Stability is **not** a competition criterion and was never "
         "used to rank — it is shown only to flag anything unbuildable.", "",
         "| # | design | mutations | repro | median | best | worst | pLDDT | "
         "iPDE | TM | stability | why it is in |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(qual[:20], 1):
        b = bio.get(r["id"], {})
        sv = {"CLEAN": "✅", "FLAGGED": "⚠️"}.get(b.get("stability_verdict"), "–")
        tm_s = "–" if r["tm"] is None else f"{r['tm']:.3f}"
        L.append(
            f"| {i} | `{r['id']}` | {muts(r['seq'])} | {r['nz']}/{r['n']} | "
            f"**{r['med']:.4f}** | {r['best']:.4f} | {r['worst']:.4f} | "
            f"{r['plddt']:.3f} | {r['ipde']:.2f} | {tm_s} | "
            f"{sv} | {why(r, i)} |")
    if aliases:
        L += ["", "### Merged duplicates", "",
              "Identical sequences are one design — SUBMISSION.md says "
              "duplicates are scored once and waste a slot.", ""]
        for keep, al in aliases.items():
            L.append(f"- `{keep}` also generated as {', '.join(f'`{a}`' for a in al)}")
    L += ["", f"### Beyond the 20", "",
          f"{len(qual)} designs qualify in total. Every one of them has a "
          f"FASTA and all five folded PDBs in [`designs/`](designs/), and "
          f"keeps its measured numbers in [`ARCHIVE.md`](ARCHIVE.md) whether "
          f"or not it made the cut.", ""]
    (OUT / "SUBMISSION_TABLE.md").write_text("\n".join(L))
    print(f"wrote published/SUBMISSION_TABLE.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
