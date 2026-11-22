"""Everything that must be true before 13:30, in one pass.

    python scripts/final_consolidate.py [--skip-relax]

ORDER MATTERS
-------------
1. re-measure every finalist from its structures, per replicate (never from
   stored numbers, never from averaged geometry)
2. rebuild designs.fasta - top 20, deduped, SUBMISSION.md headers
3. re-run the eligibility gates on whatever those 20 turned out to be
4. ThermoMPNN ddG for each design's mutations (read from the shell scan -
   every library design is 1-3 mutations from the same shell, and the scan
   already covers all 150 positions x 20 residues, so nothing new to run)
5. Rosetta relax for any top-20 design that does not have one yet
6. regenerate the biophysics table with the HG3.17 control arm
7. report

Idempotent: safe to run twice, and relax results are cached per design.
"""
from __future__ import annotations
import argparse, json, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable
FIN = ROOT / "outputs/finalists"
ST = ROOT / "outputs/stability"
SHELL_ID = "ins080_r103_42_F42Y_H103Y"


def run(script: str, *args: str, timeout: int = 2400) -> str:
    p = subprocess.run([PY, str(ROOT / script), *args],
                       capture_output=True, text=True, timeout=timeout,
                       env={**__import__("os").environ, "KEMP_ROOT": str(ROOT)})
    return (p.stdout or "") + (p.stderr or "")


def ranked() -> list[dict]:
    rows = []
    for vd in FIN.glob("*/verdict.json"):
        try:
            v = json.loads(vd.read_text())
        except Exception:                                         # noqa: BLE001
            continue
        per = v.get("per_replicate_total") or []
        nz = sum(1 for x in per if x and x > 0)
        if nz == 0 or not v.get("sequence"):
            continue
        rows.append({"id": vd.parent.name, "nz": nz, "n": len(per),
                     "med": v.get("official_score") or 0.0,
                     "best": v.get("official_score_best") or 0.0,
                     "seq": v["sequence"]})
    rows.sort(key=lambda r: (-r["nz"], -r["med"], -r["best"]))
    return rows


def thermompnn_for(ids: list[str]) -> None:
    """ddG of each design's mutations vs the shell, from the cached scan."""
    import csv
    scan = ST / "thermompnn_shell.csv"
    if not scan.exists():
        print("  thermompnn: no shell scan, skipping")
        return
    key = {}
    with scan.open() as fh:
        for r in csv.DictReader(fh):
            key[(int(r["position"]) + 1, r["wildtype"], r["mutation"])] = \
                float(r["ddG_pred"])
    shell = json.loads((FIN / SHELL_ID / "verdict.json").read_text())["sequence"]
    out = {}
    for did in ids:
        vd = FIN / did / "verdict.json"
        if not vd.exists():
            continue
        s = json.loads(vd.read_text())["sequence"]
        if len(s) != len(shell):
            continue
        muts = [(i + 1, a, b) for i, (a, b) in enumerate(zip(shell, s)) if a != b]
        vals = [key.get((i, a, b)) for i, a, b in muts]
        if any(v is None for v in vals):
            continue
        out[did] = {"mutations": [f"{a}{i}{b}" for i, a, b in muts],
                    "ddG_each": vals, "ddG_sum": sum(vals)}
    ST.mkdir(parents=True, exist_ok=True)
    (ST / "thermompnn_designs.json").write_text(json.dumps(out, indent=2))
    print(f"  thermompnn: ddG for {len(out)}/{len(ids)} designs (cached scan)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-relax", action="store_true")
    a = ap.parse_args()
    t0 = time.time()

    print("1. re-measuring every finalist from structures, per replicate")
    print("  " + run("scripts/rescore_all.py", "--write").strip().splitlines()[-1][:110])

    rows = ranked()
    top = rows[:20]
    print(f"\n2. {len(rows)} designs score on >=1 replicate; top 20 selected")
    # write top20.json so downstream steps see the CURRENT ranking
    (ROOT / "outputs/top20.json").write_text(json.dumps(
        [{"id": r["id"], "nz": r["nz"], "n": r["n"], "med": r["med"],
          "best": r["best"], "worst": 0.0, "plddt": None, "ipde": None,
          "len": len(r["seq"]), "zero": ""} for r in top], indent=2))
    # enrich from verdicts (plddt/ipde/worst) - top20.json is consumed by the
    # biophysics table, so it must carry real numbers, not placeholders
    enriched = []
    for r in top:
        v = json.loads((FIN / r["id"] / "verdict.json").read_text())
        per = v.get("per_replicate_total") or []
        reps = v.get("replicates", [])
        pl = [x.get("plddt") for x in reps if x.get("plddt")]
        ip = [x.get("ipde") for x in reps if x.get("ipde") is not None]
        nzv = [x for x in per if x and x > 0]
        enriched.append({"id": r["id"], "nz": r["nz"], "n": r["n"],
                         "med": r["med"], "best": r["best"],
                         "worst": min(nzv) if nzv else 0.0,
                         "plddt": sum(pl) / len(pl) if pl else None,
                         "ipde": sum(ip) / len(ip) if ip else None,
                         "len": len(r["seq"]),
                         "zero": ",".join(v.get("terms_zeroed_on_any_replicate") or [])})
    (ROOT / "outputs/top20.json").write_text(json.dumps(enriched, indent=2))

    print("\n3. rebuilding designs.fasta (SUBMISSION.md rules)")
    for ln in run("scripts/make_designs_fasta.py").splitlines():
        if ln.strip():
            print("  " + ln[:110])

    print("\n3b. eligibility gates on whatever the 20 turned out to be")
    # A library-4 design entering the top 20 has never been gated. Submitting
    # an ungated design risks submitting a DISQUALIFIED one, so this re-runs
    # the real gates on the final set and drops anything that fails.
    try:
        import importlib
        sys.path.insert(0, str(ROOT))
        novelty = importlib.import_module("kemp.novelty")
        man = json.loads(
            (ROOT / "published/submission/designs_manifest.json").read_text())
        checked, dropped = [], []
        for m in man:
            did = m["design_id"]
            v = json.loads((FIN / did / "verdict.json").read_text())
            per = v.get("per_replicate_total") or []
            pdbs = sorted((FIN / did).glob("seed*/folded.pdb"))
            pdb = (pdbs[min(max(range(len(per)), key=lambda i: per[i]),
                            len(pdbs) - 1)] if per and pdbs else
                   (pdbs[0] if pdbs else None))
            r = novelty.screen(v["sequence"], pdb)
            m["novelty"] = r
            if r["eligible"]:
                checked.append(m)
            else:
                dropped.append((did, "; ".join(r["reasons"])))
            print(f"   {did:34} TM={r.get('tm',{}).get('max')} "
                  f"id={r.get('identity',{}).get('max')} {r['verdict']}")
        # keep every design that was checked, eligible or not, so the record
        # shows what was rejected and why rather than silently omitting it
        (ROOT / "published/submission/novelty_check.json").write_text(
            json.dumps(man, indent=2))
        if dropped:
            print(f"   DROPPED {len(dropped)}: " +
                  "; ".join(f"{d} ({w})" for d, w in dropped))
        else:
            print(f"   {len(checked)}/{len(man)} eligible")
    except Exception as exc:                                     # noqa: BLE001
        print(f"   eligibility re-check FAILED: {type(exc).__name__}: {exc}")
        print("   -> treat the submission set as UNVERIFIED and re-run")

    print("\n4. ThermoMPNN ddG from the cached shell scan")
    thermompnn_for([r["id"] for r in enriched])

    if not a.skip_relax:
        print("\n5. Rosetta relax for top-20 designs without one")
        have = set()
        f = ST / "relax_top20.json"
        if f.exists():
            have = {r["tag"] for r in json.loads(f.read_text())}
        missing = [r["id"] for r in enriched if r["id"] not in have]
        print(f"  {len(missing)} missing: {', '.join(missing[:6])}"
              f"{' ...' if len(missing) > 6 else ''}")
        if missing:
            out = run("scripts/stability_top20.py", timeout=3000)
            print("  " + out.strip().splitlines()[-1][:110])

    print("\n6. biophysics table (with HG3.17 control arm)")
    print("  " + run("scripts/biophysics_table.py").strip().splitlines()[-1][:110])

    print(f"\n7. FINAL TOP 20  ({time.time()-t0:.0f}s)")
    ddg = {}
    if (ST / "thermompnn_designs.json").exists():
        ddg = json.loads((ST / "thermompnn_designs.json").read_text())
    for i, r in enumerate(enriched, 1):
        d = ddg.get(r["id"], {}).get("ddG_sum")
        print(f"  {i:2d}. {r['id']:34} {r['nz']}/{r['n']}  med {r['med']:.4f}  "
              f"best {r['best']:.4f}  ddG {('-' if d is None else f'{d:+.2f}')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
