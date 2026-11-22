"""Rerank and republish everything, every few minutes, until the deadline.

    python scripts/live_loop.py [--until 13:25] [--every 180]

EACH CYCLE
  1. re-measure every finalist from its structures, per replicate
  2. rerank - the top-20 view changes as new designs land
  3. rebuild designs.fasta (the 20 submitted) and designs_top30.fasta
  4. refresh the append-only archive of everything that ever scored
  5. regenerate the biophysics/thermal/biochemistry tables
  6. rewrite the results README with all of it inline
  7. copy PDBs for anything in the top 30, then commit and push

NOTHING IS EVER DELETED. Structures, GIFs and archive rows only accumulate.
A design that falls out of the top 20 keeps its PDB, its GIF and its row in
ARCHIVE.md with the numbers it was measured at - the ranking is a view, the
archive is the record.

Rosetta relax and the MD are NOT re-run here; they take minutes per structure
and their results are cached by design id. A design that has no stability
number yet shows "-" rather than being dropped or imputed.
"""
from __future__ import annotations
import argparse, json, shutil, subprocess, sys, time
from datetime import datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable
TEAM = Path("/home/mandalorian/Hackathons/team-sharing")
DEST = TEAM / "keyshav-strategies/results"
FIN = ROOT / "outputs/finalists"
LOG = ROOT / "logs/live_loop.log"


def say(m: str) -> None:
    line = f"[{datetime.now():%H:%M:%S}] {m}"
    print(line, flush=True)
    with LOG.open("a") as f:
        f.write(line + "\n")


def sh(script: str, *a: str, timeout: int = 1200) -> str:
    try:
        p = subprocess.run([PY, str(ROOT / script), *a], capture_output=True,
                           text=True, timeout=timeout,
                           env={**__import__("os").environ,
                                "KEMP_ROOT": str(ROOT)})
        return p.stdout or ""
    except Exception as exc:                                     # noqa: BLE001
        say(f"  {script}: {type(exc).__name__}")
        return ""


def ranked() -> list[dict]:
    rows = []
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
        rows.append({"id": vd.parent.name, "nz": nz, "n": len(per),
                     "med": v.get("official_score") or 0.0,
                     "best": v.get("official_score_best") or 0.0,
                     "worst": min(nzv), "len": len(v["sequence"]),
                     "plddt": sum(pl) / len(pl) if pl else None,
                     "ipde": sum(ip) / len(ip) if ip else None,
                     "zero": ",".join(v.get("terms_zeroed_on_any_replicate") or [])})
    rows.sort(key=lambda r: (-r["nz"], -r["med"], -r["best"]))
    return rows


def copy_structures(top: list[dict]) -> int:
    """Add PDBs for the top 30. Never removes one already published."""
    out = DEST / "structures"
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    for r in top[:30]:
        d = FIN / r["id"]
        for sd in sorted(d.glob("seed*/folded.pdb")):
            dst = out / f"{r['id']}__{sd.parent.name}.pdb"
            if not dst.exists():
                shutil.copy2(sd, dst)
                n += 1
    return n


def rewrite_readme() -> bool:
    import importlib.util as u
    spec = u.spec_from_file_location("lp", ROOT / "scripts/live_publish_team.py")
    m = u.module_from_spec(spec)
    sys.path.insert(0, str(ROOT))
    spec.loader.exec_module(m)
    new = "\n".join(m.biophysics_section())
    p = DEST / "README.md"
    if not p.exists() or not new.strip():
        return False
    t = p.read_text()
    a = t.find("## Stability — biochemistry")
    b = t.find("## Correction — if you saw an earlier commit")
    if a == -1 or b == -1 or b <= a:
        return False
    p.write_text(t[:a] + new + "\n" + t[b:])
    return True


def header_table(rows: list[dict]) -> None:
    """Keep the summary box at the top of the README honest."""
    p = DEST / "README.md"
    if not p.exists():
        return
    t = p.read_text().splitlines()
    best = max((r["best"] for r in rows), default=0.0)
    bn = next((r for r in rows if r["best"] == best), None)
    five = sum(1 for r in rows if r["nz"] == 5)
    for i, ln in enumerate(t[:30]):
        if ln.startswith("| **Best official score**"):
            t[i] = (f"| **Best official score** | **{best:.4f}** "
                    f"({bn['nz']}/5 seeds, `{bn['id']}`) |" if bn else ln)
        elif ln.startswith("| Scoring on ≥1 replicate"):
            t[i] = f"| Scoring on ≥1 replicate | **{len(rows)}** |"
        elif ln.startswith("| Scoring on 5/5"):
            t[i] = f"| Scoring on 5/5 seeds | **{five}** |"
        elif ln.startswith("_Auto-updated"):
            t[i] = f"_Auto-updated live. Last: {datetime.now():%Y-%m-%d %H:%M:%S}_"
    p.write_text("\n".join(t) + "\n")


def push(msg: str) -> bool:
    subprocess.run(["git", "add", "-A", "keyshav-strategies/results"],
                   cwd=TEAM, capture_output=True, text=True, timeout=300)
    c = subprocess.run(["git", "-c", "user.email=keyshav.mor@lemna.bio",
                        "-c", "user.name=Keyshav Mor", "commit", "-m", msg],
                       cwd=TEAM, capture_output=True, text=True, timeout=300)
    if "nothing to commit" in (c.stdout + c.stderr):
        return False
    p = subprocess.run(["git", "push", "origin", "kmor/four-strategies"],
                       cwd=TEAM, capture_output=True, text=True, timeout=600)
    return p.returncode == 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--until", default="13:25")
    ap.add_argument("--every", type=int, default=180)
    a = ap.parse_args()
    h, m = (int(x) for x in a.until.split(":"))
    stop = datetime.now().replace(hour=h, minute=m, second=0).timestamp()

    cycle = 0
    while time.time() < stop:
        cycle += 1
        t0 = time.time()
        say(f"cycle {cycle}")
        sh("scripts/rescore_all.py", "--write", timeout=900)
        rows = ranked()
        top = rows[:20]
        (ROOT / "outputs/top20.json").write_text(json.dumps(
            [{"id": r["id"], "nz": r["nz"], "n": r["n"], "med": r["med"],
              "best": r["best"], "worst": r["worst"], "plddt": r["plddt"],
              "ipde": r["ipde"], "len": r["len"], "zero": r["zero"]}
             for r in top], indent=2))
        sh("scripts/make_designs_fasta.py", timeout=300)
        sh("scripts/archive_results.py", timeout=300)
        sh("scripts/biophysics_table.py", timeout=600)
        sh("scripts/rankings.py", timeout=300)

        for name in ("BIOPHYSICS.md", "biophysics.csv", "biophysics.json",
                     "ARCHIVE.md", "archive.json", "top20.json",
                     "RANKINGS.md", "rankings.json"):
            src = ROOT / "published" / name
            if not src.exists():
                src = ROOT / "outputs" / name
            if src.exists():
                shutil.copy2(src, DEST / name)
        (DEST / "submission").mkdir(parents=True, exist_ok=True)
        for name in ("designs.fasta", "designs_manifest.json",
                     "designs_top30.fasta", "designs_top30.json",
                     "novelty_check.json"):
            src = ROOT / "published/submission" / name
            if src.exists():
                shutil.copy2(src, DEST / "submission" / name)
        # also at the results/ ROOT - the FASTAs are what people actually come
        # here to download, and burying them one level down hides them
        for name in ("designs.fasta", "designs_top30.fasta"):
            src = ROOT / "published/submission" / name
            if src.exists():
                shutil.copy2(src, DEST / name)
        for name in ("thermal_md.json", "relax_top20.json",
                     "thermompnn_designs.json"):
            src = ROOT / "outputs/stability" / name
            if src.exists():
                shutil.copy2(src, DEST / name)

        nstruct = copy_structures(rows)
        ok = rewrite_readme()
        header_table(rows)
        best = max((r["best"] for r in rows), default=0.0)
        five = sum(1 for r in rows if r["nz"] == 5)
        msg = (f"live: {len(rows)} scoring, {five} at 5/5, best "
               f"{best:.4f} ({datetime.now():%H:%M:%S})")
        pushed = push(msg)
        say(f"  {len(rows)} scoring ({five} at 5/5), best {best:.4f}; "
            f"+{nstruct} PDBs; readme={'ok' if ok else 'SKIPPED'}; "
            f"push={'yes' if pushed else 'no change'} ({time.time()-t0:.0f}s)")

        left = stop - time.time()
        if left <= 0:
            break
        time.sleep(min(a.every, max(5, left)))
    say("live loop finished")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
