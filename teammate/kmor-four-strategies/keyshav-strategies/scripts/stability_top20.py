"""Fold-stability metrics for the top-20 finalists, with a control arm.

Rosetta relax (protein only) on the best-scoring seed of each design, plus the
HG3.17 crystal structure as the control. A stability number with no control is
not evidence (good_practices 5b): HG3.17 has been expressed, purified,
crystallised and kinetically characterised, so it calibrates what the scale
means.

PROTEIN ONLY - the question is fold stability, and scoring a ligand without its
params file would silently drop it anyway.
"""
from __future__ import annotations
import json, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/stability"
RELAX = Path.home() / "opt/rosetta/source/bin/relax.default.linuxgccrelease"
DB = Path.home() / "opt/rosetta/database"


def best_seed_pdb(design_id: str) -> Path | None:
    d = ROOT / "outputs/finalists" / design_id
    vd = d / "verdict.json"
    pdbs = sorted(d.glob("seed*/folded.pdb"))
    if not pdbs:
        return None
    if vd.exists():
        per = json.loads(vd.read_text()).get("per_replicate_total") or []
        if per:
            i = max(range(len(per)), key=lambda k: per[k])
            return pdbs[min(i, len(pdbs) - 1)]
    return pdbs[0]


def strip_ligand(src: Path, dst: Path) -> None:
    with dst.open("w") as f:
        for ln in src.read_text().splitlines():
            if ln.startswith(("ATOM", "TER")):
                f.write(ln + "\n")
        f.write("END\n")


def relax(pdb: Path, tag: str) -> dict | None:
    d = OUT / tag
    d.mkdir(parents=True, exist_ok=True)
    sc = d / "score.sc"
    stripped = d / "in.pdb"
    if not sc.exists():
        strip_ligand(pdb, stripped)
        subprocess.run(
            [str(RELAX), "-database", str(DB), "-s", str(stripped),
             "-relax:default_repeats", "1", "-nstruct", "1",
             "-out:path:all", str(d), "-out:file:scorefile", "score.sc",
             "-overwrite", "-ignore_unrecognized_res", "-ex1", "-ex2aro",
             "-mute", "all"],
            capture_output=True, text=True, timeout=2400)
    if not sc.exists():
        return None
    lines = [l for l in sc.read_text().splitlines() if l.startswith("SCORE:")]
    rec = dict(zip(lines[0].split(), lines[-1].split()))
    nres = len(re.findall(r"^ATOM.{9}CA ", stripped.read_text(), re.M))
    tot = float(rec["total_score"])
    return {"tag": tag, "total_score": tot, "n_res": nres,
            "score_per_res": tot / max(nres, 1),
            "fa_rep": float(rec.get("fa_rep", "nan")),
            "fa_rep_per_res": float(rec.get("fa_rep", "nan")) / max(nres, 1),
            "fa_atr": float(rec.get("fa_atr", "nan")),
            "fa_sol": float(rec.get("fa_sol", "nan")),
            "hbond_sc": float(rec.get("hbond_sc", "nan")),
            "p_aa_pp": float(rec.get("p_aa_pp", "nan")),
            "source_pdb": str(pdb)}


def main() -> int:
    top = json.loads((ROOT / "outputs/top20.json").read_text())
    targets = []
    for r in top:
        p = best_seed_pdb(r["id"])
        if p:
            targets.append((p, r["id"]))
    ctrl = ROOT / "outputs/runs/controls_full/4BS0_A.pdb"
    if ctrl.exists():
        targets.append((ctrl, "REF_HG3.17_4BS0"))
    print(f"relaxing {len(targets)} structures (protein only)", flush=True)
    rows = []
    for p, t in targets:
        try:
            r = relax(p, t)
        except Exception as exc:                                  # noqa: BLE001
            print(f"  {t}: {type(exc).__name__}: {exc}", flush=True)
            continue
        if r:
            rows.append(r)
            print(f"  {t:34} {r['n_res']:4d} res  per-res {r['score_per_res']:7.3f}"
                  f"  fa_rep/res {r['fa_rep_per_res']:5.3f}", flush=True)
    (OUT / "relax_top20.json").write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {OUT/'relax_top20.json'} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
