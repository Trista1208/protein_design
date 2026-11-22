"""Rosetta relax as a fold-stability proxy for the submission set.

WHAT THIS ANSWERS
-----------------
A teammate's concern: "they lack any Cysteine bond, so they are unstable and
won't express or fold like predicted."

Disulfides are one route to stability, not the only one, and they do not form
in a reducing cytoplasm. What actually matters is whether the fold is
well-packed. Rosetta's total_score per residue measures exactly that, so we
measure it - on the designs AND on a reference whose stability is published,
because a number with no control arm is not evidence (good_practices 5b).

PROTEIN ONLY. The ligand is removed: the question is fold stability, and
scoring a ligand without its params file would silently drop it anyway.
"""
from __future__ import annotations
import re, subprocess, sys, json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/stability"
RELAX = Path.home() / "opt/rosetta/source/bin/relax.default.linuxgccrelease"
DB = Path.home() / "opt/rosetta/database"
SUB = ROOT / "published/submission"


def strip_ligand(src: Path, dst: Path) -> int:
    n = 0
    with dst.open("w") as f:
        for ln in src.read_text().splitlines():
            if ln.startswith(("ATOM", "TER")):
                f.write(ln + "\n")
                n += ln.startswith("ATOM")
        f.write("END\n")
    return n


def relax(pdb: Path, tag: str) -> dict | None:
    d = OUT / tag
    d.mkdir(parents=True, exist_ok=True)
    stripped = d / "in.pdb"
    na = strip_ligand(pdb, stripped)
    cmd = [str(RELAX), "-database", str(DB), "-s", str(stripped),
           "-relax:default_repeats", "1", "-nstruct", "1",
           "-out:path:all", str(d), "-out:file:scorefile", "score.sc",
           "-overwrite", "-ignore_unrecognized_res", "-ex1", "-ex2aro",
           "-mute", "all"]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=1500)
    sc = d / "score.sc"
    if not sc.exists():
        print(f"  {tag}: relax produced no scorefile ({p.stderr[-200:]})", flush=True)
        return None
    lines = [l for l in sc.read_text().splitlines() if l.startswith("SCORE:")]
    hdr = lines[0].split()
    row = lines[-1].split()
    rec = dict(zip(hdr, row))
    nres = len(re.findall(r"^ATOM.{9}CA ", stripped.read_text(), re.M))
    tot = float(rec["total_score"])
    return {"tag": tag, "total_score": tot, "n_res": nres,
            "score_per_res": tot / max(nres, 1),
            "fa_rep": float(rec.get("fa_rep", "nan")),
            "fa_atr": float(rec.get("fa_atr", "nan")),
            "hbond_sc": float(rec.get("hbond_sc", "nan")),
            "n_atoms": na}


def main() -> int:
    targets: list[tuple[Path, str]] = []
    sub = json.loads((SUB / "submission.json").read_text())
    for r in sub["selected"]:
        p = SUB / f"{r['design_id']}.pdb"
        if p.exists():
            targets.append((p, r["design_id"]))
    # CONTROL ARM: the reference enzyme crystal structure. Its stability is
    # published (HG3.17, Tm ~70 C), so it calibrates what a number means.
    for ref in [ROOT / "outputs/runs/controls_full/4BS0_A.pdb",
                ROOT / "outputs/runs/controls_full/4BS0_A_bstamp.pdb"]:
        if ref.exists():
            targets.append((ref, "REF_HG3.17_4BS0"))
            break
    print(f"relaxing {len(targets)} structure(s), protein only", flush=True)
    rows = []
    for p, t in targets:
        try:
            r = relax(p, t)
        except Exception as exc:                                  # noqa: BLE001
            print(f"  {t}: {type(exc).__name__}: {exc}", flush=True)
            continue
        if r:
            rows.append(r)
            print(f"  {t:28} {r['n_res']:4d} res  total {r['total_score']:9.1f}"
                  f"  per-res {r['score_per_res']:7.3f}  fa_rep {r['fa_rep']:7.2f}",
                  flush=True)
    (OUT / "relax.json").write_text(json.dumps(rows, indent=2))
    if rows:
        print("\nper-residue Rosetta score (lower = better packed):")
        for r in sorted(rows, key=lambda x: x["score_per_res"]):
            print(f"  {r['tag']:28} {r['score_per_res']:7.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
