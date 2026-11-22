"""Single-command entry points for every milestone.

    python -m kemp m0        # audit + controls      (CPU, ~2h)
    python -m kemp m1        # RFD3 smoke + VRAM     (GPU, ~2h)
    python -m kemp m2        # LigandMPNN + Boltz    (GPU, ~2h)
    python -m kemp m3        # scoring + gates tests (CPU, ~3h)
    python -m kemp m4        # scale                 (GPU, overnight)
    python -m kemp m5        # select final 5        (CPU, ~2h)
    python -m kemp all       # M0 -> M5, gated

Each milestone GATES the next (handover section 4): M4 refuses to start unless
M0's control pair ranked correctly.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import asdict
from pathlib import Path

from . import config as C
from . import gates, novelty, official, pipeline, rfd3, scoring, theozyme
from .theozyme import load, write_clean_pdb


STATE = C.OUTPUTS / "state.json"


# --------------------------------------------------------------------------- state
def _read_state() -> dict:
    if STATE.exists():
        try:
            return json.loads(STATE.read_text())
        except json.JSONDecodeError:
            pass
    return {}


def _write_state(**kw) -> dict:
    st = _read_state()
    st.update(kw)
    STATE.parent.mkdir(parents=True, exist_ok=True)
    STATE.write_text(json.dumps(st, indent=2, default=str))
    return st


def _banner(tag: str, text: str) -> None:
    print(f"\n{'=' * 72}\n{tag}  {text}\n{'=' * 72}", flush=True)


# --------------------------------------------------------------------------- M0
def m0(args) -> int:
    _banner("M0", "AUDIT & CONTROLS")
    out = C.RUNS / "m0"
    out.mkdir(parents=True, exist_ok=True)
    report: dict = {"milestone": "M0", "started": time.strftime("%F %T")}

    # 1. load both files, prove same motif / different frame
    tz_with = load(C.THEOZYME_WITH_ORI)
    tz_no = load(C.THEOZYME_NO_ORI)
    report["altlocs_withORI"] = sorted(tz_with.altlocs)
    report["ligand_atoms_survived"] = len(tz_with.ligand)
    report["frames"] = theozyme.compare_frames(tz_with, tz_no)
    print(f"  shared atoms {report['frames']['n_shared_atoms']}, "
          f"superposed RMSD {report['frames']['superposed_rmsd_A']:.4f} A, "
          f"raw RMSD {report['frames']['raw_rmsd_A']:.2f} A")

    # 2. ORI sanity in BOTH frames - the mixed-frame trap
    ori_ok = {}
    for tz in (tz_with, tz_no):
        name = tz.path.name
        try:
            d = theozyme.assert_ori_sane(tz, theozyme.ori_for(tz.path))
            ori_ok[name] = {"ok": True, "dist_A": round(d, 3)}
        except theozyme.TheozymeError as exc:
            ori_ok[name] = {"ok": False, "error": str(exc)}
    # and prove the mixed-frame failure is caught
    try:
        theozyme.assert_ori_sane(tz_no, theozyme.ori_for(C.THEOZYME_WITH_ORI))
        ori_ok["mixed_frame_detected"] = False
    except theozyme.TheozymeError:
        ori_ok["mixed_frame_detected"] = True
    report["ori"] = ori_ok
    print(f"  ORI checks: {json.dumps(ori_ok)}")

    # 3. geometry vs the verified reference table
    measured = theozyme.measure_geometry(tz_with)
    report["geometry"] = theozyme.diff_reference(measured)
    for k, v in report["geometry"].items():
        flag = " (MODELLED ref)" if v["modelled_not_observed"] else ""
        print(f"  {k:28} {v['measured']:>8} vs {v['reference']:>7}  "
              f"dev {v['deviation']:+.3f}{flag}")

    # 4. clean PDB for everything downstream
    clean = write_clean_pdb(tz_with, out / "theozyme_clean.pdb")
    report["clean_pdb"] = str(clean)
    print(f"  wrote ORI-free single-conformer PDB -> {clean}")

    # 5. catalytic burial, all three anchor conventions
    report["burial"] = scoring.catalytic_burial(tz_with)
    print(f"  burial anchors: {report['burial']['anchors']}")

    # 6. control pair - requires Boltz; reported as pending when absent
    report["controls"] = {
        "positive": C.CONTROL_POSITIVE, "negative": C.CONTROL_NEGATIVE,
        "status": "PENDING - fetch 4BS0 and 6TU6 sequences, run boltz_replicates "
                  "on both, then kemp m0-controls. M4 is gated on this.",
    }
    print("  controls: PENDING (see docs/USER_GUIDE.md section M0)")

    (out / "report.json").write_text(json.dumps(report, indent=2, default=str))
    _write_state(m0=report)
    print(f"\n  report -> {out / 'report.json'}")
    return 0


# --------------------------------------------------------------------------- M1
def m1(args) -> int:
    _banner("M1", "RFD3 SMOKE TEST (measures VRAM + wall clock)")
    t0 = time.time()
    res = rfd3.smoke(timesteps=args.timesteps)
    res["wall_clock_s"] = round(time.time() - t0, 1)
    print(json.dumps(res, indent=2, default=str))
    _write_state(m1=res)
    return 0 if res.get("returncode") == 0 else 1


# --------------------------------------------------------------------------- M3
def m3(args) -> int:
    _banner("M3", "SCORING + HARD GATES (self-test)")
    tz = load(C.THEOZYME_WITH_ORI)
    measured = theozyme.measure_geometry(tz)
    geo = scoring.geometry_score(measured)
    print("  geometry score on the reference theozyme (should be ~1.0):")
    print(f"    {geo['score']}")
    for name, t in geo["terms"].items():
        print(f"    {name:14} {t}")

    seq = "M" * 240
    g = gates.run_all(seq, motif_residues=2)
    print(f"\n  gate run on a dummy 240-aa sequence: {g['verdict']}")
    for n, d in g["gates"].items():
        print(f"    {n:20} passed={d['passed']} value={d['value']}")

    long_seq = "M" * 300
    g2 = gates.run_all(long_seq)
    assert not g2["passed"], "300-aa sequence must be disqualified"
    print(f"  300-aa sequence correctly {g2['verdict']}")

    agg = scoring.aggregate([0.81, 0.78, 0.84, 0.80, 0.79])
    print(f"\n  aggregate (median-minus-spread): {agg}")
    try:
        scoring.aggregate([0.1], mode="best")
    except scoring.ScoringError as exc:
        print(f"  'best' aggregation correctly refused: {exc}")
    _write_state(m3={"geometry_score": geo["score"], "gates_selftest": "ok"})
    return 0


# --------------------------------------------------------------------------- M4 gate
def m4(args) -> int:
    n = getattr(args, "n", 50)
    _banner("M4", f"SCALE  (n={n})")
    ctrl = _read_state().get("controls") or {}
    if not ctrl.get("ranks_correctly"):
        print("  REFUSING TO RUN.")
        print("  The control pair has not passed. Handover 3.2: if the pipeline")
        print("  cannot rank HG3.17 above HG3.17-Q50M, Geometry and Consistency")
        print("  are not measuring oxyanion stabilisation and every number")
        print("  produced afterwards is noise.")
        print("  Run:  python -m kemp controls")
        if not getattr(args, "force", False):
            return 2
        print("  --force: proceeding UNCALIBRATED.")

    recs, t0 = [], time.time()
    for rec in pipeline.scale(n, timesteps=getattr(args, "timesteps", 200),
                              seqs_per_backbone=getattr(args, "seqs", 8)):
        recs.append(rec)
        el = sum(1 for r in recs if r.eligible)
        print(f"  [{len(recs)}/{n}] {rec.design_id} stage={rec.stage} "
              f"score={rec.score} eligible={rec.eligible} "
              f"({el} eligible, {time.time()-t0:.0f}s)", flush=True)
    _write_state(m4={"n": n, "eligible": sum(1 for r in recs if r.eligible),
                     "wall_s": round(time.time()-t0, 1)})
    if getattr(args, "finalise", False):
        print(f"\n  finalising top {args.top} with Boltz x{C.REPLICATES} ...")
        fin = pipeline.finalise(recs, top_n=args.top)
        print(f"  {len(fin)} finalists scored at {C.REPLICATES} replicates")
    return 0


# --------------------------------------------------------------------------- plumbing
def m2(args) -> int:
    _banner("M2", "LIGANDMPNN + BOLTZ (catalytic residues FIXED)")
    recs = list(pipeline.scale(getattr(args, "n", 1) or 1,
                               out_root=C.RUNS / "m2",
                               timesteps=getattr(args, "timesteps", 200),
                               seqs_per_backbone=getattr(args, "seqs", 8)))
    for r in recs:
        print(f"  {r.design_id} stage={r.stage} len={r.length} "
              f"score={r.score} eligible={r.eligible} catalytic={r.catalytic}")
        if r.official:
            print(f"    official: {r.official.get('verdict')} "
                  f"subscores={r.official.get('subscores')}")
        if r.novelty:
            print(f"    novelty: {r.novelty.get('verdict')} "
                  f"{r.novelty.get('reasons', [])[:1]}")
    _write_state(m2=[asdict(r) for r in recs])
    return 0 if recs and recs[0].stage.startswith("scored") else 1


def controls(args) -> int:
    _banner("CONTROLS", "HG3.17 vs HG3.17-Q50M  (gates M4)")
    out = C.RUNS / "controls"; out.mkdir(parents=True, exist_ok=True)
    corpus = C.ROOT / "inputs" / "corpus" / "seqs"
    seqs = {}
    for name, pdb in (("positive", C.CONTROL_POSITIVE["pdb"]),
                      ("negative", C.CONTROL_NEGATIVE["pdb"])):
        f = corpus / f"{pdb}.fasta"
        if not f.exists():
            print(f"  MISSING {f} - run scripts/install_novelty.sh first")
            return 2
        best, cur = "", ""
        for line in f.read_text().splitlines():
            if line.startswith(">"):
                best = cur if len(cur) > len(best) else best; cur = ""
            else:
                cur += line.strip()
        best = cur if len(cur) > len(best) else best
        seqs[name] = best
        print(f"  {name:8} {pdb}  {len(best)} aa")

    results = {}
    for name, seq in seqs.items():
        print(f"  folding {name} x{C.REPLICATES} ...", flush=True)
        runs = pipeline.boltz_replicates(seq, out / name)
        vals = []
        for r in runs:
            if r.confidence:
                try:
                    _, v = scoring.extract_ipde(r.confidence)
                    vals.append(v)
                except scoring.ScoringError as exc:
                    print(f"    seed {r.seed}: {str(exc)[:110]}")
        if not vals:
            print(f"  no usable confidence for {name}; cannot calibrate")
            return 3
        results[name] = scoring.aggregate(
            vals, mode="median", lower_is_better=C.IPDE_LOWER_IS_BETTER)
        print(f"  {name:8} median={results[name]['median']} "
              f"spread={results[name]['spread']}")

    nf = scoring.NoiseFloor(positive=results["positive"],
                            negative=results["negative"])
    v = nf.verdict()
    print("\n  " + json.dumps(v, indent=2).replace("\n", "\n  "))
    _write_state(controls={**v, "status": "DONE"})
    (out / "verdict.json").write_text(json.dumps(v, indent=2))
    return 0 if v["ranks_correctly"] else 4


def m5(args) -> int:
    _banner("M5", "SELECT FINAL DESIGNS")
    rows = pipeline.read_ledger()
    if not rows:
        print("  ledger empty - run m4 first"); return 2
    recs = []
    for row in rows:
        row = {k: v for k, v in row.items() if k != "ts"}
        try:
            recs.append(pipeline.DesignRecord(**row))
        except TypeError:
            pass
    nf = ((_read_state().get("controls") or {}).get("noise_floor")) or 0.0
    sel = pipeline.select_final(recs, n=getattr(args, "n", C.N_LAB_SLOTS),
                                noise_floor=nf)
    print(json.dumps({k: v for k, v in sel.items() if k != "records"}, indent=2))
    out = C.RUNS / "m5"; out.mkdir(parents=True, exist_ok=True)
    (out / "selection.json").write_text(json.dumps(sel, indent=2, default=str))
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser("kemp", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("m0", help="audit theozyme + controls")
    p1 = sub.add_parser("m1", help="RFD3 smoke test")
    p1.add_argument("--timesteps", type=int, default=20)
    p2 = sub.add_parser("m2", help="LigandMPNN + Boltz smoke")
    p2.add_argument("--n", type=int, default=1)
    p2.add_argument("--timesteps", type=int, default=200)
    p2.add_argument("--seqs", type=int, default=8)
    sub.add_parser("controls", help="score the control pair (unlocks m4)")
    sub.add_parser("m3", help="scoring + gates self-test")
    p4 = sub.add_parser("m4", help="scale (gated on controls)")
    p4.add_argument("--n", type=int, default=50)
    p4.add_argument("--timesteps", type=int, default=200)
    p4.add_argument("--seqs", type=int, default=8)
    p4.add_argument("--finalise", action="store_true")
    p4.add_argument("--top", type=int, default=30)
    p4.add_argument("--force", action="store_true")
    p5 = sub.add_parser("m5", help="select final designs")
    p5.add_argument("--n", type=int, default=C.N_LAB_SLOTS)
    sub.add_parser("all", help="M0 -> M5, gated")

    args = p.parse_args(argv)
    table = {"m0": m0, "m1": m1, "m2": m2, "m3": m3, "m4": m4, "m5": m5,
             "controls": controls}
    if args.cmd == "all":
        for name in ("m0", "m3", "m1", "m2", "controls", "m4", "m5"):
            fn = table[name]
            rc = fn(args)
            if rc != 0:
                print(f"\n  {name} returned {rc}; stopping the chain.")
                return rc
        return 0
    return table[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
