"""THE MSA TEST: does an MSA move rfd4_design_045's oxyanion rotamer?

WHY THIS RUN EXISTS
-------------------
Two protocols are in circulation for the same design lineage and nothing
recorded which produced which number:

  * Theo's RFD4 pipeline: "refolded with 6NT and ColabFold MSAs"
  * scripts/refold_d045.py:50: `msa: empty`  (single sequence)

LESSONS 3.1 shows the MSA setting is not a detail: on HG3.17 it moved d_oxy
from 12.698 A to 2.952 A and the official score from 0.000 to 0.775. That is
the SAME term family that is zeroed on d045. LESSONS 3.2 counters that the
ColabFold server returns depth 1 for a de novo design, which would make the
two protocols identical in effect - but that has never been checked on THIS
sequence. The organisers' own slide (p.7) lists "MSA depth for de novo design
can impact folding confidence" as a limitation, so they know it matters too.

THE QUESTION, precisely: with an MSA, does the oxyanion dihedral leave the
230-310 dead zone? Single-sequence it sits at 233.4-239.8 deg, mean 237.2,
sd 3.06 - reproducible, and missing the 230 band edge by 3.4 deg at best.

WHAT WOULD FALSIFY WHAT
-----------------------
  depth == 1 and dihedral unchanged -> MSA is a no-op for this sequence.
      LESSONS 3.2 confirmed. The 237 deg rotamer is a property of the
      SEQUENCE, so only a sequence change can move it: the loop is on.
  depth  > 1 and dihedral moves in band -> the 3.4 deg problem was partly a
      protocol artefact and the design may already score. Re-examine before
      building anything.
  depth  > 1 and dihedral unchanged -> the rotamer is robust to the fold
      protocol. Strongest possible case for the loop.

NOTE: single-sequence pLDDT is already 0.963-0.967 here, so this is NOT a
foldability test. It asks only whether the MSA changes the Gln rotamer choice.

Everything else is held identical to scripts/refold_d045.py: same sequence,
same seeds, diffusion_samples 1, recycling 3, no pocket constraint. The ONLY
variable is the MSA. That is the point (LESSONS 3.2: record the MSA setting
alongside every score).
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C, objective, official, pipeline   # noqa: E402

SEQ = ("SSFLELMHKAVSAAKQGDFTSAISDFVEACNLLLPDDKKDDFEAMQKALLSTLLAHHDQDRAKDALTTLLTVTLILFE"
       "KHHDDVVELMQKHYDHFHDMVSDDHKIVTEYLKGELDYEDYLTSLVHNRGDALNLSEAEITEIYEEAKKMLS")
BASE_IDX, OXY_IDX = 102, 46
OUT = Path(__file__).resolve().parents[1] / "outputs" / "msa_test_d045"

#: the single-sequence baseline, measured (outputs/refold_d045/verdict.json)
BASELINE = {
    "msa_mode": "none:single_sequence",
    "dih_o": [239.4, 234.4, 239.0, 239.8, 233.4],
    "d_oxy": [3.964, 3.654, 3.816, 9.654, 4.327],
    "d_base": [3.023, 2.864, 2.975, 5.272, 2.882],
    "plddt": [0.963, 0.967, 0.967, 0.967, 0.965],
}


def msa_depth(work: Path) -> tuple[int | None, str]:
    """Number of sequences in whatever MSA Boltz actually used.

    This is the number LESSONS 3.2 turns on. Depth 1 means the server returned
    the query and nothing else, so `--use_msa_server` and `msa: empty` are the
    same run wearing different labels.
    """
    for pat in ("*.a3m", "**/*.a3m", "**/msa/*.csv", "**/*.csv"):
        for f in sorted(work.rglob(pat.replace("**/", ""))) or sorted(work.rglob(pat)):
            try:
                txt = f.read_text(errors="replace")
            except Exception:                                  # noqa: BLE001
                continue
            if f.suffix == ".a3m":
                n = sum(1 for ln in txt.splitlines() if ln.startswith(">"))
                if n:
                    return n, str(f)
            elif f.suffix == ".csv" and "msa" in str(f).lower():
                n = max(0, len([ln for ln in txt.splitlines() if ln.strip()]) - 1)
                if n:
                    return n, str(f)
    return None, "no MSA file found in the run directory"


def main() -> int:
    assert len(SEQ) == 150 and SEQ[BASE_IDX - 1] == "D" and SEQ[OXY_IDX - 1] == "Q"
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"rfd4_design_045  {len(SEQ)} aa  Asp{BASE_IDX}/Gln{OXY_IDX}")
    print("variable: MSA only.  baseline (msa: empty) dihedral_oxy = "
          f"{BASELINE['dih_o']}  -> all in the 230-310 DEAD ZONE\n")

    rows, raws, ipdes = [], [], []
    depth_seen: int | None = None
    for seed in C.BOLTZ_SEEDS[:C.REPLICATES]:
        w = OUT / f"seed{seed}"
        w.mkdir(parents=True, exist_ok=True)
        # NOTE: no `msa: empty` line here - that is the whole experiment.
        (w / "input.yaml").write_text("\n".join([
            "version: 1", "sequences:", "  - protein:", "      id: A",
            f'      sequence: "{SEQ}"',
            "  - ligand:", "      id: B", f"      ccd: {C.LIGAND_CCD}"]) + "\n")
        t0 = time.time()
        with (w / "boltz.log").open("w") as fh:
            p = subprocess.run(
                [str(C.BIN_BOLTZ), "predict", str(w / "input.yaml"),
                 "--out_dir", str(w), "--diffusion_samples", "1",
                 "--recycling_steps", "3", "--seed", str(seed),
                 "--use_msa_server", "--override"],
                stdout=fh, stderr=subprocess.STDOUT, timeout=5400)
        cifs = sorted(w.rglob("*_model_0.cif"))
        if p.returncode != 0 or not cifs:
            tail = (w / "boltz.log").read_text(errors="replace")[-400:]
            print(f"  seed {seed}: FAILED rc={p.returncode}\n    {tail}")
            rows.append({"seed": seed, "failed": True, "rc": p.returncode})
            continue

        depth, src = msa_depth(w)
        if depth_seen is None:
            depth_seen = depth
            print(f"  MSA depth = {depth}  ({src})")
            if depth == 1:
                print("  -> depth 1: the server returned the query only. "
                      "LESSONS 3.2 confirmed; MSA is a no-op for this sequence.")

        pdb = pipeline.cif_to_pdb(cifs[0], w / "folded.pdb")
        raw = official.measure_raw(str(pdb), site=official.CatalyticSite(
            base_resseq=BASE_IDX, oxy_resseq=OXY_IDX, chain="A"))
        raws.append(raw)
        conf = None
        for cj in sorted(w.rglob("confidence*.json")):
            conf = json.loads(cj.read_text())
            break
        if conf and conf.get("complex_ipde") is not None:
            ipdes.append(conf["complex_ipde"])
        gap = objective.dihedral_band_gap(raw.dihedral12_oxy)
        rows.append({
            "seed": seed, "d_base": raw.d_base, "d_oxy": raw.d_oxy,
            "dih_b": raw.dihedral12_base, "dih_o": raw.dihedral12_oxy,
            "gap_dih_oxy_deg": gap, "msa_depth": depth,
            "plddt": (conf or {}).get("complex_plddt"),
            "ipde": (conf or {}).get("complex_ipde"),
            "wall": round(time.time() - t0, 1)})
        r = rows[-1]
        print(f"  seed {seed}: d_base {r['d_base']:.3f}  d_oxy {r['d_oxy']:.3f}  "
              f"dih_o {r['dih_o']:.1f}  gap {gap:.1f} deg  "
              f"pLDDT {r['plddt']:.3f}  iPDE {r['ipde']:.2f}  ({r['wall']}s)",
              flush=True)
        (OUT / "replicates.json").write_text(json.dumps(rows, indent=2, default=str))

    if not raws:
        print("\nNo replicate folded. Nothing can be concluded about the MSA.")
        return 3

    comb = official.combine_replicates(raws, ipde_values=ipdes)
    sc = official.score(comb, strict=True)
    got = [r["dih_o"] for r in rows if not r.get("failed")]
    gaps = [r["gap_dih_oxy_deg"] for r in rows if not r.get("failed")]
    n_in_band = sum(1 for g in gaps if g == 0.0)

    print("\n" + "=" * 70)
    print(f"MSA depth: {depth_seen}")
    print(f"  with MSA   dihedral_oxy: {[round(x,1) for x in got]}")
    print(f"  no MSA     dihedral_oxy: {BASELINE['dih_o']}")
    print(f"  band gaps  with MSA:     {[round(g,1) for g in gaps]}")
    print(f"  seeds in band: {n_in_band}/{len(got)}")
    print(f"  subscores: " + " ".join(
        f"{k}={round(v,3) if v is not None else None}"
        for k, v in sc["subscores"].items()))
    print(f"  official_score = {sc.get('official_score')}  "
          f"limiting = {sc.get('limiting_term')}")

    if n_in_band:
        verdict = ("MSA MOVES THE ROTAMER IN BAND. The 3.4 deg miss was partly a "
                   "protocol artefact - re-examine before building the loop.")
    elif depth_seen == 1:
        verdict = ("MSA is a NO-OP (depth 1). The 237 deg rotamer is a property of "
                   "the SEQUENCE, so only a sequence change moves it. Loop is on.")
    else:
        verdict = (f"MSA had real depth ({depth_seen}) and the rotamer did NOT move. "
                   "The rotamer is robust to fold protocol. Strongest case for the loop.")
    print("\n  VERDICT: " + verdict)

    (OUT / "verdict.json").write_text(json.dumps({
        "protocol": objective.protocol_stamp(
            replicates=len(got), extra={"msa_mode": "msa_server",
                                        "msa_depth": depth_seen}),
        "baseline_single_sequence": BASELINE,
        "replicates": rows,
        "subscores": sc["subscores"],
        "official_score": sc.get("official_score"),
        "n_in_band": n_in_band,
        "verdict": verdict,
    }, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
