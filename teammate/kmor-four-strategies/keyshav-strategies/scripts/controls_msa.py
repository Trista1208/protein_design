"""The HG3.17 vs HG3.17-Q50M control ranking, re-run WITH an MSA.

WHY: `python -m kemp controls` folded both controls in single-sequence mode and
returned

    positive_median 3.0903   negative_median 2.6409   ranks_correctly FALSE
    "STOP - pipeline cannot rank HG3.17 above HG3.17-Q50M ... every number
     produced afterwards is noise"

and that verdict gates M4 in kemp/cli.py. scripts/control_fold_matrix.py then
showed the folding protocol itself was the problem: without an MSA, Boltz puts
HG3.17's oxyanion contact 12.7 A away and scores the crystallographically
perfect enzyme 0.0; with an MSA it scores 0.775. So the control deserves a
re-run before its STOP verdict is believed.

WHAT IS DIFFERENT FROM `kemp controls`, and nothing else:
  1. --use_msa_server instead of `msa: empty`
  2. the full official score is computed per replicate, not iPDE alone, so the
     geometry and Consistency terms get a say. iPDE is still reported on its
     own for an apples-to-apples comparison with the failed run.

SEQUENCES: taken from the deposited coordinates (4BS0 chain A, 6TU6 chain G),
300 aa each, differing at exactly ONE position - index 47, Gln vs Met. Both are
public PDB entries, so the MSA server sees nothing of ours.

THIS SCRIPT DOES NOT WRITE state.json. The M4 gate stays exactly as it was; a
machine-written verdict should not quietly unblock a pipeline on a protocol
change. Read the output, then decide.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from kemp import config as C
from kemp import official, pipeline, scoring, strain

#: Residues whose side chain can donate the oxyanion H-bond (an NE2/ND2 amide).
OXY_DONORS = {"GLN", "ASN"}


def score_with_donor_check(combined, ipdes, oxy_resname: str) -> dict:
    """Official score, with an ABSENT donor scored 0 rather than 'unmeasurable'.

    measure_raw returns d_oxy=None for the Q50M knockout because Met has no
    NE2 - and official.score() maps None to subscore None, i.e. "could not
    measure". With strict=False those terms are then EXCLUDED and the knockout
    is scored as though the oxyanion requirement did not exist: it came out at
    0.825 against the positive's 0.834, which is nonsense. An absent donor is a
    hard zero, exactly as scripts/controls_full.py:135-157 treats it.

    The total is the weighted GEOMETRIC mean the project uses
    (official.py: "S = prod s_i ** w_i"); the arithmetic sum that official.py
    keeps as a diagnostic flatters a zeroed design to ~0.61.
    """
    import math
    sc = official.score(combined, strict=False)
    donor_absent = oxy_resname not in OXY_DONORS
    if donor_absent:
        sc["subscores"]["oxy"] = 0.0
        sc["subscores"]["orient_oxy"] = 0.0
    subs = sc["subscores"]
    present = {k: v for k, v in subs.items() if v is not None}
    covered = sum(official.WEIGHTS[k] for k in present)
    if any(v <= 0.0 for v in present.values()):
        total = 0.0
    else:
        total = math.exp(sum(official.WEIGHTS[k] * math.log(v)
                             for k, v in present.items()))
    sc["corrected_geometric_total"] = round(total, 6)
    sc["weight_covered"] = round(covered, 4)
    sc["donor_absent"] = donor_absent
    sc["residue_at_oxy"] = oxy_resname
    sc["zeroed_terms"] = [k for k, v in present.items() if v <= 0.0]
    return sc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "runs" / "controls_msa"
CRYSTALS = ROOT / "outputs" / "runs" / "controls_full"

BASE_IDX = 124   # crystal Asp127
OXY_IDX = 47     # crystal Gln50 (positive) / Met50 (negative)

CASES = {
    "positive": {"pdb": "4BS0", "chain": "A", "name": C.CONTROL_POSITIVE["name"]},
    "negative": {"pdb": "6TU6", "chain": "G", "name": C.CONTROL_NEGATIVE["name"]},
}

THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
             "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
             "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
             "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def sequence_from(pdb: Path) -> str:
    atoms = strain.load_atoms(pdb)
    seen, seq = set(), []
    for a in atoms:
        key = (a.chain, a.resseq)
        if key in seen or a.resname not in THREE2ONE:
            continue
        seen.add(key)
        seq.append(THREE2ONE[a.resname])
    return "".join(seq)


def fold_one(wdir: Path, seq: str, seed: int) -> tuple[int, Path | None, dict | None]:
    wdir.mkdir(parents=True, exist_ok=True)
    yaml = wdir / "input.yaml"
    # No pocket constraint: control_fold_matrix.py showed it is worthless with
    # an MSA (0.775 -> 0.756) and cannot rescue the oxyanion without one.
    yaml.write_text("\n".join([
        "version: 1", "sequences:",
        "  - protein:", "      id: A", f'      sequence: "{seq}"',
        "  - ligand:", "      id: B", f"      ccd: {C.LIGAND_CCD}",
    ]) + "\n")
    cmd = [str(C.BIN_BOLTZ), "predict", str(yaml), "--out_dir", str(wdir),
           "--diffusion_samples", str(C.BOLTZ_SETTINGS["diffusion_samples"]),
           "--recycling_steps", str(C.BOLTZ_SETTINGS["recycling_steps"]),
           "--seed", str(seed), "--override", "--use_msa_server"]
    with (wdir / "boltz.log").open("w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                              timeout=5400)
    models = sorted(wdir.rglob("*_model_0.cif"))
    conf = None
    for cj in sorted(wdir.rglob("confidence*.json")):
        try:
            conf = json.loads(cj.read_text())
        except json.JSONDecodeError:
            pass
        break
    return proc.returncode, (models[0] if models else None), conf


def main() -> int:
    print("=" * 72)
    print("CONTROLS WITH MSA   HG3.17 vs HG3.17-Q50M")
    print("=" * 72)

    seqs = {}
    for key, meta in CASES.items():
        pdb = CRYSTALS / f"{meta['pdb']}_{meta['chain']}.pdb"
        s = sequence_from(pdb)
        seqs[key] = s
        print(f"  {key:9} {meta['pdb']}_{meta['chain']}  {len(s)} aa  "
              f"residue {OXY_IDX} = {s[OXY_IDX - 1]}  (base {s[BASE_IDX - 1]}{BASE_IDX})")
    assert seqs["positive"][OXY_IDX - 1] == "Q", "positive must have Gln50"
    assert seqs["negative"][OXY_IDX - 1] == "M", "negative must have Met50"
    diffs = [i for i, (a, b) in enumerate(zip(seqs["positive"], seqs["negative"]), 1)
             if a != b]
    print(f"  sequences differ at {len(diffs)} position(s): {diffs}")
    print()

    results: dict[str, dict] = {}
    for key, seq in seqs.items():
        print(f"--- folding {key} x{C.REPLICATES} with MSA", flush=True)
        per_rep, ipdes, rows = [], [], []
        for seed in C.BOLTZ_SEEDS[:C.REPLICATES]:
            t0 = time.time()
            rc, model, conf = fold_one(OUT / key / f"seed{seed}", seq, seed)
            if rc != 0 or model is None:
                print(f"    seed {seed}: FAILED rc={rc}")
                continue
            pdb = pipeline.cif_to_pdb(model, model.parent / "folded.pdb")
            site = official.CatalyticSite(base_resseq=BASE_IDX,
                                          oxy_resseq=OXY_IDX, chain="A")
            raw = official.measure_raw(str(pdb), site=site)
            per_rep.append(raw)
            ipde = None
            if conf:
                try:
                    _, ipde = scoring.extract_ipde(conf)
                    ipdes.append(ipde)
                except scoring.ScoringError as exc:
                    print(f"    seed {seed}: ipde: {str(exc)[:90]}")
            row = {"seed": seed,
                   "d_base": None if raw.d_base is None else round(raw.d_base, 3),
                   "d_oxy": None if raw.d_oxy is None else round(raw.d_oxy, 3),
                   "ipde": None if ipde is None else round(ipde, 4),
                   "plddt": round(conf.get("complex_plddt"), 4) if conf else None,
                   "wall_s": round(time.time() - t0, 1)}
            rows.append(row)
            print(f"    seed {seed}: d_base={row['d_base']} d_oxy={row['d_oxy']} "
                  f"ipde={row['ipde']} plddt={row['plddt']} ({row['wall_s']}s)",
                  flush=True)

        if not per_rep:
            print(f"  no usable replicates for {key}")
            return 3

        combined = official.combine_replicates(per_rep, ipde_values=ipdes)
        oxy_resname = strain._resname(
            strain.load_atoms(sorted((OUT / key).rglob("folded.pdb"))[0]),
            ("A", OXY_IDX))
        scored = score_with_donor_check(combined, ipdes, oxy_resname)
        ipde_agg = (scoring.aggregate(ipdes, mode="median",
                                      lower_is_better=C.IPDE_LOWER_IS_BETTER)
                    if ipdes else None)
        results[key] = {
            "name": CASES[key]["name"],
            "replicates": rows,
            "d_base_mean": None if combined.d_base is None else round(combined.d_base, 3),
            "d_oxy_mean": None if combined.d_oxy is None else round(combined.d_oxy, 3),
            "sd_base": None if combined.sd_base is None else round(combined.sd_base, 4),
            "subscores": scored.get("subscores"),
            "official_score": scored.get("corrected_geometric_total"),
            "prefilter_score": scored.get("prefilter_score"),
            "residue_at_oxy": scored.get("residue_at_oxy"),
            "donor_absent": scored.get("donor_absent"),
            "zeroed_terms": scored.get("zeroed_terms"),
            "limiting_term": scored.get("limiting_term"),
            "verdict": scored.get("verdict"),
            "ipde": ipde_agg,
        }
        print(f"  {key:9} score={results[key]['official_score']} "
              f"d_base={results[key]['d_base_mean']} d_oxy={results[key]['d_oxy_mean']} "
              f"sd_base={results[key]['sd_base']} ipde_median="
              f"{(ipde_agg or {}).get('median')}\n", flush=True)

    # --- verdict 1: iPDE alone, exactly the comparison that failed before
    ipde_verdict = None
    if results["positive"]["ipde"] and results["negative"]["ipde"]:
        nf = scoring.NoiseFloor(positive=results["positive"]["ipde"],
                                negative=results["negative"]["ipde"])
        ipde_verdict = nf.verdict()

    # --- verdict 2: the full official score
    p = results["positive"]["official_score"]
    n = results["negative"]["official_score"]
    score_verdict = {
        "positive_score": p, "negative_score": n,
        "ranks_correctly": bool(p is not None and n is not None and p > n),
        "separation": (None if p is None or n is None else round(p - n, 6)),
    }

    out = {"protocol": "boltz + MSA server, no pocket constraint, "
                       f"{C.REPLICATES} seeds",
           "cases": results,
           "ipde_verdict": ipde_verdict,
           "score_verdict": score_verdict}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "verdict.json").write_text(json.dumps(out, indent=2, default=str))

    print("=" * 72)
    print("VERDICT")
    print("=" * 72)
    print("\n  iPDE alone (the comparison that previously returned STOP):")
    print("  " + json.dumps(ipde_verdict, indent=2).replace("\n", "\n  "))
    print("\n  Full official score:")
    print("  " + json.dumps(score_verdict, indent=2).replace("\n", "\n  "))
    print("\n  For reference, the same comparison WITHOUT an MSA was:")
    print("    positive_median 3.0903  negative_median 2.6409  "
          "ranks_correctly False")
    print(f"\n  wrote {OUT / 'verdict.json'}")
    print("  state.json NOT written - the M4 gate is unchanged.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
