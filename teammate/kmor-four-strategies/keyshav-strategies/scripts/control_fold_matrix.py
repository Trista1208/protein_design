"""Can Boltz reproduce HG3.17's catalytic site from HG3.17's own sequence?

The crystal coordinates score 0.959 (outputs/runs/controls_full/verdict.json),
so the 2.0-3.6 A window is satisfiable on both contacts at once by a real
backbone. What is NOT yet known is whether the FOLDING protocol can recover
that site. Until it is, a de novo design failing by ~1 A cannot be attributed
to the design rather than to the evaluator.

Four folds, one variable pair:

    constraint x {none, pocket@5.0}  x  MSA x {empty, server}

Everything else is the project's pinned setting (recycling_steps 3,
diffusion_samples 1, seed 1) so the numbers are comparable to the design runs.

Outcomes, both decisive:
  * Boltz recovers ~2.5 / ~3.1  -> protocol validated; the ~1 A gap is the
    de novo backbone, and a physically valid scaffold is worth pursuing.
  * Boltz cannot, even from HG3.17's own sequence -> the gap is an evaluator
    limitation and no amount of design work closes it.

NOTE ON THE MSA ARM: --use_msa_server sends the query sequence to the public
ColabFold MSA server. The sequence used here is HG3.17 from PDB 4BS0, which is
already public; no design of ours is sent anywhere.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from kemp import config as C
from kemp import official, pipeline, strain

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "controlfold"
CRYSTAL = ROOT / "outputs" / "runs" / "controls_full" / "4BS0_A.pdb"

# HG3.17 (4BS0 chain A), 300 aa. Boltz renumbers 1..N, so the catalytic
# residues are at SEQUENCE indices, not crystal resseq:
#   crystal Asp127 -> index 124      crystal Gln50 -> index 47
# Verified by walking 4BS0_A.pdb; asserted below rather than trusted.
BASE_IDX = 124
OXY_IDX = 47

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


def write_yaml(path: Path, seq: str, *, pocket: bool, msa: bool) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["version: 1", "sequences:",
             "  - protein:", "      id: A", f'      sequence: "{seq}"']
    if not msa:
        lines.append("      msa: empty")
    lines += ["  - ligand:", "      id: B", f"      ccd: {C.LIGAND_CCD}"]
    if pocket:
        lines += ["constraints:",
                  "  - pocket:", "      binder: B",
                  f"      contacts: [[A, {BASE_IDX}], [A, {OXY_IDX}]]",
                  "      max_distance: 5.0"]
    path.write_text("\n".join(lines) + "\n")
    return path


def run_fold(tag: str, seq: str, *, pocket: bool, msa: bool) -> dict:
    wdir = OUT / tag
    wdir.mkdir(parents=True, exist_ok=True)
    inp = write_yaml(wdir / "input.yaml", seq, pocket=pocket, msa=msa)
    cmd = [str(C.BIN_BOLTZ), "predict", str(inp), "--out_dir", str(wdir),
           "--diffusion_samples", str(C.BOLTZ_SETTINGS["diffusion_samples"]),
           "--recycling_steps", str(C.BOLTZ_SETTINGS["recycling_steps"]),
           "--seed", "1", "--override"]
    if msa:
        cmd.append("--use_msa_server")
    t0 = time.time()
    log = wdir / "boltz.log"
    with log.open("w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                              timeout=5400)
    wall = time.time() - t0
    rec = {"tag": tag, "pocket": pocket, "msa": msa, "rc": proc.returncode,
           "wall_s": round(wall, 1)}

    cifs = [p for p in sorted(wdir.rglob("*_model_0.cif"))]
    if proc.returncode != 0 or not cifs:
        rec["error"] = f"no structure (rc={proc.returncode}); see {log}"
        return rec
    model = cifs[0]
    rec["model"] = str(model)

    # measure with the project's own code, so the numbers sit on the same
    # scale as outputs/runs/controls_full/verdict.json.
    # measure_raw reads PDB, not CIF - convert first, exactly as the pipeline
    # does. Passing the CIF straight in raises TheozymeError on the first atom.
    pdb = pipeline.cif_to_pdb(model, wdir / "folded.pdb")
    site = official.CatalyticSite(base_resseq=BASE_IDX, oxy_resseq=OXY_IDX,
                                  chain="A")
    raw = official.measure_raw(str(pdb), site=site)
    scored = official.score(raw, strict=False)
    rec.update({
        "d_base": None if raw.d_base is None else round(raw.d_base, 3),
        "d_oxy": None if raw.d_oxy is None else round(raw.d_oxy, 3),
        "base": (scored.get("subscores") or {}).get("base"),
        "oxy": (scored.get("subscores") or {}).get("oxy"),
        "prefilter_score": scored.get("prefilter_score"),
        "limiting_term": scored.get("limiting_term"),
        "notes": raw.notes[:3],
    })
    integ = strain.check_chain_integrity(strain.load_atoms(model),
                                        motif_resseq={BASE_IDX, OXY_IDX})
    rec["cn_mean_A"] = integ["cn_mean_A"]
    rec["n_nonphysical"] = integ["n_nonphysical"]

    for cj in sorted(wdir.rglob("confidence*.json")):
        try:
            conf = json.loads(cj.read_text())
            rec["plddt"] = round(conf.get("complex_plddt") or -1, 4)
            rec["ipde"] = round(conf.get("complex_ipde") or -1, 4)
        except json.JSONDecodeError:
            pass
        break
    return rec


def main() -> int:
    seq = sequence_from(CRYSTAL)
    assert len(seq) == 300, f"expected 300 aa, got {len(seq)}"
    assert seq[BASE_IDX - 1] == "D", f"index {BASE_IDX} is {seq[BASE_IDX - 1]}, not D"
    assert seq[OXY_IDX - 1] == "Q", f"index {OXY_IDX} is {seq[OXY_IDX - 1]}, not Q"
    print(f"HG3.17 (4BS0_A): {len(seq)} aa, "
          f"Asp at {BASE_IDX}, Gln at {OXY_IDX} - verified\n")

    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for pocket in (False, True):
        for msa in (False, True):
            tag = f"{'pocket' if pocket else 'none'}_{'msa' if msa else 'nomsa'}"
            print(f"--- {tag}: constraint={'pocket@5.0' if pocket else 'none'}, "
                  f"msa={'server' if msa else 'empty'}", flush=True)
            try:
                rec = run_fold(tag, seq, pocket=pocket, msa=msa)
            except subprocess.TimeoutExpired:
                rec = {"tag": tag, "pocket": pocket, "msa": msa,
                       "error": "timeout"}
            rows.append(rec)
            print("    " + json.dumps(rec, default=str), flush=True)
            (OUT / "results.json").write_text(json.dumps(rows, indent=2,
                                                         default=str))

    print("\n\n================ CONTROL FOLD MATRIX ================")
    print("crystal 4BS0_A (reference): d_base 2.532  d_oxy 3.102  "
          "base 0.989  oxy 0.996  score 0.959")
    print(f"{'condition':18} {'d_base':>8} {'d_oxy':>8} {'base':>7} {'oxy':>7} "
          f"{'score':>7} {'in window':>22}")
    for r in rows:
        if r.get("error"):
            print(f"{r['tag']:18} {'ERROR: ' + r['error'][:60]}")
            continue
        db, do = r.get("d_base"), r.get("d_oxy")
        inw = []
        for nm, v in (("base", db), ("oxy", do)):
            inw.append(f"{nm}:{'YES' if v is not None and 2.0 <= v <= 3.6 else 'no'}")
        print(f"{r['tag']:18} {str(db):>8} {str(do):>8} {str(r.get('base')):>7} "
              f"{str(r.get('oxy')):>7} {str(r.get('prefilter_score')):>7} "
              f"{'  '.join(inw):>22}")
    print(f"\nwrote {OUT / 'results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
