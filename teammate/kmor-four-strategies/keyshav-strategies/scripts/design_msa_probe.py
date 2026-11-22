"""Does a de novo design have ANY homologs? One sequence, one question.

The control matrix showed MSA depth is the difference between scoring 0.0 and
0.775 on HG3.17 - but HG3.17 is a natural xylanase derivative with thousands of
homologs. A de novo RFD3 design may have none, in which case --use_msa_server
returns an empty alignment and the finding does not transfer to our designs at
all. That is the question here, and it is cheap to answer.

SEQUENCE UNDER TEST: diag_8seq_fixed sequence 4 - the LigandMPNN sequence that
reached d_base 2.557 and was used for the whole foldtest comparison. Chosen so
the result is directly comparable to numbers already measured WITHOUT an MSA:

    none      d_base  9.149   d_oxy 14.738
    pocket5   d_base  4.373   d_oxy  4.761

so any movement is attributable to the MSA and nothing else.

DISCLOSURE: this sends ONE design sequence to the public ColabFold MSA server,
which the user explicitly approved for this single test. No other design is
sent. If the alignment comes back shallow or empty, the question is closed and
no further disclosure is needed.

GPU: waits for any other boltz process to finish first, so this never contends
for VRAM with a running campaign or control job.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from kemp import config as C
from kemp import official, pipeline

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "msaprobe"
SRC_YAML = ROOT / "outputs" / "foldtest" / "pocket5.yaml"

BASE_RESSEQ = 131      # diag_8seq_fixed diffused_index_map G2 -> A131 (Asp)
OXY_RESSEQ = 100       # ...                              G1 -> A100 (Gln)

#: measured previously on this exact sequence, single-sequence mode
NO_MSA = {"none": (9.149, 14.738), "pocket": (4.373, 4.761)}


def sequence_from_yaml(path: Path) -> str:
    m = re.search(r'sequence:\s*"([A-Z]+)"', path.read_text())
    if not m:
        raise SystemExit(f"no sequence found in {path}")
    return m.group(1)


def wait_for_gpu(poll: int = 20, limit: int = 3600) -> None:
    waited = 0
    while waited < limit:
        r = subprocess.run(["pgrep", "-f", "boltz predic[t]"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            return
        print(f"    (another boltz run active; waiting {poll}s)", flush=True)
        time.sleep(poll)
        waited += poll
    print("    WARNING: gave up waiting; proceeding anyway", flush=True)


def msa_depth(wdir: Path) -> tuple[int, str]:
    """Number of sequences in whatever alignment boltz fetched."""
    for a3m in sorted(wdir.rglob("*.a3m")) + sorted(wdir.rglob("*.csv")):
        try:
            text = a3m.read_text(errors="replace")
        except OSError:
            continue
        n = text.count(">") if a3m.suffix == ".a3m" else max(
            0, len(text.strip().splitlines()) - 1)
        return n, str(a3m)
    return -1, "no alignment file found"


def run(tag: str, seq: str, *, pocket: bool) -> dict:
    wdir = OUT / tag
    wdir.mkdir(parents=True, exist_ok=True)
    lines = ["version: 1", "sequences:",
             "  - protein:", "      id: A", f'      sequence: "{seq}"',
             "  - ligand:", "      id: B", f"      ccd: {C.LIGAND_CCD}"]
    if pocket:
        lines += ["constraints:", "  - pocket:", "      binder: B",
                  f"      contacts: [[A, {BASE_RESSEQ}], [A, {OXY_RESSEQ}]]",
                  "      max_distance: 5.0"]
    yaml = wdir / "input.yaml"
    yaml.write_text("\n".join(lines) + "\n")

    wait_for_gpu()
    cmd = [str(C.BIN_BOLTZ), "predict", str(yaml), "--out_dir", str(wdir),
           "--diffusion_samples", str(C.BOLTZ_SETTINGS["diffusion_samples"]),
           "--recycling_steps", str(C.BOLTZ_SETTINGS["recycling_steps"]),
           "--seed", "1", "--override", "--use_msa_server"]
    t0 = time.time()
    with (wdir / "boltz.log").open("w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                              timeout=5400)
    rec = {"tag": tag, "pocket": pocket, "rc": proc.returncode,
           "wall_s": round(time.time() - t0, 1)}
    depth, where = msa_depth(wdir)
    rec["msa_depth"] = depth
    rec["msa_file"] = where

    models = sorted(wdir.rglob("*_model_0.cif"))
    if proc.returncode != 0 or not models:
        rec["error"] = f"no structure (rc={proc.returncode})"
        return rec
    pdb = pipeline.cif_to_pdb(models[0], wdir / "folded.pdb")
    site = official.CatalyticSite(base_resseq=BASE_RESSEQ,
                                  oxy_resseq=OXY_RESSEQ, chain="A")
    raw = official.measure_raw(str(pdb), site=site)
    scored = official.score(raw, strict=False)
    rec.update({
        "d_base": None if raw.d_base is None else round(raw.d_base, 3),
        "d_oxy": None if raw.d_oxy is None else round(raw.d_oxy, 3),
        "subscores": scored.get("subscores"),
        "prefilter_score": scored.get("prefilter_score"),
        "limiting_term": scored.get("limiting_term"),
    })
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
    seq = sequence_from_yaml(SRC_YAML)
    print(f"design sequence under test: {len(seq)} aa")
    print(f"  catalytic: Asp{BASE_RESSEQ} / Gln{OXY_RESSEQ}")
    print(f"  residue {BASE_RESSEQ} = {seq[BASE_RESSEQ - 1]} (expect D), "
          f"residue {OXY_RESSEQ} = {seq[OXY_RESSEQ - 1]} (expect Q)")
    if seq[BASE_RESSEQ - 1] != "D" or seq[OXY_RESSEQ - 1] != "Q":
        print("  WARNING: catalytic residues are not where expected - "
              "indices may be off by the LigandMPNN numbering")
    print()

    rows = []
    for tag, pocket in (("none_msa", False), ("pocket_msa", True)):
        print(f"--- {tag}", flush=True)
        rows.append(run(tag, seq, pocket=pocket))
        print("    " + json.dumps(rows[-1], default=str), flush=True)
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "results.json").write_text(json.dumps(rows, indent=2, default=str))

    print("\n" + "=" * 72)
    print("DOES A DE NOVO DESIGN GET AN MSA?")
    print("=" * 72)
    d = rows[0].get("msa_depth")
    print(f"\n  alignment depth returned by the server: {d} sequences")
    if d is not None and d <= 1:
        print("  -> EMPTY. The design has no detectable homologs, so the MSA")
        print("     finding does NOT transfer to de novo designs. Single-")
        print("     sequence mode is the only honest protocol for them.")
    else:
        print("  -> NON-EMPTY. Compare the geometry below before concluding.")

    print(f"\n{'condition':14} {'d_base':>8} {'d_oxy':>8} "
          f"{'(no-MSA baseline)':>22} {'score':>9}")
    for r in rows:
        key = "pocket" if r["pocket"] else "none"
        b, o = NO_MSA[key]
        print(f"{r['tag']:14} {str(r.get('d_base')):>8} {str(r.get('d_oxy')):>8} "
              f"{f'{b} / {o}':>22} {str(r.get('prefilter_score')):>9}")
    print(f"\nwrote {OUT / 'results.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
