"""Does biasing LigandMPNN away from Ala and toward aromatics fix foldability?

THE HYPOTHESIS: our RFD3+LigandMPNN sequences are 24-27% Ala, 1.5-2% aromatic,
zero Trp, and fold at pLDDT 0.41-0.47. Designs that fold well (HG3.17 11% Ala /
9% arom; rfd4_design_045 8.7% / 7.3%) are not. Aromatics pack cores; a
near-poly-alanine chain gives a predictor nothing to lock onto.

THE CONTROL: identical backbone, identical fixed catalytic residues, identical
seed and temperature. ONLY --bias_AA changes. Baseline sequences already exist
(outputs/diag_8seq_fixed/lmpnn) and were folded (boltz_s0..s7), so the
comparison is like-for-like.

WHY IT MATTERS EVEN IF WE ABANDON RFD3: the answer is generator-independent.
RFD4's own inverse folding reports ala_fraction 0.2674 - the same pathology - so
"bias or filter on composition" would transfer to that pipeline too.

OUTCOME: if pLDDT rises toward 0.8, composition was the cause and this is a
one-parameter fix. If it does not move, the backbone is the problem and no
sequence rescues it - equally worth knowing, because it says stop tuning MPNN.
"""
from __future__ import annotations
import collections, json, math, subprocess, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C, runners

ROOT = Path(__file__).resolve().parents[1]
BACKBONE = ROOT / "outputs/diag_8seq_fixed/backbone.pdb"
OUT = ROOT / "outputs/mpnn_bias"
FIXED = "A131 A100"          # catalytic Asp131 / Gln100, same as the baseline run
BASE_IDX, OXY_IDX = 131, 100
# push away from Ala/Gly, toward aromatics and larger hydrophobics
BIAS = {"A": -2.0, "G": -1.0, "W": 1.5, "F": 1.5, "Y": 1.5, "L": 0.5, "I": 0.5, "M": 0.3}
AROM = set("FWY")

def compose(s):
    c = collections.Counter(s); L = len(s)
    ent = -sum((v/L)*math.log(v/L, 20) for v in c.values())
    return dict(n=L, ala=100*c["A"]/L, gly=100*c["G"]/L,
                arom=100*sum(c[a] for a in AROM)/L, trp=c["W"], ent=ent)

def read_fa(p: Path):
    out, hdr = [], None
    for line in p.read_text().splitlines():
        line = line.strip()
        if line.startswith(">"): hdr = line
        elif line: out.append((hdr, line))
    return out[1:]          # first entry is the input backbone sequence

def main() -> int:
    if not BACKBONE.exists():
        print(f"missing backbone {BACKBONE}"); return 2
    OUT.mkdir(parents=True, exist_ok=True)
    bias_str = ",".join(f"{k}:{v}" for k, v in BIAS.items())
    print(f"backbone: {BACKBONE}\nfixed: {FIXED}\nbias: {bias_str}\n")

    ckpt = str(C.LIGANDMPNN_DIR / "model_params" / "ligandmpnn_v_32_010_25.pt")
    wdir = OUT / "biased"; wdir.mkdir(parents=True, exist_ok=True)
    cmd = [str(C.PY_AIZYMES), str(C.LIGANDMPNN_DIR / "run.py"),
           "--model_type", "ligand_mpnn", "--checkpoint_ligand_mpnn", ckpt,
           "--pdb_path", str(BACKBONE), "--out_folder", str(wdir),
           "--batch_size", "8", "--number_of_batches", "1",
           "--temperature", "0.1", "--seed", "37",
           "--fixed_residues", FIXED,
           # LigandMPNN wants "A:-2.0,G:-1.0,..." - NOT json. Passing json
           # raises ValueError: could not convert string to float: ' 0.3}'
           "--bias_AA", ",".join(f"{k}:{v}" for k, v in BIAS.items())]
    with (wdir/"ligandmpnn.log").open("w") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                           cwd=str(C.LIGANDMPNN_DIR), timeout=1800)
    if p.returncode != 0:
        print(f"LigandMPNN failed rc={p.returncode}; see {wdir/'ligandmpnn.log'}"); return 3
    fa = sorted(wdir.rglob("*.fa"))
    if not fa: print("no sequences produced"); return 3
    biased = read_fa(fa[0])
    baseline = read_fa(ROOT/"outputs/diag_8seq_fixed/lmpnn/seqs/backbone.fa")
    print(f"{'set':10} {'#':>2} {'len':>4} {'Ala%':>6} {'Gly%':>6} {'arom%':>6} {'Trp':>4} {'entropy':>8}")
    for lbl, seqs in (("baseline", baseline), ("biased", biased)):
        for i,(h,s) in enumerate(seqs):
            c = compose(s)
            print(f"{lbl:10} {i:>2} {c['n']:>4} {c['ala']:6.1f} {c['gly']:6.1f} {c['arom']:6.1f} {c['trp']:>4} {c['ent']:8.3f}")
    for lbl, seqs in (("baseline", baseline), ("biased", biased)):
        cs = [compose(s) for _, s in seqs]
        print(f"\n{lbl} mean: Ala {sum(c['ala'] for c in cs)/len(cs):.1f}%  "
              f"arom {sum(c['arom'] for c in cs)/len(cs):.1f}%  "
              f"Trp {sum(c['trp'] for c in cs)/len(cs):.1f}  "
              f"entropy {sum(c['ent'] for c in cs)/len(cs):.3f}")
    # verify the catalytic residues survived the bias
    bad = [i for i,(h,s) in enumerate(biased)
           if s[BASE_IDX-1] != "D" or s[OXY_IDX-1] != "Q"]
    print(f"\ncatalytic residues preserved in all biased sequences: {not bad}"
          + (f"  (violations: {bad})" if bad else ""))

    # fold the first 4 biased sequences, same protocol as the baseline folds
    print("\nfolding 4 biased sequences (Boltz, single-seq, seed 1)...", flush=True)
    rows=[]
    for i,(h,s) in enumerate(biased[:4]):
        w = OUT/f"fold{i}"; w.mkdir(parents=True, exist_ok=True)
        (w/"input.yaml").write_text("\n".join([
            "version: 1","sequences:","  - protein:","      id: A",
            f'      sequence: "{s}"',"      msa: empty",
            "  - ligand:","      id: B",f"      ccd: {C.LIGAND_CCD}"])+"\n")
        t0=time.time()
        with (w/"boltz.log").open("w") as fh:
            q=subprocess.run([str(C.BIN_BOLTZ),"predict",str(w/"input.yaml"),"--out_dir",str(w),
                "--diffusion_samples","1","--recycling_steps","3","--seed","1","--override"],
                stdout=fh,stderr=subprocess.STDOUT,timeout=5400)
        conf=None
        for cj in sorted(w.rglob("confidence*.json")):
            conf=json.loads(cj.read_text()); break
        c=compose(s)
        rows.append(dict(i=i, ala=round(c['ala'],1), arom=round(c['arom'],1),
                         plddt=(conf or {}).get("complex_plddt"),
                         ipde=(conf or {}).get("complex_ipde"), wall=round(time.time()-t0,1)))
        print(f"  seq {i}: Ala {c['ala']:.1f}%  arom {c['arom']:.1f}%  "
              f"pLDDT {(rows[-1]['plddt'] or 0):.3f}  iPDE {(rows[-1]['ipde'] or 0):.2f}  ({rows[-1]['wall']}s)",
              flush=True)
    (OUT/"results.json").write_text(json.dumps({"bias":BIAS,"folds":rows},indent=2,default=str))
    pl=[r["plddt"] for r in rows if r["plddt"]]
    if pl:
        print(f"\nbiased pLDDT: mean {sum(pl)/len(pl):.3f}  range {min(pl):.3f}-{max(pl):.3f}")
        print("baseline pLDDT for the SAME backbone: 0.411-0.467")
        print("rfd4_design_045 (folds well): 0.963   HG3.17+MSA: 0.983")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
