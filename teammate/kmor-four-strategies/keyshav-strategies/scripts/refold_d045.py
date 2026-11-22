"""5-seed Boltz-2 refold of rfd4_design_045, the team's best design.

WHY: its only zeroed term is orient_oxy, because the oxyanion dihedral sits at
86.1 deg, in the dead zone between the 130-230 and 310-50 bands. A rotamer scan
showed 360 clash-free Gln rotamers that satisfy BOTH the distance window and the
dihedral band, the best reaching orient_oxy 0.745 and a 5-term geometric mean of
0.721. So the pocket is geometrically capable.

But that scan edited coordinates, not sequence. If the organisers fold the
SEQUENCE (handover: "Replicates 5 per sequence"), what matters is which rotamer
a predictor picks. This run asks exactly that, and supplies the two terms the
OpenFold3 structure cannot: complex_ipde, and sd_base for cons.

PROTOCOL: the project's pinned settings - diffusion_samples 1, recycling 3,
msa: empty, 5 seeds. NO pocket constraint: control_fold_matrix.py showed it is
worthless with an MSA (0.775 -> 0.756) and cannot rescue the oxyanion without
one, so including it would only deviate from the official protocol for nothing.

CAVEAT: the design was produced and folded with OpenFold3. Boltz-2 is a
different predictor, so a disagreement here is predictor difference as much as
pose instability. Read the spread, not any single seed.
"""
from __future__ import annotations
import json, subprocess, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C, official, pipeline, strain
from kemp.official import LAUKO

SEQ = ("SSFLELMHKAVSAAKQGDFTSAISDFVEACNLLLPDDKKDDFEAMQKALLSTLLAHHDQDRAKDALTTLLTVTLILFE"
       "KHHDDVVELMQKHYDHFHDMVSDDHKIVTEYLKGELDYEDYLTSLVHNRGDALNLSEAEITEIYEEAKKMLS")
BASE_IDX, OXY_IDX = 102, 46
OUT = Path(__file__).resolve().parents[1] / "outputs" / "refold_d045"

def inband(v):
    for lo, hi in (LAUKO["dihedral_sp2_180"], LAUKO["dihedral_sp2_0"]):
        if (lo <= v <= hi) if lo <= hi else (v >= lo or v <= hi):
            return True
    return False

def main() -> int:
    assert len(SEQ) == 150 and SEQ[BASE_IDX-1] == "D" and SEQ[OXY_IDX-1] == "Q"
    print(f"rfd4_design_045: {len(SEQ)} aa, Asp{BASE_IDX} / Gln{OXY_IDX}")
    print("OpenFold3 reference: d_base 2.858  d_oxy 3.094  dihedral_oxy 86.1 (OUT of band)\n")
    raws, ipdes, rows = [], [], []
    for seed in C.BOLTZ_SEEDS[:C.REPLICATES]:
        w = OUT / f"seed{seed}"; w.mkdir(parents=True, exist_ok=True)
        (w/"input.yaml").write_text("\n".join([
            "version: 1","sequences:","  - protein:","      id: A",
            f'      sequence: "{SEQ}"',"      msa: empty",
            "  - ligand:","      id: B",f"      ccd: {C.LIGAND_CCD}"])+"\n")
        t0=time.time()
        with (w/"boltz.log").open("w") as fh:
            p=subprocess.run([str(C.BIN_BOLTZ),"predict",str(w/"input.yaml"),"--out_dir",str(w),
                "--diffusion_samples","1","--recycling_steps","3","--seed",str(seed),"--override"],
                stdout=fh, stderr=subprocess.STDOUT, timeout=5400)
        cifs=sorted(w.rglob("*_model_0.cif"))
        if p.returncode!=0 or not cifs:
            print(f"  seed {seed}: FAILED rc={p.returncode}"); continue
        pdb=pipeline.cif_to_pdb(cifs[0], w/"folded.pdb")
        raw=official.measure_raw(str(pdb), site=official.CatalyticSite(
            base_resseq=BASE_IDX, oxy_resseq=OXY_IDX, chain="A"))
        raws.append(raw)
        conf=None
        for cj in sorted(w.rglob("confidence*.json")):
            conf=json.loads(cj.read_text()); break
        if conf and conf.get("complex_ipde") is not None: ipdes.append(conf["complex_ipde"])
        b = inband(raw.dihedral12_base) if raw.dihedral12_base is not None else False
        o = inband(raw.dihedral12_oxy) if raw.dihedral12_oxy is not None else False
        rows.append(dict(seed=seed, d_base=raw.d_base, d_oxy=raw.d_oxy,
            dih_b=raw.dihedral12_base, dih_o=raw.dihedral12_oxy,
            band=("B" if b else "-")+("O" if o else "-"),
            plddt=(conf or {}).get("complex_plddt"), ipde=(conf or {}).get("complex_ipde"),
            wall=round(time.time()-t0,1)))
        r=rows[-1]
        print(f"  seed {seed}: d_base {r['d_base']:.3f}  d_oxy {r['d_oxy']:.3f}  "
              f"dih_b {r['dih_b']:.1f}  dih_o {r['dih_o']:.1f}  [{r['band']}]  "
              f"pLDDT {r['plddt']:.3f}  iPDE {r['ipde']:.2f}  ({r['wall']}s)", flush=True)
        OUT.mkdir(parents=True,exist_ok=True)
        (OUT/"replicates.json").write_text(json.dumps(rows,indent=2,default=str))
    if not raws: return 3
    comb=official.combine_replicates(raws, ipde_values=ipdes)
    sc=official.score(comb, strict=True)
    print("\n" + "="*68)
    print(f"AGGREGATE over {comb.n_replicates} replicates")
    print(f"  d_base {comb.d_base:.3f}  d_oxy {comb.d_oxy:.3f}  sd_base {comb.sd_base:.4f}  ipde {comb.ipde}")
    print("  subscores: " + " ".join(f"{k}={round(v,3) if v is not None else None}"
                                     for k,v in sc['subscores'].items()))
    print(f"  official_score = {sc.get('official_score')}   limiting = {sc.get('limiting_term')}")
    if sc.get("verdict","").startswith("INCOMPLETE"): print(f"  {sc['verdict'][:120]}")
    n=sum(1 for r in rows if r["band"]=="BO")
    print(f"\n  seeds with BOTH dihedrals in band: {n}/{len(rows)}")
    print("  -> " + ("Boltz DOES sample the in-band rotamer; the fix is reachable without a design change."
                     if n else "no replicate found the in-band rotamer; a design change is needed."))
    (OUT/"verdict.json").write_text(json.dumps({"replicates":rows,
        "aggregate":{k:(round(v,4) if isinstance(v,float) else v) for k,v in sc['subscores'].items()},
        "official_score":sc.get("official_score"),"sd_base":comb.sd_base,
        "n_in_band":n},indent=2,default=str))
    return 0
if __name__=="__main__": raise SystemExit(main())
