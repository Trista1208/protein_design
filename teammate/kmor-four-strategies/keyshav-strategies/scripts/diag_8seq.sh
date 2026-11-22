#!/bin/bash
# Which axis is actually failing?
#
# The sweep samples BACKBONES (100 RFD3 seeds) but folds only 1 of the 8
# LigandMPNN sequences per backbone, with 1 Boltz seed. If the active site dies
# after a good backbone, more backbones cannot fix it.
#
# Take ONE backbone that RFD3 built well and fold ALL 8 sequences:
#   some reproduce the site  -> failure is SEQUENCE sampling. Fold more
#                               sequences per backbone, not more backbones.
#   all 8 fail               -> failure is Boltz itself: the site is real but
#                               not recoverable from sequence alone. No amount
#                               of sampling on any axis helps; needs a loop.
#
# ~11 min (8 x 80 s) and worth more than the next 50 designs.
set -u
D=/home/mandalorian/Hackathons/Enzymes-BSSE2026
PY=/home/mandalorian/anaconda3/envs/mm/bin/python
OUT=$D/outputs/diag_8seq_fixed
rm -rf "$OUT"; mkdir -p "$OUT"
cd "$D" || exit 1
export KEMP_ROOT="$D"

"$PY" - <<'PYEOF' 2>&1
import json, sys, pathlib
sys.path.insert(0, "/home/mandalorian/Hackathons/Enzymes-BSSE2026")
from kemp import apnovo, config as C, official, pipeline, rfd3
from kemp.runners import ligandmpnn
from kemp.theozyme import load, write_clean_pdb

OUT = pathlib.Path("/home/mandalorian/Hackathons/Enzymes-BSSE2026/outputs/diag_8seq_fixed")

# 1. one backbone, canonical, production settings
tz = load(C.THEOZYME_DESIGN)
clean = write_clean_pdb(tz, OUT / "theozyme_clean.pdb")
tzc = load(clean)
spec = rfd3.build_spec("canonical", tzc, length="190-210")
d = spec.to_dict()   # KEEP the ligand: RFD3 must hold it
js = OUT / "inputs.json"; js.write_text(json.dumps({"canonical": d}, indent=2))
res = rfd3.run(js, OUT / "rfd3", batch_size=1, n_batches=1, timesteps=200, seed=909)
print(f"  RFD3 rc={res['returncode']} structures={res['n_structures']}", flush=True)
if res["returncode"] != 0 or res["n_structures"] == 0:
    sys.exit("  RFD3 failed - cannot run the diagnostic")

site = pipeline.catalytic_site_from_rfd3(OUT / "rfd3", base_tag="G2", oxy_tag="G1")
if site is None:
    # clean pdb renumbers GLN->1, ASP->2
    site = pipeline.catalytic_site_from_rfd3(OUT / "rfd3", base_tag="G1", oxy_tag="G2")
print(f"  catalytic site: base={site.base_resseq} oxy={site.oxy_resseq} "
      f"chain={site.chain}", flush=True)

cifs = sorted((OUT / "rfd3").rglob("*.cif.gz")) or sorted((OUT / "rfd3").rglob("*.cif"))
bb = pipeline.cif_to_pdb(cifs[0], OUT / "backbone.pdb")

m = json.loads(sorted((OUT / "rfd3").rglob("*_model_*.json"))[0].read_text())["metrics"]
print(f"  backbone quality: insertion_rmsd={m.get('insertion_rmsd'):.3f} "
      f"chainbreaks={m.get('n_chainbreaks')} clashes="
      f"{(m.get('n_clashing.interresidue_clashes_w_backbone') or 0)}", flush=True)

# 2. ALL 8 sequences from that one backbone
fixed = f"{site.chain}{site.base_resseq} {site.chain}{site.oxy_resseq}"
mp = ligandmpnn(bb, OUT / "lmpnn", n_seqs=8, fixed_residues=fixed)
seqs = mp.sequences[1:]
print(f"  LigandMPNN produced {len(seqs)} sequences", flush=True)

pocket = [site.base_resseq, site.oxy_resseq]
rows = []
for i, (hdr, seq) in enumerate(seqs):
    rec = pipeline.DesignRecord(design_id=f"s{i}", variant="canonical",
                                sequence=seq, length=len(seq))
    rec = pipeline.score_design(rec, bb, replicates=1,
                                out_dir=OUT / f"boltz_s{i}", site=site)
    off = rec.official or {}
    sub = off.get("subscores") or {}
    det = off.get("detail") or {}
    row = {"seq": i, "len": len(seq),
           "d_base": (det.get("base") or {}).get("d_base_A"),
           "d_oxy": (det.get("oxy") or {}).get("d_oxy_A"),
           "base": sub.get("base"), "oxy": sub.get("oxy"),
           "ipde": sub.get("ipde"), "bury": sub.get("bury"),
           "S": off.get("official_score") or off.get("prefilter_score")}
    rows.append(row)
    print(f"  seq {i}: d_base={row['d_base']} d_oxy={row['d_oxy']} "
          f"base={row['base']} oxy={row['oxy']} S={row['S']}", flush=True)

(OUT / "results.json").write_text(json.dumps(rows, indent=2, default=str))

nz_base = [r for r in rows if (r["base"] or 0) > 0]
nz_oxy = [r for r in rows if (r["oxy"] or 0) > 0]
print("\n  ================= VERDICT =================")
print(f"  sequences with non-zero base: {len(nz_base)}/8")
print(f"  sequences with non-zero oxy : {len(nz_oxy)}/8")
if nz_oxy:
    print("  => SEQUENCE sampling matters. Fold more sequences per backbone,")
    print("     not more backbones. The sweep is sampling the wrong axis.")
elif nz_base:
    print("  => the base can be reproduced but the oxyanion never is.")
    print("     The donor placement is the systematic failure.")
else:
    print("  => NO sequence from a well-built backbone reproduces the site.")
    print("     Boltz cannot recover the active site from sequence alone.")
    print("     More sampling on ANY axis will not help - needs a feedback loop.")
PYEOF
echo "### DIAG DONE ###"
