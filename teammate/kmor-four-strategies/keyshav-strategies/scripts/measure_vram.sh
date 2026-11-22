#!/bin/bash
# Measure PEAK VRAM + wall clock per GPU tool, one at a time, for docs/CONCURRENCY.md.
set -u
E=/home/mandalorian/anaconda3/envs
D=/home/mandalorian/Hackathons/Enzymes-BSSE2026
V=$D/outputs/vram; rm -rf $V; mkdir -p $V
SEQ="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQKLGSGSYATWFKQALDNA"   # 50 aa

sample() { # sample <tag> -- <cmd...>
  local tag="$1"; shift; shift
  ( while true; do
      nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null >> "$V/$tag.mem"
      sleep 0.4
    done ) & local SP=$!
  local t0=$(date +%s)
  "$@" > "$V/$tag.log" 2>&1; local rc=$?
  local t1=$(date +%s)
  kill $SP 2>/dev/null; wait $SP 2>/dev/null
  local peak=$(sort -n "$V/$tag.mem" 2>/dev/null | tail -1)
  local base=$(sort -n "$V/$tag.mem" 2>/dev/null | head -1)
  printf "%-14s rc=%-3s peak=%-6s MiB  baseline=%-6s MiB  delta=%-6s MiB  wall=%ss\n" \
     "$tag" "$rc" "${peak:-?}" "${base:-?}" "$(( ${peak:-0} - ${base:-0} ))" "$((t1-t0))" \
     | tee -a "$V/summary.txt"
}

echo "idle GPU: $(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader)" | tee "$V/summary.txt"
echo | tee -a "$V/summary.txt"

# --- ESMFold
cat > $V/esmfold.py <<PY
import torch, sys
from transformers import AutoTokenizer, EsmForProteinFolding
m = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True).cuda().eval()
m.esm = m.esm.half()
with torch.no_grad():
    pdb = m.infer_pdb("$SEQ")
print("atoms", pdb.count("ATOM"))
print("torch peak MiB", torch.cuda.max_memory_allocated()//1048576)
PY
sample esmfold -- $E/aizymes/bin/python $V/esmfold.py

# --- LigandMPNN
sample ligandmpnn -- $E/aizymes/bin/python $D/tools/LigandMPNN/run.py \
   --model_type ligand_mpnn \
   --checkpoint_ligand_mpnn $D/tools/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt \
   --pdb_path $D/tools/LigandMPNN/inputs/1BC8.pdb \
   --out_folder $V/lmpnn_out --batch_size 8 --number_of_batches 1 --seed 37

# --- Chai-1
cat > $V/chai.py <<PY
from pathlib import Path
import torch
from chai_lab.chai1 import run_inference
Path("$V/chai_in.fasta").write_text(">protein|x\n$SEQ\n")
run_inference(fasta_file=Path("$V/chai_in.fasta"), output_dir=Path("$V/chai_out"),
              num_trunk_recycles=1, num_diffn_timesteps=20, seed=1,
              device="cuda:0", use_esm_embeddings=True)
print("torch peak MiB", torch.cuda.max_memory_allocated()//1048576)
PY
sample chai1 -- $E/chai1/bin/python $V/chai.py

# --- Boltz-2
printf '>A|protein\n%s\n' "$SEQ" > $V/boltz_in.fasta
sample boltz2 -- $E/her2_boltz2/bin/boltz predict $V/boltz_in.fasta \
   --out_dir $V/boltz_out --diffusion_samples 1 --recycling_steps 1 \
   --use_msa_server --override

# --- RF3 (needs non-fasta input: give it the clean theozyme PDB)
CLEAN=$D/outputs/runs/m0/theozyme_clean.pdb
sample rf3 -- $E/foundry/bin/rf3 fold inputs=$CLEAN out_dir=$V/rf3_out seed=1

# --- RFD3 via the project's own spec builder (correct schema)
cat > $V/rfd3_run.py <<PY
import sys, json
sys.path.insert(0, "$D")
from kemp import rfd3, config as C
from kemp.theozyme import load
tz = load("$CLEAN")
spec = rfd3.build_spec("canonical", tz, length="190-200")
js = rfd3.write_inputs_json([spec], "$V/rfd3_inputs.json")
print(json.dumps(json.loads(open(js).read()), indent=2))
res = rfd3.run(js, "$V/rfd3_out", batch_size=1, n_batches=1, timesteps=20)
print(json.dumps(res, indent=2, default=str))
PY
sample rfd3 -- $E/foundry/bin/python $V/rfd3_run.py

echo | tee -a "$V/summary.txt"
echo "=== SUMMARY ===" | tee -a "$V/summary.txt"
cat "$V/summary.txt"
echo "### VRAM DONE ###"
