#!/bin/bash
# Per-design cost at PRODUCTION settings, so 5/10/50/100/500 can be costed honestly.
#
# Measures the two dominant stages:
#   RFD3  : length 190-250, 200 timesteps, batch 1   (the real generator setting)
#   Boltz-2: ~240 aa + 6NT ligand                    (the real scoring setting)
# LigandMPNN was already measured at 0.55 GB / 4 s for 8 sequences.
#
# RFD3 is run TWICE: the first run JIT-compiles Triton kernels and is not
# representative. Only run 2 is reported.
set -u
E=/home/mandalorian/anaconda3/envs
D=/home/mandalorian/Hackathons/Enzymes-BSSE2026
B=$D/outputs/bench; rm -rf $B; mkdir -p $B
CLEAN=$D/outputs/runs/m0/theozyme_clean.pdb
RESULTS=$B/results.tsv
printf "stage\trun\twall_s\tpeak_mib\tdelta_mib\tnote\n" > $RESULTS

sample() { # sample <tag> <note> -- <cmd...>
  local tag="$1" note="$2"; shift 3
  ( while true; do
      nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits 2>/dev/null >> "$B/$tag.mem"
      sleep 0.5
    done ) & local SP=$!
  local t0=$(date +%s)
  "$@" > "$B/$tag.log" 2>&1; local rc=$?
  local t1=$(date +%s)
  kill $SP 2>/dev/null; wait $SP 2>/dev/null
  local peak=$(sort -n "$B/$tag.mem" 2>/dev/null | tail -1)
  local base=$(sort -n "$B/$tag.mem" 2>/dev/null | head -1)
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$tag" "$rc" "$((t1-t0))" "${peak:-0}" \
     "$(( ${peak:-0} - ${base:-0} ))" "$note" >> $RESULTS
  printf "  %-22s rc=%s  wall=%4ss  peak=%6s MiB  delta=%6s MiB  %s\n" \
     "$tag" "$rc" "$((t1-t0))" "${peak:-?}" "$(( ${peak:-0} - ${base:-0} ))" "$note"
}

echo "=== idle: $(nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader) ==="
echo

# ---------------------------------------------------------------- RFD3 production
cat > $B/rfd3_prod.py <<PY
import sys, json, pathlib
sys.path.insert(0, "$D")
from kemp import rfd3
from kemp.theozyme import load
tz = load("$CLEAN")
spec = rfd3.build_spec("canonical", tz, length="190-250")
d = spec.to_dict(); d.pop("ligand", None)      # ligand already present in input
p = pathlib.Path("$B/rfd3_inputs.json"); p.write_text(json.dumps({"canonical": d}, indent=2))
out = sys.argv[1]
r = rfd3.run(p, out, batch_size=1, n_batches=1, timesteps=200)
print("rc", r["returncode"], "n_structures", r["n_structures"])
print("structures", r["structures"][:4])
PY
echo "--- RFD3 @ length 190-250, 200 timesteps, batch 1"
sample rfd3_run1 "FIRST run - includes Triton JIT, NOT representative" -- \
   $E/foundry/bin/python $B/rfd3_prod.py $B/rfd3_out1
sample rfd3_run2 "second run - THIS is the per-design cost" -- \
   $E/foundry/bin/python $B/rfd3_prod.py $B/rfd3_out2

# batch of 4 to see whether batching amortises
cat > $B/rfd3_batch.py <<PY
import sys, json, pathlib
sys.path.insert(0, "$D")
from kemp import rfd3
from kemp.theozyme import load
tz = load("$CLEAN")
spec = rfd3.build_spec("canonical", tz, length="190-250")
d = spec.to_dict(); d.pop("ligand", None)
p = pathlib.Path("$B/rfd3_inputs_b.json"); p.write_text(json.dumps({"canonical": d}, indent=2))
r = rfd3.run(p, "$B/rfd3_out_b4", batch_size=4, n_batches=1, timesteps=200)
print("rc", r["returncode"], "n_structures", r["n_structures"])
PY
sample rfd3_batch4 "batch_size=4 - all 4 share one length" -- \
   $E/foundry/bin/python $B/rfd3_batch.py

# ---------------------------------------------------------------- Boltz-2 production
SEQ240=$($E/mm/bin/python -c "
import random; random.seed(1)
aa='ACDEFGHIKLMNPQRSTVWY'
print('M'+''.join(random.choice(aa) for _ in range(239)))")
printf '>A|protein\n%s\n>B|ccd\n6NT\n' "$SEQ240" > $B/boltz240.fasta
echo
echo "--- Boltz-2 @ 240 aa + 6NT ligand, 1 seed"
sample boltz_240 "240 aa + ligand, recycling=3, 1 sample" -- \
   $E/her2_boltz2/bin/boltz predict $B/boltz240.fasta --out_dir $B/boltz_out \
     --diffusion_samples 1 --recycling_steps 3 --use_msa_server --override

echo
echo "=== RESULTS ==="
column -t $RESULTS
echo
echo "### BENCH DONE ###"
