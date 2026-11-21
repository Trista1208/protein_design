#!/usr/bin/env bash
#SBATCH --job-name=kemp-seq
#SBATCH --partition=gpu
#SBATCH --constraint=g6.xlarge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=01:30:00

set -e -o pipefail

source "${PDH_WORKSPACE:?Submit from the project directory}/env.sh"
export PDH_RUN_MODE=project
PROJECT="$PDH_WORKSPACE/src"
PARENT="$PDH_WORKSPACE/results/kemp-18"
RUN_ROOT="$PDH_WORKSPACE/results/kemp-seq-$SLURM_JOB_ID"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
mkdir -p "$RUN_ROOT/ligandmpnn"

echo "Job: $SLURM_JOB_ID; node: $(hostname); backbone: $PARENT/bridge/kemp.pdb"
cp "$PARENT/bridge/kemp.pdb" "$RUN_ROOT/backbone.pdb"
cp "$PARENT/bridge/motif.json" "$RUN_ROOT/motif.json"
cp "$PARENT/bridge/fixed_residues.txt" "$RUN_ROOT/fixed_residues.txt"

echo "START 1/2: LigandMPNN, 4 sequences, Asp41 and Gln65 fixed"
"$ENV_ROOT/env-ligandmpnn/bin/python" "$LIGANDMPNN_ROOT/run.py" \
    --model_type ligand_mpnn \
    --checkpoint_ligand_mpnn "$LIGANDMPNN_CHECKPOINT" \
    --pdb_path "$RUN_ROOT/backbone.pdb" \
    --out_folder "$RUN_ROOT/ligandmpnn" \
    --fixed_residues "$(cat "$RUN_ROOT/fixed_residues.txt")" \
    --ligand_mpnn_use_atom_context 1 --ligand_mpnn_use_side_chain_context 1 \
    --batch_size 1 --number_of_batches 4 --pack_side_chains 0 --seed 2 --verbose 1 \
    2>&1 | tee "$RUN_ROOT/ligandmpnn.log"

"$ENV_ROOT/env-rfd3/bin/python" - "$RUN_ROOT" << 'PY'
import json
import sys
from pathlib import Path

run = Path(sys.argv[1])
fastas = sorted((run / "ligandmpnn/seqs").glob("*.fa"))
if len(fastas) != 1:
    raise SystemExit(f"expected one LigandMPNN fasta, found {fastas}")
text = fastas[0].read_text().splitlines()
designs = []
header, chunks = None, []

def flush():
    if header and "id=" in header:
        designs.append((header, "".join(chunks)))

for line in text:
    if line.startswith(">"):
        flush()
        header, chunks = line[1:], []
    else:
        chunks.append(line.strip())
flush()
kept = []
for header, seq in designs:
    if len(seq) != 171 or seq[40] != "D" or seq[64] != "Q":
        print("skip", header[:80], "len", len(seq))
        continue
    fields = {}
    for part in header.split(", "):
        if "=" in part:
            key, value = part.split("=", 1)
            fields[key] = value
    kept.append({"header": header, "sequence": seq, "scores": fields})
if not kept:
    raise SystemExit("no usable LigandMPNN sequences")
(run / "designs.json").write_text(json.dumps(kept, indent=2) + "\n")
print(f"kept {len(kept)} sequences")
PY

echo "START 2/2: one Boltz-2 model per sequence"
"$ENV_ROOT/env-rfd3/bin/python" - "$RUN_ROOT" "$PROJECT" << 'PY'
import json
import os
import subprocess
import sys
from pathlib import Path

run = Path(sys.argv[1])
project = Path(sys.argv[2])
designs = json.loads((run / "designs.json").read_text())
assets = Path(os.environ["BOLTZ_ASSETS"])
boltz = os.environ["ENV_ROOT"] + "/env-boltz2/bin/boltz"
python = os.environ["ENV_ROOT"] + "/env-rfd3/bin/python"
for index, design in enumerate(designs, start=1):
    seq_run = run / f"seq{index}"
    bridge = seq_run / "bridge"
    bridge.mkdir(parents=True, exist_ok=True)
    (bridge / "kemp.pdb").write_bytes((run / "backbone.pdb").read_bytes())
    (bridge / "motif.json").write_bytes((run / "motif.json").read_bytes())
    (bridge / "sequence.json").write_text(json.dumps({
        "sequence": design["sequence"],
        "length": len(design["sequence"]),
        "fixed": [
            {"original": "G50", "id": "A65", "position": 65, "amino_acid": "Q"},
            {"original": "G127", "id": "A41", "position": 41, "amino_acid": "D"},
        ],
        "ligandmpnn_scores": design["scores"],
    }, indent=2) + "\n")
    msa = seq_run / "msa"
    msa.mkdir(parents=True, exist_ok=True)
    (msa / "kemp.a3m").write_text(">kemp_design\n" + design["sequence"] + "\n")
    yaml = seq_run / "kemp.yaml"
    yaml.write_text(
        "version: 1\nsequences:\n"
        "- protein:\n    id: A\n    sequence: " + design["sequence"] + "\n"
        "    msa: " + str((msa / "kemp.a3m").resolve()) + "\n"
        "- ligand:\n    id: B\n    ccd: 6NT\n"
    )
    cache = seq_run / "boltz-cache"
    cache.mkdir(parents=True, exist_ok=True)
    for name in ("boltz2_conf.ckpt", "boltz2_aff.ckpt", "mols.tar", "mols"):
        link = cache / name
        if not link.exists():
            link.symlink_to(assets / name)
    log = seq_run / "boltz.log"
    print(f"Boltz sequence {index}", flush=True)
    with log.open("w") as handle:
        proc = subprocess.run(
            [
                boltz, "predict", str(yaml),
                "--model", "boltz2", "--cache", str(cache), "--out_dir", str(seq_run / "boltz"),
                "--accelerator", "gpu", "--devices", "1", "--no_kernels",
                "--sampling_steps", "200", "--recycling_steps", "3", "--diffusion_samples", "1",
                "--num_workers", "0", "--preprocessing-threads", "1", "--output_format", "mmcif",
                "--seed", str(index),
            ],
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    if proc.returncode != 0:
        print(f"Boltz sequence {index} failed, see {log}")
        continue
    subprocess.run([python, str(project / "verify_boltz.py"), str(seq_run)], check=False)
print("sequence sweep done")
PY
