#!/usr/bin/env bash
#SBATCH --job-name=kemp-pt5
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
PARENT="$PDH_WORKSPACE/results/kemp-18/bridge/kemp.pdb"
RUN="$PDH_WORKSPACE/results/kemp-pt5-$SLURM_JOB_ID"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
mkdir -p "$RUN/rfd3" "$RUN/bridge" "$RUN/ligandmpnn"

echo "Job: $SLURM_JOB_ID; partial diffusion of $PARENT at 5 A"
cp "$PARENT" "$RUN/bridge/input.pdb"
"$ENV_ROOT/env-rfd3/bin/python" - "$RUN/bridge/input.pdb" "$RUN/partial.json" << 'PY'
import json, sys
json.dump({
    "kemp18_partial": {
        "input": sys.argv[1],
        "ligand": "6NT",
        "partial_t": 5.0,
        "allow_ligand_on_existing_chain": True,
        "select_fixed_atoms": {
            "A41": "OD1,CG,OD2",
            "A65": "OE1,CD,NE2",
            "6NT": "ALL",
        },
    }
}, open(sys.argv[2], "w"), indent=2)
print("wrote", sys.argv[2])
PY

echo "START 1/3: RFD3 partial diffusion"
"$ENV_ROOT/env-rfd3/bin/rfd3" \
    inputs="$RUN/partial.json" \
    'json_keys_subset=[kemp18_partial]' \
    ckpt_path="$RFD3_CHECKPOINT" \
    out_dir="$RUN/rfd3" \
    diffusion_batch_size=1 n_batches=1 seed=1 \
    skip_existing=False prevalidate_inputs=True \
    cleanup_guideposts=True output_full_json=True dump_prediction_metadata_json=True

echo "START 2/3: LigandMPNN, Asp41 and Gln65 fixed"
"$ENV_ROOT/env-rfd3/bin/python" - "$RUN" << 'PY'
import gzip
from pathlib import Path
import sys
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
run = Path(sys.argv[1])
cifs = sorted((run / "rfd3").glob("*.cif.gz"))
if len(cifs) != 1:
    raise SystemExit(f"expected one RFD3 structure, found {cifs}")
with gzip.open(cifs[0], "rt") as handle:
    atoms = pdbx.get_structure(pdbx.CIFFile.read(handle), model=1)
out = run / "bridge" / "kemp.pdb"
file = pdb.PDBFile()
pdb.set_structure(file, atoms)
file.write(out)
seq = {}
for line in out.read_text().splitlines():
    if line.startswith("ATOM") and line[12:16].strip() == "CA":
        seq[int(line[22:26])] = line[17:20].strip()
if seq.get(41) != "ASP" or seq.get(65) != "GLN":
    raise SystemExit(f"catalytic identities moved: 41={seq.get(41)} 65={seq.get(65)}")
print("partial structure", out, "length", len(seq))
PY

"$ENV_ROOT/env-ligandmpnn/bin/python" "$LIGANDMPNN_ROOT/run.py" \
    --model_type ligand_mpnn \
    --checkpoint_ligand_mpnn "$LIGANDMPNN_CHECKPOINT" \
    --pdb_path "$RUN/bridge/kemp.pdb" \
    --out_folder "$RUN/ligandmpnn" \
    --fixed_residues "A41 A65" \
    --ligand_mpnn_use_atom_context 1 --ligand_mpnn_use_side_chain_context 1 \
    --batch_size 1 --number_of_batches 1 --pack_side_chains 0 --seed 3 --verbose 1 \
    2>&1 | tee "$RUN/ligandmpnn.log"

echo "START 3/3: one Boltz-2 model"
"$ENV_ROOT/env-rfd3/bin/python" - "$RUN" "$PROJECT" << 'PY'
import json, os, subprocess, sys
from pathlib import Path
run, project = Path(sys.argv[1]), Path(sys.argv[2])
fastas = sorted((run / "ligandmpnn/seqs").glob("*.fa"))
if len(fastas) != 1:
    raise SystemExit(f"expected one fasta, found {fastas}")
header, chunks = None, []
design = None
for line in fastas[0].read_text().splitlines():
    if line.startswith(">"):
        if header and "id=" in header:
            design = (header, "".join(chunks))
        header, chunks = line[1:], []
    else:
        chunks.append(line.strip())
if header and "id=" in header:
    design = (header, "".join(chunks))
if design is None:
    raise SystemExit("no designed sequence")
header, seq = design
if len(seq) < 65 or seq[40] != "D" or seq[64] != "Q":
    raise SystemExit(f"sequence does not keep D41 and Q65: 41={seq[40:41]} 65={seq[64:65]} len={len(seq)}")
scores = {}
for part in header.split(", "):
    if "=" in part:
        key, value = part.split("=", 1)
        scores[key] = value
bridge = run / "bridge"
(bridge / "sequence.json").write_text(json.dumps({
    "sequence": seq,
    "length": len(seq),
    "fixed": [
        {"original": "G50", "id": "A65", "position": 65, "amino_acid": "Q"},
        {"original": "G127", "id": "A41", "position": 41, "amino_acid": "D"},
    ],
    "ligandmpnn_scores": scores,
}, indent=2) + "\n")
(bridge / "motif.json").write_text(json.dumps({
    "partial_t_angstrom": 5.0,
    "source": "kemp-18/bridge/kemp.pdb",
    "fixed": ["A41", "A65"],
}, indent=2) + "\n")
msa = run / "msa"
msa.mkdir(exist_ok=True)
(msa / "kemp.a3m").write_text(">kemp_design\n" + seq + "\n")
(msa / "verification.json").write_text(json.dumps({
    "kemp_design": {"depth": 1, "nonquery_rows": 0, "query_length": len(seq)}
}, indent=2) + "\n")
yaml = run / "kemp.yaml"
yaml.write_text(
    "version: 1\nsequences:\n- protein:\n    id: A\n    sequence: " + seq
    + "\n    msa: " + str((msa / "kemp.a3m").resolve())
    + "\n- ligand:\n    id: B\n    ccd: 6NT\n"
)
cache = run / "boltz-cache"
cache.mkdir(exist_ok=True)
assets = Path(os.environ["BOLTZ_ASSETS"])
for name in ("boltz2_conf.ckpt", "boltz2_aff.ckpt", "mols.tar", "mols"):
    link = cache / name
    if not link.exists():
        link.symlink_to(assets / name)
log = run / "boltz.log"
with log.open("w") as handle:
    proc = subprocess.run([
        os.environ["ENV_ROOT"] + "/env-boltz2/bin/boltz", "predict", str(yaml),
        "--model", "boltz2", "--cache", str(cache), "--out_dir", str(run / "boltz"),
        "--accelerator", "gpu", "--devices", "1", "--no_kernels",
        "--sampling_steps", "200", "--recycling_steps", "3", "--diffusion_samples", "1",
        "--num_workers", "0", "--preprocessing-threads", "1", "--output_format", "mmcif",
        "--seed", "1",
    ], stdout=handle, stderr=subprocess.STDOUT)
if proc.returncode != 0:
    raise SystemExit(f"Boltz failed, see {log}")
subprocess.run([os.environ["ENV_ROOT"] + "/env-rfd3/bin/python", str(project / "verify_boltz.py"), str(run)], check=True)
PY
