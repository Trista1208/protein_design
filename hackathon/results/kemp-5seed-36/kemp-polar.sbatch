#!/usr/bin/env bash
#SBATCH --job-name=kemp-polar
#SBATCH --partition=gpu
#SBATCH --constraint=g6.xlarge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=00:40:00

# Polar-pocket design on the passing Boltz model.
# No cysteine. The leucine cage around Gln65 is biased toward Ser/Thr/Asn.
# Lys/Arg next to Asp41 are biased down so they do not capture the base.

set -e -o pipefail
source "${PDH_WORKSPACE:?Submit from the project directory}/env.sh"
export PDH_RUN_MODE=project
SRC="$PDH_WORKSPACE/results/kemp-5seed-36/boltz/boltz_results_kemp/predictions/kemp/kemp_model_2.cif"
RUN="$PDH_WORKSPACE/results/kemp-polar-$SLURM_JOB_ID"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
mkdir -p "$RUN/bridge" "$RUN/ligandmpnn"

echo "Job: $SLURM_JOB_ID; polar pocket, no cysteine, from passing model 2"
"$ENV_ROOT/env-rfd3/bin/python" - "$SRC" "$RUN" << 'PY'
import json, sys
from pathlib import Path
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
src, run = Path(sys.argv[1]), Path(sys.argv[2])
atoms = pdbx.get_structure(pdbx.CIFFile.read(str(src)), model=1, use_author_fields=True)
out = run / "bridge" / "kemp.pdb"
file = pdb.PDBFile()
pdb.set_structure(file, atoms)
file.write(out)
bias = {}
polar = {"S": 1.5, "T": 1.2, "N": 1.2, "Q": 0.6, "L": -1.5, "I": -1.5, "V": -1.5, "F": -1.5, "C": -5.0}
calm = {"K": -2.0, "R": -2.0, "S": 1.2, "T": 1.0, "N": 1.0, "C": -5.0}
for residue in (17, 51, 52, 53, 63, 64, 66, 67, 76, 77, 78):
    bias[f"A{residue}"] = polar
for residue in (47, 49, 50, 163):
    bias[f"A{residue}"] = calm
(run / "shell_bias.json").write_text(json.dumps(bias, indent=2) + "\n")
print("pdb", out, "biased", len(bias))
PY

echo "START 1/2: four polar sequences, Asp41 and Gln65 fixed, cysteine omitted"
"$ENV_ROOT/env-ligandmpnn/bin/python" "$LIGANDMPNN_ROOT/run.py" \
    --model_type ligand_mpnn \
    --checkpoint_ligand_mpnn "$LIGANDMPNN_CHECKPOINT" \
    --pdb_path "$RUN/bridge/kemp.pdb" \
    --out_folder "$RUN/ligandmpnn" \
    --fixed_residues "A41 A65" \
    --omit_AA C \
    --bias_AA_per_residue "$RUN/shell_bias.json" \
    --ligand_mpnn_use_atom_context 1 --ligand_mpnn_use_side_chain_context 1 \
    --batch_size 1 --number_of_batches 4 --temperature 0.3 \
    --pack_side_chains 0 --seed 7 --verbose 1 \
    2>&1 | tee "$RUN/ligandmpnn.log"

echo "START 2/2: one Boltz model per sequence"
"$ENV_ROOT/env-rfd3/bin/python" - "$RUN" << 'PY'
import json, os, subprocess, sys
from pathlib import Path
import numpy as np
from biotite.structure.io import pdbx

run = Path(sys.argv[1])
fastas = sorted((run / "ligandmpnn/seqs").glob("*.fa"))
if len(fastas) != 1:
    raise SystemExit(f"expected one fasta, found {fastas}")
designs, header, chunks = [], None, []

def flush():
    if header and "id=" in header:
        designs.append((header, "".join(chunks)))

for line in fastas[0].read_text().splitlines():
    if line.startswith(">"):
        flush()
        header, chunks = line[1:], []
    else:
        chunks.append(line.strip())
flush()
kept = [(h, s) for h, s in designs if len(s) == 171 and s[40] == "D" and s[64] == "Q" and "C" not in s]
print(f"kept {len(kept)} cysteine-free sequences")
if not kept:
    raise SystemExit("no cysteine-free sequence kept D41 and Q65")
(run / "sequences.json").write_text(json.dumps([{"header": h, "sequence": s} for h, s in kept], indent=2) + "\n")

def contacts(cif):
    atoms = pdbx.get_structure(pdbx.CIFFile.read(str(cif)), model=1, use_author_fields=True)
    protein = atoms[(~atoms.hetero) & (atoms.chain_id == "A")]
    ligand = atoms[(atoms.res_name == "6NT") & (atoms.element != "H")]
    residue_ids = []
    for residue in protein[protein.atom_name == "CA"].res_id:
        if int(residue) not in residue_ids:
            residue_ids.append(int(residue))
    def point(index, name):
        selected = protein[(protein.res_id == residue_ids[index]) & (protein.atom_name == name)]
        return selected.coord[0]
    def lig(name):
        return ligand[ligand.atom_name == name].coord[0]
    base = float(np.linalg.norm(point(40, "OD2") - lig("N3")))
    donor = float(np.linalg.norm(point(64, "NE2") - lig("N1")))
    return base, donor

assets = Path(os.environ["BOLTZ_ASSETS"])
rows = []
for index, (_, seq) in enumerate(kept, start=1):
    seq_run = run / f"seq{index}"
    msa = seq_run / "msa"
    msa.mkdir(parents=True)
    (msa / "kemp.a3m").write_text(">kemp_design\n" + seq + "\n")
    yaml = seq_run / "kemp.yaml"
    yaml.write_text(
        "version: 1\nsequences:\n- protein:\n    id: A\n    sequence: " + seq
        + "\n    msa: " + str((msa / "kemp.a3m").resolve())
        + "\n- ligand:\n    id: B\n    ccd: 6NT\n"
    )
    cache = seq_run / "boltz-cache"
    cache.mkdir()
    for name in ("boltz2_conf.ckpt", "boltz2_aff.ckpt", "mols.tar", "mols"):
        (cache / name).symlink_to(assets / name)
    log = seq_run / "boltz.log"
    print(f"Boltz sequence {index}", flush=True)
    with log.open("w") as handle:
        proc = subprocess.run([
            os.environ["ENV_ROOT"] + "/env-boltz2/bin/boltz", "predict", str(yaml),
            "--model", "boltz2", "--cache", str(cache), "--out_dir", str(seq_run / "boltz"),
            "--accelerator", "gpu", "--devices", "1", "--no_kernels",
            "--sampling_steps", "200", "--recycling_steps", "3", "--diffusion_samples", "1",
            "--num_workers", "0", "--preprocessing-threads", "1", "--output_format", "mmcif",
            "--seed", "1",
        ], stdout=handle, stderr=subprocess.STDOUT)
    cif = seq_run / "boltz/boltz_results_kemp/predictions/kemp/kemp_model_0.cif"
    if proc.returncode != 0 or not cif.is_file():
        rows.append({"sequence": index, "error": "boltz failed"})
        continue
    base, donor = contacts(cif)
    rows.append({
        "sequence": index,
        "sequence_string": seq,
        "ASP_OD2_to_6NT_N3": base,
        "GLN_NE2_to_6NT_N1": donor,
        "under_5": base < 5 and donor < 5,
        "cysteines": seq.count("C"),
    })
    print(f"seq {index}: base {base:.2f} donor {donor:.2f} cys {seq.count('C')}")
(run / "contacts.json").write_text(json.dumps(rows, indent=2) + "\n")
print("wrote", run / "contacts.json")
PY
