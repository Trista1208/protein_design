#!/usr/bin/env bash
#SBATCH --job-name=kemp-des
#SBATCH --partition=gpu
#SBATCH --constraint=g6.xlarge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=00:40:00

set -e -o pipefail
source "${PDH_WORKSPACE:?Submit from the project directory}/env.sh"
export PDH_RUN_MODE=project
SRC="$PDH_WORKSPACE/results/kemp-5seed-36/boltz/boltz_results_kemp/predictions/kemp/kemp_model_2.cif"
RUN="$PDH_WORKSPACE/results/kemp-des-$SLURM_JOB_ID"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
mkdir -p "$RUN/bridge" "$RUN/ligandmpnn"

echo "Job: $SLURM_JOB_ID; design from passing Boltz model 2"
"$ENV_ROOT/env-rfd3/bin/python" - "$SRC" "$RUN/bridge/kemp.pdb" << 'PY'
import sys
from pathlib import Path
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
src, out = Path(sys.argv[1]), Path(sys.argv[2])
atoms = pdbx.get_structure(pdbx.CIFFile.read(str(src)), model=1, use_author_fields=True)
file = pdb.PDBFile()
pdb.set_structure(file, atoms)
out.parent.mkdir(parents=True, exist_ok=True)
file.write(out)
seq = {}
ligand = set()
for line in out.read_text().splitlines():
    if line.startswith("ATOM") and line[12:16].strip() == "CA":
        seq[int(line[22:26])] = line[17:20].strip()
    if line.startswith("HETATM") and line[17:20].strip() == "6NT":
        ligand.add(line[21] + line[22:26].strip())
if seq.get(41) != "ASP" or seq.get(65) != "GLN":
    raise SystemExit(f"passing model lost catalytic identities: 41={seq.get(41)} 65={seq.get(65)}")
print("pdb", out, "length", len(seq), "6NT", sorted(ligand))
PY

echo "START 1/2: four sequences, Asp41 and Gln65 fixed"
"$ENV_ROOT/env-ligandmpnn/bin/python" "$LIGANDMPNN_ROOT/run.py" \
    --model_type ligand_mpnn \
    --checkpoint_ligand_mpnn "$LIGANDMPNN_CHECKPOINT" \
    --pdb_path "$RUN/bridge/kemp.pdb" \
    --out_folder "$RUN/ligandmpnn" \
    --fixed_residues "A41 A65" \
    --ligand_mpnn_use_atom_context 1 --ligand_mpnn_use_side_chain_context 1 \
    --batch_size 1 --number_of_batches 4 --pack_side_chains 0 --seed 4 --verbose 1 \
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
designs = []
header, chunks = None, []

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
kept = []
for header, seq in designs:
    if len(seq) != 171 or seq[40] != "D" or seq[64] != "Q":
        print("skip", len(seq), seq[40:41], seq[64:65])
        continue
    kept.append((header, seq))
if not kept:
    raise SystemExit("no sequences kept D41 and Q65")
print(f"kept {len(kept)}")

def contacts(cif):
    atoms = pdbx.get_structure(pdbx.CIFFile.read(str(cif)), model=1, use_author_fields=True)
    protein = atoms[(~atoms.hetero) & (atoms.chain_id == "A")]
    ligand = atoms[(atoms.res_name == "6NT") & (atoms.element != "H")]
    residue_ids = []
    for residue in protein[protein.atom_name == "CA"].res_id:
        if int(residue) not in residue_ids:
            residue_ids.append(int(residue))
    def xyz(res_index, atom_name, table):
        res_id = residue_ids[res_index] if table is protein else None
        selected = table[(table.res_id == (residue_ids[res_index] if table is protein else table.res_id)) & (table.atom_name == atom_name)] if False else None
        return selected
    def point(res_index, atom_name):
        res_id = residue_ids[res_index]
        selected = protein[(protein.res_id == res_id) & (protein.atom_name == atom_name)]
        return selected.coord[0]
    def lig(name):
        selected = ligand[ligand.atom_name == name]
        return selected.coord[0]
    base = float(np.linalg.norm(point(40, "OD2") - lig("N3")))
    donor = float(np.linalg.norm(point(64, "NE2") - lig("N1")))
    return base, donor

assets = Path(os.environ["BOLTZ_ASSETS"])
rows = []
for index, (header, seq) in enumerate(kept, start=1):
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
    rows.append({"sequence": index, "ASP_OD2_to_6NT_N3": base, "GLN_NE2_to_6NT_N1": donor, "under_5": base < 5 and donor < 5})
    print(f"seq {index}: base {base:.2f} donor {donor:.2f}")
(run / "contacts.json").write_text(json.dumps(rows, indent=2) + "\n")
print("wrote", run / "contacts.json")
PY
