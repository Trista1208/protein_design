#!/usr/bin/env bash
#SBATCH --job-name=kemp-5s
#SBATCH --partition=gpu
#SBATCH --constraint=g6.xlarge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=00:30:00

set -e -o pipefail
source "${PDH_WORKSPACE:?Submit from the project directory}/env.sh"
export PDH_RUN_MODE=project
RUN="$PDH_WORKSPACE/results/kemp-5seed-$SLURM_JOB_ID"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
mkdir -p "$RUN/msa" "$RUN/boltz-cache"

echo "Job: $SLURM_JOB_ID; five Boltz seeds of kemp-seq-21 sequence 4"
"$ENV_ROOT/env-rfd3/bin/python" - "$RUN" << 'PY'
import json, os, subprocess, sys
from pathlib import Path
import numpy as np
from biotite.structure.io import pdbx

run = Path(sys.argv[1])
seq = "DKRGKRVEMTEEEKEAFEKGCKKANLEELLNKEWTVYESEDVIGGKKVRKIIAVTEYNGKIYLLQLLLVDGELVAVLVDVVDTKGTPKIETTTLEYRIPAAFGSKAPAIAAVYAAMMEKPLLFNLVSQDGDVSIFRGTLEVDGRTLNFTMEVIWENGVPVAIRVTDDEGIG"
assert len(seq) == 171 and seq[40] == "D" and seq[64] == "Q"
(run / "msa/kemp.a3m").write_text(">kemp_design\n" + seq + "\n")
yaml = run / "kemp.yaml"
yaml.write_text(
    "version: 1\nsequences:\n- protein:\n    id: A\n    sequence: " + seq
    + "\n    msa: " + str((run / "msa/kemp.a3m").resolve())
    + "\n- ligand:\n    id: B\n    ccd: 6NT\n"
)
cache = run / "boltz-cache"
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
        "--sampling_steps", "200", "--recycling_steps", "3", "--diffusion_samples", "5",
        "--num_workers", "0", "--preprocessing-threads", "1", "--output_format", "mmcif",
        "--seed", "1",
    ], stdout=handle, stderr=subprocess.STDOUT)
if proc.returncode != 0:
    raise SystemExit(f"Boltz failed, see {log}")

rows = []
for model in range(5):
    cif = run / "boltz/boltz_results_kemp/predictions/kemp" / f"kemp_model_{model}.cif"
    if not cif.is_file():
        rows.append({"model": model, "error": "missing"})
        continue
    atoms = pdbx.get_structure(pdbx.CIFFile.read(str(cif)), model=1, use_author_fields=True)
    protein = atoms[(~atoms.hetero) & (atoms.chain_id == "A")]
    ligand = atoms[(atoms.res_name == "6NT") & (atoms.element != "H")]
    def pick(res_index, atom_name):
        residue_ids = []
        for residue in protein[protein.atom_name == "CA"].res_id:
            if residue not in residue_ids:
                residue_ids.append(int(residue))
        res_id = residue_ids[res_index]
        selected = protein[(protein.res_id == res_id) & (protein.atom_name == atom_name)]
        return selected.coord[0]
    def ligand_atom(name):
        selected = ligand[ligand.atom_name == name]
        return selected.coord[0]
    base = float(np.linalg.norm(pick(40, "OD2") - ligand_atom("N3")))
    donor = float(np.linalg.norm(pick(64, "NE2") - ligand_atom("N1")))
    rows.append({"model": model, "ASP_OD2_to_6NT_N3": base, "GLN_NE2_to_6NT_N1": donor, "under_5": base < 5 and donor < 5})
    print(f"model {model}: base {base:.2f} donor {donor:.2f}")
(run / "contacts.json").write_text(json.dumps(rows, indent=2) + "\n")
print("wrote", run / "contacts.json")
PY
