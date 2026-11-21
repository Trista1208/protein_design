#!/usr/bin/env bash
#SBATCH --job-name=kemp-finish
#SBATCH --partition=gpu
#SBATCH --constraint=g6.xlarge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=01:00:00

# Finish one Kemp design from the only sequence that has already
# placed both catalytic contacts under 5 A.
# Parent: five new Boltz samples.
# C21S and C21A: the same sequence with the unpaired cysteine removed.

set -e -o pipefail
source "${PDH_WORKSPACE:?Submit from the project directory}/env.sh"
export PDH_RUN_MODE=project
RUN="$PDH_WORKSPACE/results/kemp-finish-$SLURM_JOB_ID"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"
mkdir -p "$RUN"

echo "Job: $SLURM_JOB_ID; parent plus cysteine-free point designs, five Boltz samples each"
"$ENV_ROOT/env-rfd3/bin/python" - "$RUN" << 'PY'
import json, os, shutil, subprocess, sys
from pathlib import Path
import numpy as np
from biotite.structure.io import pdbx

run = Path(sys.argv[1])
parent = "DKRGKRVEMTEEEKEAFEKGCKKANLEELLNKEWTVYESEDVIGGKKVRKIIAVTEYNGKIYLLQLLLVDGELVAVLVDVVDTKGTPKIETTTLEYRIPAAFGSKAPAIAAVYAAMMEKPLLFNLVSQDGDVSIFRGTLEVDGRTLNFTMEVIWENGVPVAIRVTDDEGIG"
assert len(parent) == 171 and parent[20] == "C" and parent[40] == "D" and parent[64] == "Q"
variants = {
    "parent": parent,
    "C21S": parent[:20] + "S" + parent[21:],
    "C21A": parent[:20] + "A" + parent[21:],
}
assets = Path(os.environ["BOLTZ_ASSETS"])
rows = []

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

for name, seq in variants.items():
    seq_run = run / name
    msa = seq_run / "msa"
    msa.mkdir(parents=True)
    (msa / "kemp.a3m").write_text(f">{name}\n{seq}\n")
    (seq_run / "kemp.fasta").write_text(f">{name}\n{seq}\n")
    yaml = seq_run / "kemp.yaml"
    yaml.write_text(
        "version: 1\nsequences:\n- protein:\n    id: A\n    sequence: " + seq
        + "\n    msa: " + str((msa / "kemp.a3m").resolve())
        + "\n- ligand:\n    id: B\n    ccd: 6NT\n"
    )
    cache = seq_run / "boltz-cache"
    cache.mkdir()
    for asset in ("boltz2_conf.ckpt", "boltz2_aff.ckpt", "mols.tar", "mols"):
        (cache / asset).symlink_to(assets / asset)
    log = seq_run / "boltz.log"
    print(f"Boltz {name}", flush=True)
    with log.open("w") as handle:
        proc = subprocess.run([
            os.environ["ENV_ROOT"] + "/env-boltz2/bin/boltz", "predict", str(yaml),
            "--model", "boltz2", "--cache", str(cache), "--out_dir", str(seq_run / "boltz"),
            "--accelerator", "gpu", "--devices", "1", "--no_kernels",
            "--sampling_steps", "200", "--recycling_steps", "3", "--diffusion_samples", "5",
            "--num_workers", "0", "--preprocessing-threads", "1", "--output_format", "mmcif",
            "--seed", "8",
        ], stdout=handle, stderr=subprocess.STDOUT)
    pred = seq_run / "boltz/boltz_results_kemp/predictions/kemp"
    if proc.returncode != 0:
        rows.append({"variant": name, "error": "boltz failed"})
        print(f"{name}: boltz failed")
        continue
    for model in range(5):
        cif = pred / f"kemp_model_{model}.cif"
        if not cif.is_file():
            rows.append({"variant": name, "model": model, "error": "missing"})
            continue
        base, donor = contacts(cif)
        confidence_path = pred / f"confidence_kemp_model_{model}.json"
        confidence = json.loads(confidence_path.read_text()) if confidence_path.is_file() else {}
        row = {
            "variant": name,
            "model": model,
            "sequence": seq,
            "position_21": seq[20],
            "ASP_OD2_to_6NT_N3": base,
            "GLN_NE2_to_6NT_N1": donor,
            "under_5": base < 5 and donor < 5,
            "complex_plddt": confidence.get("complex_plddt"),
            "complex_ipde": confidence.get("complex_ipde"),
            "ligand_iptm": confidence.get("ligand_iptm"),
            "cif": str(cif),
        }
        rows.append(row)
        print(f"{name} model {model}: base {base:.2f} donor {donor:.2f} plddt {row['complex_plddt']}")

passing = [row for row in rows if row.get("under_5")]
passing.sort(key=lambda row: (row["sequence"].count("C"), row["ASP_OD2_to_6NT_N3"] + row["GLN_NE2_to_6NT_N1"]))
(run / "contacts.json").write_text(json.dumps(rows, indent=2) + "\n")
if passing:
    best = passing[0]
    finished = run / "finished"
    finished.mkdir()
    shutil.copy(best["cif"], finished / "kemp_model.cif")
    (finished / "kemp.fasta").write_text(">" + best["variant"] + "\n" + best["sequence"] + "\n")
    (finished / "summary.json").write_text(json.dumps(best, indent=2) + "\n")
    print("FINISHED", best["variant"], "model", best["model"], "base", round(best["ASP_OD2_to_6NT_N3"], 2), "donor", round(best["GLN_NE2_to_6NT_N1"], 2))
else:
    print("NO PASS")
print("wrote", run / "contacts.json")
PY
