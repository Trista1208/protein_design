#!/usr/bin/env bash
#SBATCH --job-name=kemp-bb
#SBATCH --partition=gpu
#SBATCH --constraint=g6.xlarge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=00:20:00

set -e -o pipefail
source "${PDH_WORKSPACE:?Submit from the project directory}/env.sh"
export PDH_RUN_MODE=project
PROJECT="$PDH_WORKSPACE/src"
PARENT="$PDH_WORKSPACE/results/kemp-18"
RUN="$PDH_WORKSPACE/results/kemp-bb-$SLURM_JOB_ID"
export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

echo "Job: $SLURM_JOB_ID; Boltz of the RFD3 backbone sequence"
"$ENV_ROOT/env-rfd3/bin/python" - "$PARENT" "$RUN" "$PROJECT" << 'PY'
import json, os, subprocess, sys
from pathlib import Path
parent, run, project = map(Path, sys.argv[1:4])
seq = "DLHGAEVSLSAATESAITTGLAKASLSSLLNQAFTVYESADVIGGVAATLKTAITAYKGKTYLLQLSLISNELISLTLSVVDLEANPQSVIQVLSASIPAALGAHAVAVQQLFDAAQQTPIVYTLASEDGSVKTYTGQITVDGETFSFTKLEYYVDGVPVSIIVTDDKNIG"
assert len(seq) == 171 and seq[40] == "D" and seq[64] == "Q"
bridge = run / "bridge"
bridge.mkdir(parents=True)
(bridge / "kemp.pdb").write_bytes((parent / "bridge/kemp.pdb").read_bytes())
(bridge / "motif.json").write_bytes((parent / "bridge/motif.json").read_bytes())
(bridge / "sequence.json").write_text(json.dumps({
    "sequence": seq,
    "length": 171,
    "fixed": [
        {"original": "G50", "id": "A65", "position": 65, "amino_acid": "Q"},
        {"original": "G127", "id": "A41", "position": 41, "amino_acid": "D"},
    ],
    "ligandmpnn_scores": {"id": "backbone", "note": "RFD3 sequence, not redesigned"},
}, indent=2) + "\n")
msa = run / "msa"
msa.mkdir()
(msa / "kemp.a3m").write_text(">kemp_design\n" + seq + "\n")
(msa / "verification.json").write_text(json.dumps({"kemp_design": {"depth": 1, "nonquery_rows": 0, "query_length": 171}}, indent=2) + "\n")
(run / "kemp.yaml").write_text(
    "version: 1\nsequences:\n- protein:\n    id: A\n    sequence: " + seq + "\n    msa: "
    + str((msa / "kemp.a3m").resolve()) + "\n- ligand:\n    id: B\n    ccd: 6NT\n"
)
cache = run / "boltz-cache"
cache.mkdir()
assets = Path(os.environ["BOLTZ_ASSETS"])
for name in ("boltz2_conf.ckpt", "boltz2_aff.ckpt", "mols.tar", "mols"):
    (cache / name).symlink_to(assets / name)
log = run / "boltz.log"
with log.open("w") as handle:
    proc = subprocess.run([
        os.environ["ENV_ROOT"] + "/env-boltz2/bin/boltz", "predict", str(run / "kemp.yaml"),
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
