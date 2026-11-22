#!/bin/bash
# APNovo provides the ORGANISERS' ancestral scorer. Using it removes a whole
# class of divergence from reimplementing angle2/dihedral2/burial by hand.
# Heavy deps (alphafold3 from git + jax[cuda12]) so it gets its own uv venv.
set -u
D=/home/mandalorian/Hackathons/Enzymes-BSSE2026
A=$D/tools/alphaprotein-novo
V=$D/.venv-apnovo
export UV_CACHE_DIR=/home/mandalorian/.cache/uv
echo "[$(date '+%F %T')] creating uv venv (py3.12)"
uv venv --python 3.12 "$V" 2>&1 | tail -3
echo "[$(date '+%F %T')] installing alphaprotein-novo (this pulls alphafold3 + jax)"
VIRTUAL_ENV="$V" uv pip install --python "$V/bin/python" -e "$A" 2>&1 | tail -25
echo "rc=$?"
echo "[$(date '+%F %T')] import check"
"$V/bin/python" - <<'PY' 2>&1 | tail -12
try:
    from alphaprotein_novo.metrics import kemp_metrics
    print("  OK kemp_metrics")
    print("  entry:", hasattr(kemp_metrics, "kemp_structure_metrics"))
    from alphafold3 import structure
    print("  OK alphafold3.structure")
except Exception as e:
    print(f"  FAIL {type(e).__name__}: {str(e)[:300]}")
PY
echo "### APNOVO INSTALL DONE ###"
