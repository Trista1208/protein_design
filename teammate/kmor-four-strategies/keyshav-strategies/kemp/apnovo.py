"""Bridge to the organisers' ancestral scorer (alphaprotein-novo).

Why this exists, in one paragraph: reimplementing `angle2` / `dihedral2` by hand
produced plausible-looking numbers that were simply the wrong metric. Measured
on the supplied theozyme, a hand-rolled version gave `dihedral2` = 328.1 deg at
the base and 50.9 deg at the oxyanion; APNovo's own code gives 162.8 and 168.5.
The hand-rolled oxyanion value fell 0.9 deg outside the sp2 scoring band and
scored 0.000 - which, under a weighted geometric mean, would have zeroed EVERY
design in the campaign. The reference atoms are not guessable: APNovo picks
`angle_atom2 = L1C3A` at the base and `L1C7A` at the oxyanion via
`get_reference_atoms_and_hybridization`, not the N2 a chemist would reach for.

So: measure with their code, score with our window layer on top.

APNovo lives in its own uv venv (`.venv-apnovo`) because it drags in alphafold3
and jax[cuda12], which must not touch the pipeline envs. This module shells out
to that interpreter and returns the metrics dict.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from . import config as C

VENV = C.ROOT / ".venv-apnovo"
PY = VENV / "bin/python"
PROBE = C.ROOT / "scripts" / "apnovo_metrics.py"

#: Mapping from APNovo metric keys to our RawMetrics fields.
FIELD_MAP = {
    "d_base": "acid_base/distance",
    "d_oxy": "oxyanion_oah/distance",
    "angle2_base": "acid_base/angle2",
    "dihedral12_base": "acid_base/dihedral2",
    "angle2_oxy": "oxyanion_oah/angle2",
    "dihedral12_oxy": "oxyanion_oah/dihedral2",
    "ca_neighbors_mean": "cat_base_tip_ca_neighbors_mean",
}


class ApNovoError(RuntimeError):
    pass


def available() -> bool:
    ccd = (VENV / "lib/python3.12/site-packages/alphafold3/constants/converters"
                  "/ccd.msgpack")
    return PY.exists() and PROBE.exists() and ccd.exists()


def metrics(pdb: Path | str, *, motif_residues: list[int] | str | None = None,
            ligand_res_name: str = C.LIGAND_CCD,
            timeout: int = 900) -> dict:
    """Run `kemp_structure_metrics` on ``pdb`` and return the raw metric dict.

    ``motif_residues`` must ALWAYS be passed for real work. When it is None,
    APNovo infers the motif from ``atom_b_factor > 0``; the supplied Gln50 has
    B = 0.00 on all 17 atoms and is therefore silently dropped, taking `oxy`
    (0.20) and `orient_oxy` (0.08) with it. Under a geometric mean a dropped
    term that scores 0 zeroes the design outright.
    """
    if not available():
        raise ApNovoError(
            f"alphaprotein-novo bridge unavailable (need {PY}, {PROBE} and "
            "alphafold3's ccd.msgpack - run scripts/install_apnovo.sh then "
            "`.venv-apnovo/bin/build_data`)")
    cmd = [str(PY), str(PROBE), "--pdb", str(pdb),
           "--ligand", ligand_res_name, "--json"]
    if motif_residues is not None:
        spec = (motif_residues if isinstance(motif_residues, str)
                else ",".join(str(i) for i in motif_residues))
        cmd += ["--motif", spec]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    # Parse stdout FIRST: the worker exits 1 on a measurement failure but puts
    # the real reason in the JSON payload. Checking the return code first threw
    # that detail away and surfaced an empty "apnovo metrics failed:".
    payload = None
    if proc.stdout.strip():
        try:
            payload = json.loads(proc.stdout.strip().splitlines()[-1])
        except (json.JSONDecodeError, IndexError):
            payload = None
    if payload is not None and "error" in payload:
        raise ApNovoError(payload["error"])
    if proc.returncode != 0 or payload is None:
        raise ApNovoError(
            f"apnovo metrics failed (rc={proc.returncode}); "
            f"stderr={proc.stderr[-300:]!r} stdout={proc.stdout[-300:]!r}")
    return payload


def to_raw(apnovo_metrics: dict, *, ipde: float | None = None):
    """Translate APNovo metric keys into our RawMetrics."""
    from .official import RawMetrics

    m = RawMetrics()
    for ours, theirs in FIELD_MAP.items():
        v = apnovo_metrics.get(theirs)
        if isinstance(v, (int, float)):
            setattr(m, ours, float(v))
    m.ipde = ipde
    m.notes.append("metrics measured by alphaprotein-novo "
                   "kemp_structure_metrics (authoritative), not reimplemented")
    # record which reference atoms their code chose, so a surprising
    # orientation score can be traced to the atom selection
    for k in ("acid_base/angle_atom2", "acid_base/dihedral_atom1",
              "acid_base/dihedral_atom2", "oxyanion_oah/angle_atom2",
              "oxyanion_oah/dihedral_atom1", "oxyanion_oah/dihedral_atom2",
              "acid_base/hyb_atom2", "oxyanion_oah/hyb_atom2"):
        if k in apnovo_metrics:
            m.notes.append(f"{k}={apnovo_metrics[k]}")
    for k in ("acid_base/hbond", "oxyanion_oah/hbond"):
        if k in apnovo_metrics:
            m.notes.append(f"{k}={apnovo_metrics[k]}")
    return m
