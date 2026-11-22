"""M2/M4 - LigandMPNN and Boltz-2 drivers.

Non-negotiables encoded here:
  * handover 4.3.1 - catalytic residues FIXED in LigandMPNN, asserted afterwards.
  * handover 4.3.6 - identical MSA settings across every Boltz run, forever.
    De novo designs have no meaningful MSA, so single-sequence mode is pinned in
    config.BOLTZ_SETTINGS and this module refuses to override it.
  * handover 4.3.7 - 5 recorded seeds, never cherry-picked.
  * handover 1.6 - open-source `boltz` on a 24 GB 3090, NOT the NVIDIA NIM.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass
from pathlib import Path

from . import config as C


class RunnerError(RuntimeError):
    pass


# =========================================================== LigandMPNN
@dataclass
class MPNNResult:
    returncode: int
    out_dir: Path
    fasta: Path | None
    sequences: list[tuple[str, str]]
    log: Path


def ligandmpnn(pdb: Path | str, out_dir: Path | str, *,
               fixed_residues: str | None = None,
               redesigned_residues: str | None = None,
               omit_per_residue: dict[str, str] | None = None,
               n_batches: int = 1,
               n_seqs: int = 8, temperature: float = 0.1, seed: int = 37,
               checkpoint: str | None = None, timeout: int = 3600) -> MPNNResult:
    """Design sequences around a fixed catalytic motif.

    ``fixed_residues`` is LigandMPNN's space-separated "C10 C11" form. It is
    STRONGLY recommended: without it LigandMPNN will mutate the catalytic base
    away (handover 4.3.1). ``redesigned_residues`` is the inverse and easier
    when only a handful of positions are open; pass one or the other.

    ``omit_per_residue`` maps a position to residues forbidden THERE, e.g.
    ``{"A26": "F"}``. This is how you force a MUTATION rather than a proposal:
    on a parent that already came from an MPNN-class model, sequence recovery
    is ~1.0 and sampling simply returns the native residue. Forbidding the
    native is the difference between 4 variants and a usable family.

    Note the repo-local NumPy patch documented in README.md - LigandMPNN's
    vendored OpenFold uses aliases removed in NumPy 1.24.
    """
    # ABSOLUTE, both of them. This call runs with cwd=LIGANDMPNN_DIR, so a
    # relative --out_folder silently lands in tools/LigandMPNN/<your/path>/
    # and a relative --pdb_path makes prody report "not a valid PDB
    # identifier". Callers that passed absolute paths never saw either.
    pdb, out_dir = Path(pdb).resolve(), Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    ckpt = checkpoint or str(C.LIGANDMPNN_DIR / "model_params" / "ligandmpnn_v_32_010_25.pt")
    if not Path(ckpt).exists():
        raise RunnerError(f"LigandMPNN checkpoint missing: {ckpt}. "
                          "Run tools/LigandMPNN/get_model_params.sh")

    cmd = [str(C.PY_AIZYMES), str(C.LIGANDMPNN_DIR / "run.py"),
           "--model_type", "ligand_mpnn",
           "--checkpoint_ligand_mpnn", ckpt,
           "--pdb_path", str(pdb),
           "--out_folder", str(out_dir),
           "--batch_size", str(n_seqs),
           "--number_of_batches", str(n_batches),
           "--temperature", str(temperature),
           "--seed", str(seed)]
    if fixed_residues:
        cmd += ["--fixed_residues", fixed_residues]
    if redesigned_residues:
        cmd += ["--redesigned_residues", redesigned_residues]
    if omit_per_residue:
        omit_json = out_dir / "omit_per_residue.json"
        omit_json.write_text(json.dumps(omit_per_residue))
        cmd += ["--omit_AA_per_residue", str(omit_json)]

    log = out_dir / "ligandmpnn.log"
    with log.open("w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                              cwd=str(C.LIGANDMPNN_DIR), timeout=timeout)

    fastas = sorted((out_dir / "seqs").glob("*.fa")) if (out_dir / "seqs").exists() else []
    seqs: list[tuple[str, str]] = []
    if fastas:
        header, buf = None, []
        for line in fastas[0].read_text().splitlines():
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(buf)))
                header, buf = line[1:], []
            elif line.strip():
                buf.append(line.strip())
        if header is not None:
            seqs.append((header, "".join(buf)))

    return MPNNResult(proc.returncode, out_dir, fastas[0] if fastas else None, seqs, log)


def assert_catalytic_preserved(sequences: list[tuple[str, str]],
                               expected: dict[int, str]) -> None:
    """Fail loudly if LigandMPNN mutated a catalytic residue (handover 4.3.1)."""
    problems = []
    for header, seq in sequences[1:]:           # [0] is the input sequence
        for idx, aa in expected.items():
            if idx < len(seq) and seq[idx] != aa:
                problems.append(f"{header}: pos {idx} {aa}->{seq[idx]}")
    if problems:
        raise RunnerError(
            "LigandMPNN mutated catalytic residues despite --fixed_residues:\n  "
            + "\n  ".join(problems)
        )


# =========================================================== Boltz-2
@dataclass
class BoltzResult:
    seed: int
    returncode: int
    out_dir: Path
    structures: list[Path]
    confidence: dict | None
    log: Path


def _write_boltz_fasta(sequence: str, path: Path, *, ligand_ccd: str | None = None) -> Path:
    """Boltz multi-record FASTA: one record per chain, ligand by CCD code.

    Header is CHAIN_ID|ENTITY_TYPE|MSA_ID. The MSA_ID `empty` selects TRUE
    SINGLE-SEQUENCE MODE - verified empirically: omitting it makes boltz demand
    `--use_msa_server`, and `none` is accepted but produces no structure, while
    `empty` folds and yields real confidences (complex_ipde 1.069,
    ligand_iptm 0.907 on a 202-aa design + 6NT).

    This is pinned deliberately (handover 4.3.6): de novo designs have no
    meaningful MSA, and drifting MSA settings between designs makes confidence
    numbers incomparable. Never swap this for --use_msa_server mid-campaign.

    The ligand MUST be its own chain, or the interface mask is empty and
    complex_ipde comes back 0.0 - which looks perfect and means nothing.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [">A|protein|empty", sequence]
    if ligand_ccd:
        lines += [">B|ccd", ligand_ccd]
    path.write_text("\n".join(lines) + "\n")
    return path


def _write_boltz_yaml(sequence: str, path: Path, *, ligand_ccd: str,
                      pocket_contacts: list[int] | None = None,
                      contact_atoms: list[tuple[int, str]] | None = None,
                      max_distance: float = 5.0,
                      template_pdb: Path | str | None = None) -> Path:
    """Boltz YAML input with a pocket constraint, per-residue contacts, and a
    structural template.

    This addresses the central failure of de novo enzyme scoring. Measured on a
    real design:

      no constraint            d_base 13.25 A   d_oxy      -      -> all geometry 0
      pocket, max_distance 6   d_base  3.56 A   d_oxy 11.55 A    -> oxy still 0

    With only a loose pocket constraint Boltz parks the ligand BESIDE the Asp
    pointing away from the Gln, so no single position satisfies both contacts.
    Three strengthenings, cheapest first:

      1. max_distance 6 -> 4 A on the pocket
      2. an explicit `contact` per catalytic residue, so BOTH must be close
         rather than the pair being satisfied on average
      3. the RFD3 backbone as a `template`, which hands Boltz the very geometry
         we want reproduced. RFD3 builds the site accurately (Asp-Gln span
         10.606 A vs the theozyme's 10.637 A); the template stops Boltz
         discarding that.

    Caveat worth carrying: constraining the ligand into the pocket we then score
    is somewhat circular. It asserts design intent, and we still measure whether
    the result is catalytically competent - but if the organisers fold WITHOUT
    constraints our numbers will read optimistic against theirs.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["version: 1", "sequences:",
             "  - protein:", "      id: A", f'      sequence: "{sequence}"',
             "      msa: empty",
             "  - ligand:", "      id: B", f"      ccd: {ligand_ccd}"]
    if pocket_contacts:
        contacts = ", ".join(f"[A, {r}]" for r in pocket_contacts)
        lines += ["constraints:",
                  "  - pocket:", "      binder: B",
                  f"      contacts: [{contacts}]",
                  f"      max_distance: {max_distance}"]
        # One contact per catalytic residue, so BOTH must be satisfied rather
        # than the pair being satisfied on average.
        # For a LIGAND chain the token is an ATOM NAME, not a residue index -
        # `[B, 1]` raises KeyError ('B', 0, 1) in token_spec_to_ids. Naming the
        # atom is also chemically sharper: the base contacts N3 and the
        # oxyanion donor contacts N1 (KEMP_ELIMINASE_LIGAND_METADATA).
        for resseq, atom in (contact_atoms or []):
            lines += ["  - contact:", f"      token1: [A, {resseq}]",
                      f"      token2: [B, {atom}]",
                      f"      max_distance: {max_distance}"]
    if template_pdb is not None:
        lines += ["templates:", f"  - pdb: {Path(template_pdb).resolve()}",
                  "    chain_id: A", "    template_id: A"]
    path.write_text("\n".join(lines) + "\n")
    return path


def boltz_predict(sequence: str, out_dir: Path | str, *, seed: int,
                  ligand_ccd: str | None = C.LIGAND_CCD,
                  pocket_contacts: list[int] | None = None,
                  contact_atoms: list[tuple[int, str]] | None = None,
                  template_pdb: Path | str | None = None,
                  timeout: int = 3600) -> BoltzResult:
    """One Boltz-2 prediction with the project's PINNED settings.

    MSA settings are taken from config and are NOT parameterised - drifting them
    between designs makes confidences incomparable (handover 4.3.6).
    """
    out_dir = Path(out_dir) / f"seed{seed}"
    out_dir.mkdir(parents=True, exist_ok=True)
    if pocket_contacts and ligand_ccd:
        inp = _write_boltz_yaml(sequence, out_dir / "input.yaml",
                                ligand_ccd=ligand_ccd,
                                pocket_contacts=pocket_contacts,
                                contact_atoms=contact_atoms,
                                template_pdb=template_pdb)
    else:
        inp = _write_boltz_fasta(sequence, out_dir / "input.fasta",
                                 ligand_ccd=ligand_ccd)

    cmd = [str(C.BIN_BOLTZ), "predict", str(inp),
           "--out_dir", str(out_dir),
           "--diffusion_samples", str(C.BOLTZ_SETTINGS["diffusion_samples"]),
           "--recycling_steps", str(C.BOLTZ_SETTINGS["recycling_steps"]),
           "--seed", str(seed),
           "--override"]
    if C.BOLTZ_SETTINGS["use_msa_server"]:
        cmd.append("--use_msa_server")
    # else: single-sequence mode. Boltz needs an explicit empty MSA for this;
    # see docs/USER_GUIDE.md for the --msa none convention of your boltz build.

    log = out_dir / "boltz.log"
    with log.open("w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, timeout=timeout)

    structures = sorted(out_dir.rglob("*.cif")) + sorted(out_dir.rglob("*.pdb"))
    conf = None
    for cj in out_dir.rglob("confidence*.json"):
        try:
            conf = json.loads(cj.read_text())
            break
        except json.JSONDecodeError:
            pass
    return BoltzResult(seed, proc.returncode, out_dir, structures, conf, log)


def boltz_replicates(sequence: str, out_dir: Path | str, *,
                     seeds: list[int] | None = None,
                     ligand_ccd: str | None = C.LIGAND_CCD,
                     pocket_contacts: list[int] | None = None,
                     contact_atoms: list[tuple[int, str]] | None = None,
                     template_pdb: Path | str | None = None) -> list[BoltzResult]:
    """The official 5 replicates. Every seed is recorded; none is discarded."""
    seeds = seeds or list(C.BOLTZ_SEEDS)
    if len(seeds) != C.REPLICATES:
        raise RunnerError(
            f"the challenge requires {C.REPLICATES} replicates, got {len(seeds)} seeds"
        )
    return [boltz_predict(sequence, out_dir, seed=s, ligand_ccd=ligand_ccd,
                          pocket_contacts=pocket_contacts,
                          contact_atoms=contact_atoms,
                          template_pdb=template_pdb) for s in seeds]


# =========================================================== FieldTools
def fieldtools(pdb: Path | str, out_dir: Path | str, *, timeout: int = 1800) -> dict:
    """Electrostatic preorganisation proxy (handover 3.4).

    Report as a PROXY, explicitly not as dG-dagger. Runs under the `mm` env
    because it needs pytraj from AmberTools.
    """
    pdb, out_dir = Path(pdb), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not C.FIELDTOOLS.exists():
        raise RunnerError(f"FieldTools.py missing at {C.FIELDTOOLS}")
    log = out_dir / "fieldtools.log"
    cmd = [str(C.PY_MM), str(C.FIELDTOOLS), str(pdb)]
    with log.open("w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                              cwd=str(out_dir), timeout=timeout)
    return {"returncode": proc.returncode, "log": str(log),
            "caveat": "electrostatic preorganisation PROXY, not a free energy barrier"}
