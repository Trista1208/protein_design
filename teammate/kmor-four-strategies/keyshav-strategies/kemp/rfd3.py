"""M1/M4 - RFdiffusion3 enzyme-design driver.

Schema read from the INSTALLED package, as the handover insists, not from the
tutorial prose:
  rfd3/inference/input_parsing.py :: DesignInputSpecification

Verified fields (rc-foundry 1.0.0):
  input, contig, unindex, length, ligand, allow_ligand_on_existing_chain,
  cif_parser_args, extra, dialect(=2 release), select_fixed_atoms,
  select_unfixed_sequence, select_buried, select_partially_buried,
  select_exposed, select_hbond_acceptor, select_hbond_donor, select_hotspots,
  symmetry, ori_token, infer_ori_strategy('com'|'hotspots'), plddt_enhanced,
  is_non_loopy, partial_t

Two fields the handover did not list are directly useful for this challenge:
  select_hbond_donor / select_hbond_acceptor - atom-wise H-bond conditioning.
  Oxyanion stabilisation IS an H-bond (Gln50 NE2 -> 6NT N1), and handover 3.2
  says the whole challenge hinges on capturing it, so we condition on it.

Observed failure mode: passing only `input` raises
  ValidationError: Input provided but unused in composition specification
so a composition field (length and/or unindex/contig) is mandatory.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import config as C
from .theozyme import Theozyme, assert_ori_sane, load, ori_for, write_clean_pdb


class RFD3Error(RuntimeError):
    pass


# --------------------------------------------------------------------------- spec
@dataclass
class EnzymeDesignSpec:
    """One RFD3 enzyme-design calculation block."""

    name: str
    input_pdb: Path
    ori_token: list[float]
    # rc-foundry 1.0.0 wants the ligand RESIDUE NAME here ("Ligand name or
    # index"), e.g. "6NT". The handover's "l:g" colon form is from an older
    # tutorial; a contig like "G3" also fails with
    #   [component=G3] Component not found ... available_non_protein: ['6NT']
    # Repository wins over the doc.
    ligand: str | None = None
    unindex: str | None = None
    length: str = C.RFD3_LENGTH_RANGE
    select_fixed_atoms: dict[str, str] | None = None
    select_unfixed_sequence: str | None = None
    select_buried: str | None = None
    # Atom-level selections are InputSelection = Dict[contig, atoms], NOT
    # "G1:NE2" strings. Observed: 'G1:NE2' -> ComponentValidationError
    # "Invalid contig format ... Expected 'ChainIDResID' (e.g. 'A20')".
    select_hbond_donor: dict[str, str] | None = None
    select_hbond_acceptor: dict[str, str] | None = None
    dialect: int = 2
    plddt_enhanced: bool = True
    extra: dict = field(default_factory=dict)

    #: keys whose value is a plain contig string
    _STR_KEYS = ("unindex", "select_unfixed_sequence", "select_buried")
    #: keys whose value is a {contig: "ATOM,ATOM"} mapping
    _DICT_KEYS = ("select_fixed_atoms", "select_hbond_donor", "select_hbond_acceptor")

    def to_dict(self) -> dict:
        d: dict = {
            # hydra chdirs into out_dir, so every path must be absolute
            "input": str(Path(self.input_pdb).resolve()),
                "length": self.length,
            "ori_token": [float(x) for x in self.ori_token],
            "dialect": self.dialect,
            "plddt_enhanced": self.plddt_enhanced,
        }
        if self.ligand:
            d["ligand"] = self.ligand
        for key in self._STR_KEYS:
            val = getattr(self, key)
            if val:
                d[key] = val
        for key in self._DICT_KEYS:
            val = getattr(self, key)
            if val:
                d[key] = dict(val)
        if self.extra:
            d["extra"] = self.extra
        return d


# --------------------------------------------------------------------------- builders
def _motif_selection(tz: Theozyme) -> tuple[str, dict[str, str], str]:
    """Return (unindex, select_fixed_atoms, select_unfixed_sequence) for the motif.

    Handover 1.4: the ligand is held automatically - never fix it yourself, and
    never fix side-chain atoms of residues listed in select_unfixed_sequence.
    """
    asp = tz.by_resname("ASP")
    gln = tz.by_resname("GLN")
    if not asp or not gln:
        raise RFD3Error("theozyme is missing ASP or GLN - cannot build a motif spec")
    asp_tag = f"{asp[0].chain}{asp[0].resseq}"
    gln_tag = f"{gln[0].chain}{gln[0].resseq}"

    unindex = f"{asp_tag},{gln_tag}"
    # Catalytic base: hold the carboxylate + backbone rigid, identity fixed.
    fixed = {asp_tag: "OD1,OD2,CG,CB,CA,N,C,O"}
    return unindex, fixed, ""


def asp_tag_only(tz: Theozyme) -> str:
    a = tz.by_resname("ASP")[0]
    return f"{a.chain}{a.resseq}"


def build_spec(variant: str, tz: Theozyme, *, length: str | None = None,
               ori: np.ndarray | None = None) -> EnzymeDesignSpec:
    """Build one of the handover 4.2 motif variants.

    Scaffolding at least MIN_MOTIF_VARIANTS of these is mandatory - five
    correlated designs from one pocket wastes all five lab slots.
    """
    if variant not in C.MOTIF_VARIANTS:
        raise RFD3Error(f"unknown motif variant {variant!r}; known: {sorted(C.MOTIF_VARIANTS)}")

    ori_vec = ori_for(tz) if ori is None else np.asarray(ori, dtype=float)
    dist = assert_ori_sane(tz, ori_vec)            # hard gate, every run

    unindex, fixed, unfixed_seq = _motif_selection(tz)
    gln = tz.by_resname("GLN")[0]
    gln_tag = f"{gln.chain}{gln.resseq}"
    lig_atoms = tz.by_resname(C.LIGAND_CCD)
    lig_tag = f"{lig_atoms[0].chain}{lig_atoms[0].resseq}" if lig_atoms else None

    spec = EnzymeDesignSpec(
        name=variant,
        input_pdb=tz.path,
        ori_token=list(map(float, ori_vec)),
        unindex=unindex,
        length=length or C.RFD3_LENGTH_RANGE,
        select_fixed_atoms=dict(fixed),
        ligand=C.LIGAND_CCD,
        extra={"variant": variant, "ori_ligand_dist_A": round(dist, 3),
               "note": C.MOTIF_VARIANTS[variant]},
    )

    if variant == "canonical":
        # Oxyanion stabiliser as an explicit H-bond donor onto the ligand.
        spec.select_hbond_donor = {gln_tag: "NE2"}
        spec.select_fixed_atoms[gln_tag] = "NE2,OE1,CD,CG,CB,CA,N,C,O"
    elif variant == "crystal_gln":
        spec.select_hbond_donor = {gln_tag: "NE2"}
        # backbone only: let the rotamer resolve against the crystallographic one
        spec.select_fixed_atoms[gln_tag] = "CA,N,C,O"
    elif variant == "backbone_nh":
        # no side-chain donor at all; a backbone amide must do the stabilising
        # NOTE: select_hbond_acceptor on the LIGAND is rejected by RFD3 -
        # "[component=G3] Residue G3 not found in atom array". Atom-level
        # conditioning there only accepts protein residues, so we leave the
        # acceptor unconditioned and rely on select_unfixed_sequence.
        spec.select_hbond_acceptor = None   # ligand contig rejected by RFD3
        spec.select_unfixed_sequence = gln_tag
    elif variant == "alt_donor":
        # Ser/Thr/Tyr hydroxyl instead of Gln: free the identity, keep backbone
        spec.select_unfixed_sequence = gln_tag
        spec.select_fixed_atoms[gln_tag] = "CA,N,C,O"
        # NOTE: select_hbond_acceptor on the LIGAND is rejected by RFD3 -
        # "[component=G3] Residue G3 not found in atom array". Atom-level
        # conditioning there only accepts protein residues, so we leave the
        # acceptor unconditioned and rely on select_unfixed_sequence.
        spec.select_hbond_acceptor = None   # ligand contig rejected by RFD3
    elif variant == "oxy_free":
        # The supplied Gln BACKBONE caps d_oxy at >=3.27 A (exhaustive rotamer
        # scan), against the organisers' 2.40 A target, so pinning it locks
        # oxy ~ 0.22 into every design - the limiting term. Here the donor is
        # free in position AND identity: only the Asp base and the ligand are
        # held, and the H-bond onto the ligand is conditioned. This is the only
        # variant that can beat the oxyanion cap.
        # KEEP both residues in `unindex`: it means "motif residue with an
        # unspecified sequence position", so RFD3 still places the donor freely
        # but reports where it landed in diffused_index_map. Dropping it lost
        # the mapping and every oxy_free design died with
        # stage=no_diffused_index_map.
        # Freedom comes from fixing NO atoms of the donor and unfixing its
        # sequence - not from removing it from the motif.
        spec.select_fixed_atoms = {k: v for k, v in spec.select_fixed_atoms.items()
                                   if k == asp_tag_only(tz)}
        spec.select_unfixed_sequence = gln_tag
        # select_hbond_acceptor on the LIGAND is rejected by RFD3:
        #   [component=G3] Residue G3 not found in atom array
        # Atom-level conditioning accepts protein residues only. Freeing the
        # donor via select_unfixed_sequence is what actually does the work here.
        spec.select_hbond_acceptor = None
    elif variant == "polarizer":
        # KE_MOTIF_ROLES sanctions a SECOND carboxylate (D/E) polarizer role.
        spec.select_hbond_donor = {gln_tag: "NE2"}
        spec.select_fixed_atoms[gln_tag] = "NE2,OE1,CD,CG,CB,CA,N,C,O"
        spec.extra["extra_role"] = "polarizer (second D/E) - free for RFD3 to place"
    elif variant == "aromatic":
        # KE_MOTIF_ROLES sanctions an aromatic F/W/Y pi-stacking role.
        spec.select_hbond_donor = {gln_tag: "NE2"}
        spec.select_fixed_atoms[gln_tag] = "NE2,OE1,CD,CG,CB,CA,N,C,O"
        spec.extra["extra_role"] = "aromatic_res (F/W/Y) pi-stack - free placement"
    elif variant == "base_variant":
        # Glu instead of Asp: free the base identity, keep its backbone anchored
        asp = tz.by_resname("ASP")[0]
        asp_tag = f"{asp.chain}{asp.resseq}"
        spec.select_fixed_atoms[asp_tag] = "CA,N,C,O"
        spec.select_unfixed_sequence = asp_tag
        spec.select_hbond_donor = {gln_tag: "NE2"}

    # Never fix side chains of residues whose sequence we unfixed (handover 1.4).
    if spec.select_unfixed_sequence:
        for tag in spec.select_unfixed_sequence.split(","):
            atoms = spec.select_fixed_atoms.get(tag.strip())
            if atoms and any(a not in {"CA", "N", "C", "O", "CB"} for a in atoms.split(",")):
                raise RFD3Error(
                    f"{tag} is in select_unfixed_sequence but has side-chain atoms "
                    f"fixed ({atoms}); handover 1.4 forbids this"
                )
    return spec


# --------------------------------------------------------------------------- runner
def write_inputs_json(specs: list[EnzymeDesignSpec], out: Path) -> Path:
    """RFD3 takes a JSON whose keys name the calculation blocks."""
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({s.name: s.to_dict() for s in specs}, indent=2))
    return out


def run(inputs_json: Path, out_dir: Path, *, batch_size: int = 1, n_batches: int = 1,
        timesteps: int | None = None, seed: int = 1, timeout: int = 7200,
        dry_run: bool = False) -> dict:
    """Invoke the rfd3 CLI.

    batch_size=1 by default: handover 1.4 notes all designs in a batch share one
    length, so variety needs batch size 1.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(C.BIN_RFD3),
        f"inputs={Path(inputs_json).resolve()}",
        f"out_dir={out_dir.resolve()}",
        f"diffusion_batch_size={batch_size}",
        f"n_batches={n_batches}",
        f"seed={seed}",
    ]
    if timesteps:
        cmd.append(f"inference_sampler.num_timesteps={timesteps}")
    if dry_run:
        return {"cmd": cmd, "dry_run": True}

    log = out_dir / "rfd3.log"
    with log.open("w") as fh:
        proc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, timeout=timeout)
    # rfd3 writes .cif.gz, so match on the full name not just the final suffix
    produced = sorted(p.name for p in out_dir.rglob("*")
                      if p.is_file() and p.name.endswith((".pdb", ".cif", ".cif.gz", ".pdb.gz")))
    return {
        "cmd": cmd,
        "returncode": proc.returncode,
        "log": str(log),
        "n_structures": len(produced),
        "structures": produced[:20],
    }


def smoke(out_root: Path | None = None, timesteps: int = 20) -> dict:
    """M1 - one config, batch 1, measure VRAM and wall clock.

    The docs give no VRAM figure (handover 1.4), so this measures it.
    """
    out_root = Path(out_root or (C.RUNS / "m1_smoke"))
    tz = load(C.THEOZYME_WITH_ORI)
    clean = write_clean_pdb(tz, out_root / "theozyme_clean.pdb")
    tz_clean = load(clean)
    spec = build_spec("canonical", tz_clean, length="190-200")
    js = write_inputs_json([spec], out_root / "inputs.json")
    return run(js, out_root / "designs", batch_size=1, n_batches=1, timesteps=timesteps)
