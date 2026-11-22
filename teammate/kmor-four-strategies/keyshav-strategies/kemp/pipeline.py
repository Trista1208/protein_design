"""M2 / M4 / M5 - the parts that actually chain the tools together.

Design notes:

* RFD3 emits ``.cif.gz``; LigandMPNN wants a PDB. ``cif_to_pdb`` handles that
  with whichever of gemmi/biotite is importable, rather than assuming one.
* Every stage writes a row to a JSONL ledger so the live reporter can publish
  partial results while the sweep is still running.
* M4 is a generator: it yields after each design so results stream out instead
  of appearing only at the end.
"""

from __future__ import annotations

import gzip
import json
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterator

from . import config as C
from . import apnovo, gates, novelty, objective, official, rfd3, scoring, strain
from .runners import (RunnerError, assert_catalytic_preserved, boltz_predict,
                      boltz_replicates, fieldtools, ligandmpnn)
from .theozyme import load, write_clean_pdb

LEDGER = C.OUTPUTS / "ledger.jsonl"


class PipelineError(RuntimeError):
    pass


# --------------------------------------------------------------------------- ledger
def append_ledger(row: dict, path: Path | None = None) -> None:
    """Append one result row. The live reporter tails this file."""
    path = path or LEDGER
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {"ts": time.strftime("%F %T"), **row}
    with path.open("a") as fh:
        fh.write(json.dumps(row, default=str) + "\n")


def read_ledger(path: Path | None = None) -> list[dict]:
    path = path or LEDGER
    if not path.exists():
        return []
    rows = []
    for line in path.read_text().splitlines():
        if line.strip():
            try:
                rows.append(json.loads(line))
            except json.JSONDecodeError:
                pass
    return rows


# --------------------------------------------------------------------------- format
def cif_to_pdb(src: Path | str, dst: Path | str) -> Path:
    """Convert (optionally gzipped) mmCIF to PDB."""
    src, dst = Path(src), Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)

    plain = src
    if src.name.endswith(".gz"):
        plain = dst.with_suffix(".cif")
        with gzip.open(src, "rb") as fi, plain.open("wb") as fo:
            shutil.copyfileobj(fi, fo)

    # try gemmi, then biotite, in whichever env has one
    for env in ("foundry", "rfd4", "mm"):
        py = C.ENVS / env / "bin/python"
        if not py.exists():
            continue
        code = (
            "import sys\n"
            "src, dst = sys.argv[1], sys.argv[2]\n"
            "try:\n"
            "    import gemmi\n"
            "    st = gemmi.read_structure(src); st.setup_entities()\n"
            "    st.write_pdb(dst); print('gemmi')\n"
            "except ImportError:\n"
            "    from biotite.structure.io.pdbx import CIFFile, get_structure\n"
            "    from biotite.structure.io.pdb import PDBFile\n"
            "    arr = get_structure(CIFFile.read(src), model=1)\n"
            "    f = PDBFile(); f.set_structure(arr); f.write(dst); print('biotite')\n"
        )
        proc = subprocess.run([str(py), "-c", code, str(plain), str(dst)],
                              capture_output=True, text=True, timeout=600)
        if proc.returncode == 0 and dst.exists():
            return dst
    raise PipelineError(f"could not convert {src} to PDB (tried gemmi and biotite)")


def sequence_from_pdb(pdb: Path | str) -> str:
    """One-letter sequence from CA records, chain order preserved."""
    three = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
        "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
        "TYR": "Y", "VAL": "V",
    }
    out = []
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out.append(three.get(line[17:20].strip(), "X"))
    return "".join(out)


def catalytic_site_from_rfd3(out_dir: Path, *, base_tag: str,
                             oxy_tag: str) -> official.CatalyticSite | None:
    """Recover where the motif landed from RFD3's output JSON.

    RFD3 writes ``diffused_index_map`` e.g. {"G2": "A95", "G1": "A39"}. Without
    this, geometry on a real design cannot be measured at all (65% of the
    official score), so a missing map is a hard failure, not a warning.
    """
    for js in sorted(out_dir.rglob("*_model_*.json")):
        try:
            d = json.loads(js.read_text())
        except json.JSONDecodeError:
            continue
        dim = d.get("diffused_index_map")
        if isinstance(dim, dict) and base_tag in dim and oxy_tag in dim:
            return official.CatalyticSite.from_diffused_index_map(
                dim, base_tag=base_tag, oxy_tag=oxy_tag)
    return None


# --------------------------------------------------------------------------- M2
@dataclass
class DesignRecord:
    design_id: str
    variant: str
    backbone: str = ""
    sequence: str = ""
    length: int = 0
    gates: dict = field(default_factory=dict)
    novelty: dict = field(default_factory=dict)
    geometry: dict = field(default_factory=dict)
    burial: dict = field(default_factory=dict)
    ipde: dict = field(default_factory=dict)
    consistency: dict = field(default_factory=dict)
    fieldtools: dict = field(default_factory=dict)
    official: dict = field(default_factory=dict)
    catalytic: dict = field(default_factory=dict)
    backbone_gate: dict = field(default_factory=dict)
    strain: dict = field(default_factory=dict)
    #: fold protocol that produced `official`. A score without it is not a
    #: measurement: the same sequence and seed scores 0.000 single-sequence and
    #: 0.775 with an MSA (LESSONS 3.1-3.2).
    protocol: dict = field(default_factory=dict)
    #: dense continuous defects per official term, defined where the official
    #: sub-scores are a flat 0. Internal ranking only (kemp/objective.py).
    dense: dict = field(default_factory=dict)
    #: mean/min complex pLDDT and the C.PLDDT_MIN verdict. A design below it is
    #: unmeasurable rather than bad, and is excluded from ranking.
    foldability: dict = field(default_factory=dict)
    stage: str = "generated"
    eligible: bool | None = None
    score: float | None = None
    wall_s: float = 0.0


def design_sequences(backbone_pdb: Path, out_dir: Path, *, n_seqs: int = 8,
                     fixed_residues: str | None = None,
                     catalytic: dict[int, str] | None = None) -> list[tuple[str, str]]:
    """LigandMPNN with the catalytic motif pinned, then assert it survived."""
    res = ligandmpnn(backbone_pdb, out_dir, n_seqs=n_seqs,
                     fixed_residues=fixed_residues)
    if res.returncode != 0:
        raise PipelineError(f"LigandMPNN failed (rc={res.returncode}); see {res.log}")
    if catalytic:
        assert_catalytic_preserved(res.sequences, catalytic)
    return res.sequences[1:]          # [0] is the input sequence


def prefilter_sequence(rec: "DesignRecord", backbone: Path, *,
                       catalytic: dict[int, str] | None = None) -> dict:
    """CPU-only eligibility gates, run BEFORE paying 80 s of GPU for a fold.

    Both novelty gates need inputs that already exist at this point:
      * MMseqs2 identity needs only the SEQUENCE (LigandMPNN just produced it)
      * Foldseek TM needs only a STRUCTURE (the RFD3 backbone serves)

    Running them after folding - as the first version did - meant a design that
    is 85% identical to HG3.17 still consumed a full Boltz run before being
    thrown away.
    """
    out: dict = {}
    out["gates"] = gates.run_all(rec.sequence, motif_residues=2, catalytic=catalytic)
    try:
        out["novelty"] = novelty.screen(rec.sequence, backbone)
    except Exception as exc:                        # noqa: BLE001
        out["novelty"] = {"eligible": False, "verdict": "DISQUALIFIED",
                          "reasons": [f"novelty screen error: {exc}"]}
    out["passed"] = bool(out["gates"].get("passed")) and bool(
        out["novelty"].get("eligible"))
    return out


def score_design(rec: DesignRecord, structure: Path, *, replicates: int = 1,
                 out_dir: Path, site: official.CatalyticSite | None = None,
                 ligand_ccd: str | None = C.LIGAND_CCD) -> DesignRecord:
    """Fold, then score with the ORGANISERS' weights.

    Geometry is measured on EVERY replicate, because the official Consistency
    term is sd_base - the standard deviation of d_base across replicates. A
    pipeline that folds 5 seeds but measures geometry once cannot compute it.
    """
    t0 = time.time()
    # Stamp the fold protocol BEFORE folding, so it is recorded even when the
    # scoring below fails. LESSONS 3.2: "Record the MSA setting alongside every
    # score. This is how the control and the designs silently became
    # incomparable."
    rec.protocol = objective.protocol_stamp(replicates=replicates)
    # Constrain the ligand to the DESIGNED pocket. Without this Boltz docks it
    # blind and lands ~13 A away, zeroing all four geometry terms.
    pocket = ([site.base_resseq, site.oxy_resseq] if site is not None else None)
    # Per-residue atom contacts were tested and REJECTED: with both present
    # Boltz trades one off against the other (d_base 3.56 -> 13.68 A while
    # d_oxy 11.55 -> 6.65 A). A single pocket constraint naming both residues
    # is one coherent instruction and gave the best result. The RFD3 backbone
    # as a template also failed - Boltz's template parser raises IndexError on
    # it (isolated: pocket OK, pocket+atoms OK, pocket+template IndexError).
    catoms = None
    tmpl = None   # template path breaks Boltz; see note above
    runs = ([boltz_predict(rec.sequence, out_dir, seed=C.BOLTZ_SEEDS[0],
                           ligand_ccd=ligand_ccd, pocket_contacts=pocket,
                           contact_atoms=catoms, template_pdb=tmpl)]
            if replicates == 1 else
            boltz_replicates(rec.sequence, out_dir, ligand_ccd=ligand_ccd,
                             pocket_contacts=pocket, contact_atoms=catoms,
                             template_pdb=tmpl))

    ipde_values: list[float] = []
    plddt_values: list[float] = []
    per_replicate: list[official.RawMetrics] = []
    for r in runs:
        if r.confidence:
            try:
                _, val = scoring.extract_ipde(r.confidence)
                ipde_values.append(val)
            except scoring.ScoringError as exc:
                rec.ipde.setdefault("warnings", []).append(str(exc)[:160])
            pl = r.confidence.get("complex_plddt")
            if pl is not None:
                plddt_values.append(float(pl))
        folded = next((p for p in r.structures), None)
        if folded is None:
            continue
        try:
            pdb = (folded if folded.suffix == ".pdb"
                   else cif_to_pdb(folded, r.out_dir / "folded.pdb"))
            # Measure with APNovo's own code, not a reimplementation. Hand-
            # rolled angle2/dihedral2 gave 328.1/50.9 deg where theirs gives
            # 162.8/168.5 - the 50.9 fell outside the sp2 band and scored
            # 0.000, which under a geometric mean would have zeroed every
            # design. Always pass motif_residues: inference keys on
            # atom_b_factor > 0 and the Gln has B=0.00.
            if apnovo.available() and site is not None:
                try:
                    am = apnovo.metrics(
                        pdb, motif_residues=[site.base_resseq, site.oxy_resseq])
                    per_replicate.append(apnovo.to_raw(am))
                except apnovo.ApNovoError as exc:
                    per_replicate.append(official.RawMetrics(
                        notes=[f"apnovo failed on seed {r.seed}: {str(exc)[:200]}"]))
            else:
                tz = load(pdb, require_ligand=False)
                per_replicate.append(official.measure_raw(tz, site=site))
        except Exception as exc:                 # noqa: BLE001
            per_replicate.append(official.RawMetrics(
                notes=[f"replicate {r.seed} unmeasurable: {type(exc).__name__}: {exc}"]))

    if per_replicate:
        raw = official.combine_replicates(per_replicate, ipde_values=ipde_values)
        rec.official = official.score(raw)
        # Dense defects, recorded for EVERY design whatever the official score.
        # LESSONS 7.3: "Report every term, every design, every replicate -
        # especially the dihedrals. Six of eight designs died on a term no
        # table showed." The per-replicate dihedrals go in too, because the
        # Consistency term proves single-shot geometry is not enough.
        rec.dense = {
            **objective.diagnostics(raw),
            "badness": objective.badness(raw),
            "per_replicate_dihedral_oxy_deg": [
                None if r.dihedral12_oxy is None else round(r.dihedral12_oxy, 2)
                for r in per_replicate],
            "per_replicate_dihedral_base_deg": [
                None if r.dihedral12_base is None else round(r.dihedral12_base, 2)
                for r in per_replicate],
            "per_replicate_d_base_A": [
                None if r.d_base is None else round(r.d_base, 3)
                for r in per_replicate],
        }
        rec.score = (rec.official.get("official_score")
                     if rec.official.get("official_score") is not None
                     else rec.official.get("prefilter_score"))
        rec.geometry = {k: rec.official["detail"].get(k) for k in
                        ("base", "oxy", "orient_base", "orient_oxy")}
        rec.burial = rec.official["detail"].get("bury", {})
        rec.consistency = rec.official["detail"].get("cons", {})
        rec.ipde = {**rec.ipde, **(rec.official["detail"].get("ipde") or {})}
    else:
        rec.official = {"official_score": None,
                        "verdict": "INCOMPLETE - no measurable folded replicate"}

    # gates/novelty were already run by prefilter_sequence() before folding;
    # only fill them in if this was called directly
    if not rec.gates:
        rec.gates = gates.run_all(rec.sequence, motif_residues=2)
    if not rec.novelty:
        try:
            rec.novelty = novelty.screen(rec.sequence, structure)
        except Exception as exc:                 # noqa: BLE001
            rec.novelty = {"eligible": False, "verdict": "DISQUALIFIED",
                           "reasons": [f"novelty screen error: {exc}"]}

    # Foldability gate. A design below the threshold is UNMEASURABLE, not bad
    # (LESSONS 7.1), so it is recorded with a named reason rather than silently
    # dropped - rejection reasons must be countable into a histogram
    # (LESSONS 4.5), or you learn the campaign died without learning why.
    if plddt_values:
        mean_plddt = sum(plddt_values) / len(plddt_values)
        rec.foldability = {
            "mean_plddt": round(mean_plddt, 4),
            "min_plddt": round(min(plddt_values), 4),
            "per_replicate": [round(p, 4) for p in plddt_values],
            "threshold": C.PLDDT_MIN,
            "passed": mean_plddt >= C.PLDDT_MIN,
            "rationale": ("below this the geometry terms are noise, not "
                          "measurement; see config.PLDDT_MIN"),
        }
    else:
        rec.foldability = {"mean_plddt": None, "threshold": C.PLDDT_MIN,
                           "passed": False,
                           "rationale": "no complex_plddt reported - cannot "
                                        "establish the fold is measurable"}
    if not rec.foldability["passed"]:
        rec.stage = (f"unmeasurable:plddt_{rec.foldability['mean_plddt']}"
                     f"_below_{C.PLDDT_MIN}")

    rec.eligible = (bool(rec.gates.get("passed"))
                    and bool(rec.novelty.get("eligible"))
                    and bool(rec.foldability.get("passed"))
                    and rec.score is not None)
    rec.wall_s = round(time.time() - t0, 1)
    if rec.foldability["passed"]:
        rec.stage = f"scored_x{replicates}"
    return rec


# --------------------------------------------------------------------------- M4
def scale(n_designs: int, *, variants: list[str] | None = None,
          lengths: list[str] | None = None, seqs_per_backbone: int = 8,
          prefilter_only: bool = True, out_root: Path | None = None,
          timesteps: int = 200,
          accept_all_backbones: bool = False,
          checkpoint_after: int = objective.CAMPAIGN_CHECKPOINT_AFTER,
          ) -> Iterator[DesignRecord]:
    """Generate, design and prefilter ``n_designs`` backbones.

    Yields after each design so the live reporter can publish partial results.
    Boltz x5 is deliberately NOT run here - that is for finalists only
    (handover 4.1). Use ``finalise`` for that.

    ``checkpoint_after`` scored designs in, the sweep audits its own score
    columns and ABORTS if any is constant across designs - the LESSONS 2.5
    failure, where RFD3's inherited geometry was ranked on for a whole
    campaign. Set it to 0 to disable (a deliberate choice, e.g. a 3-design
    pilot where no distribution exists to judge).
    """
    out_root = Path(out_root or (C.RUNS / "m4"))
    out_root.mkdir(parents=True, exist_ok=True)
    # oxy_free added deliberately: every design so far has been zeroed by
    # oxy=0.000, and oxy_free is the only variant that frees the donor backbone
    # instead of pinning the supplied one (which caps d_oxy at >=3.27 A).
    variants = variants or ["canonical", "oxy_free", "crystal_gln",
                            "backbone_nh", "alt_donor", "base_variant"]
    lengths = lengths or ["190-210", "210-230", "230-250"]

    tz = load(C.THEOZYME_DESIGN if C.THEOZYME_DESIGN.exists()
              else C.THEOZYME_WITH_ORI)
    clean = write_clean_pdb(tz, out_root / "theozyme_clean.pdb")
    tz_clean = load(clean)

    made = 0
    scored_badness: list[dict] = []
    checkpoint_done = checkpoint_after <= 0
    for i in range(n_designs):
        variant = variants[i % len(variants)]
        length = lengths[(i // len(variants)) % len(lengths)]
        did = f"d{i:05d}_{variant}_{length.replace('-', '_')}"
        asp0 = tz_clean.by_resname("ASP")[0]
        gln0 = tz_clean.by_resname("GLN")[0]
        base_tag = f"{asp0.chain}{asp0.resseq}"
        oxy_tag = f"{gln0.chain}{gln0.resseq}"
        rec = DesignRecord(design_id=did, variant=variant)
        wdir = out_root / did
        try:
            spec = rfd3.build_spec(variant, tz_clean, length=length)
            d = spec.to_dict()
            # KEEP `ligand`: with the ligand on its own chain RFD3 can hold it in
            # place. Popping it made RFD3 build backbone THROUGH the ligand
            # volume (0.46 A closest approach), so Boltz had to dock 6NT
            # elsewhere and every geometry term scored 0.
            js = wdir / "inputs.json"
            js.parent.mkdir(parents=True, exist_ok=True)
            js.write_text(json.dumps({variant: d}, indent=2))

            res = rfd3.run(js, wdir / "rfd3", batch_size=1, n_batches=1,
                           timesteps=timesteps, seed=1000 + i)
            if res["returncode"] != 0 or res["n_structures"] == 0:
                rec.stage = "rfd3_failed"
                append_ledger(asdict(rec))
                yield rec
                continue

            cifs = sorted((wdir / "rfd3").rglob("*.cif.gz")) or \
                   sorted((wdir / "rfd3").rglob("*.cif"))
            bb = cif_to_pdb(cifs[0], wdir / "backbone.pdb")
            rec.backbone = str(bb)

            # FREE pre-Boltz gate on RFD3's own metrics. The official scalar is a
            # geometric mean, so a misplaced motif scores a hard 0.000 however
            # good the rest is - reject it before paying for MPNN + folding.
            rfd3_metrics = {}
            for js in sorted((wdir / "rfd3").rglob("*_model_*.json")):
                try:
                    rfd3_metrics = json.loads(js.read_text()).get("metrics", {}) or {}
                    break
                except json.JSONDecodeError:
                    pass
            bg = official.gate_backbone(rfd3_metrics)
            rec.backbone_gate = bg
            if not bg["passed"] and not accept_all_backbones:
                rec.stage = "backbone_rejected:" + ",".join(bg["failed"])
                append_ledger(asdict(rec))
                yield rec
                continue

            site = catalytic_site_from_rfd3(wdir / "rfd3", base_tag=base_tag,
                                            oxy_tag=oxy_tag)
            if site is None:
                rec.stage = "no_diffused_index_map"
                append_ledger(asdict(rec))
                yield rec
                continue
            rec.catalytic = {"base": site.base_resseq, "oxy": site.oxy_resseq,
                             "chain": site.chain}

            # MEASURED structural validation of the RFD3 output, before MPNN.
            # gate_backbone above trusts RFD3's self-reported metrics; this
            # measures the coordinates. It catches what that cannot:
            #   - chain cuts RFD3 reports as n_chainbreaks=0 (observed: a 2.68 A
            #     C-N bond reported as zero breaks)
            #   - bond/angle strain in the catalytic residues (never checked)
            #   - d_base / d_oxy on the design itself
            # Run on the CIF, not backbone.pdb: the ligand chain survives only
            # in the CIF, and without it d_base/d_oxy cannot be measured.
            val = strain.validate(cifs[0], base_resseq=site.base_resseq,
                                  oxy_resseq=site.oxy_resseq, chain=site.chain)
            rec.strain = val.as_dict()
            if not val.passed and not accept_all_backbones:
                rec.stage = "strain_rejected:" + ",".join(val.failed)
                append_ledger(asdict(rec))
                yield rec
                continue

            fixed = f"{site.chain}{site.base_resseq} {site.chain}{site.oxy_resseq}"
            seqs = design_sequences(bb, wdir / "lmpnn", n_seqs=seqs_per_backbone,
                                    fixed_residues=fixed)
            if not seqs:
                rec.stage = "mpnn_failed"
                append_ledger(asdict(rec))
                yield rec
                continue

            header, seq = seqs[0]
            rec.sequence, rec.length = seq, len(seq)

            # CPU gates FIRST - no GPU spent on a design that cannot qualify
            cat = {site.base_resseq - 1: "D", site.oxy_resseq - 1: "Q"}
            pre = prefilter_sequence(rec, bb, catalytic=cat)
            rec.gates, rec.novelty = pre["gates"], pre["novelty"]
            if not pre["passed"]:
                why = (rec.gates.get("verdict"), rec.novelty.get("reasons", [])[:1])
                rec.stage = f"prefilter_rejected:{why}"[:200]
                rec.eligible = False
                append_ledger(asdict(rec))
                yield rec
                continue

            rec = score_design(rec, bb, replicates=1, out_dir=wdir / "boltz",
                               site=site)
            made += 1
            if rec.dense.get("badness"):
                scored_badness.append(rec.dense["badness"])
        except Exception as exc:          # noqa: BLE001 - a sweep must not die on one design
            rec.stage = f"error: {type(exc).__name__}: {exc}"[:300]
        append_ledger(asdict(rec))
        yield rec

        # Early-abort checkpoint, fired ONCE. A degenerate score column means
        # the pipeline is reporting its input rather than measuring designs
        # (LESSONS 2.5), and every further design is wasted GPU time. Checked
        # here rather than before launch because at launch there is nothing to
        # judge. Advisory findings are logged and the sweep continues.
        if not checkpoint_done and len(scored_badness) >= checkpoint_after:
            checkpoint_done = True
            cp = objective.campaign_checkpoint(scored_badness,
                                               after=checkpoint_after)
            (out_root / "checkpoint.json").write_text(json.dumps(cp, indent=2))
            if cp["advisory"]:
                print(f"[checkpoint] advisory, sweep continues: "
                      f"{cp['advisory']} often missing", flush=True)
            if not cp["passed"]:
                raise objective.DegenerateScoresError(
                    f"campaign aborted after {cp['n_scored']} scored designs: "
                    f"term(s) {cp['fatal']} are constant across designs. "
                    f"{cp['reason']}. See {out_root / 'checkpoint.json'}. "
                    "Fix the measurement before spending more GPU time.")
            print(f"[checkpoint] {cp['n_scored']} designs scored, "
                  "no degenerate term; sweep continues", flush=True)


def finalise(records: list[DesignRecord], *, top_n: int = 30,
             out_root: Path | None = None) -> list[DesignRecord]:
    """Boltz x5 on the best prefiltered, eligible designs only."""
    out_root = Path(out_root or (C.RUNS / "m4_final"))
    # a zeroed geometry term means the official geometric mean is 0.000, so
    # finalising such a design would waste 5 folds on a guaranteed zero
    eligible = [r for r in records if r.eligible and r.score is not None
                and not (r.official or {}).get("zeroed_terms")]
    eligible.sort(key=lambda r: r.score, reverse=True)
    out = []
    for rec in eligible[:top_n]:
        wdir = out_root / rec.design_id
        rec = score_design(rec, Path(rec.backbone), replicates=C.REPLICATES,
                           out_dir=wdir)
        append_ledger(asdict(rec))
        out.append(rec)
    return out


# --------------------------------------------------------------------------- M5
def select_final(records: list[DesignRecord], *, n: int = C.N_LAB_SLOTS,
                 noise_floor: float = 0.0) -> dict:
    """Diversity first, then expressibility and electrostatics as tie-breaks.

    Handover 3.5: of 400 proteins tested in a comparable competition, 95%
    expressed but only 14% bound. Diversity across the slots beats depth.
    """
    pool = [r for r in records if r.eligible and r.score is not None]
    if not pool:
        return {"selected": [], "reason": "no eligible scored designs"}

    pool.sort(key=lambda r: r.score, reverse=True)

    # one design per motif variant first - that is the diversity axis we control
    chosen: list[DesignRecord] = []
    seen_variants: set[str] = set()
    for rec in pool:
        if rec.variant not in seen_variants:
            chosen.append(rec)
            seen_variants.add(rec.variant)
        if len(chosen) == n:
            break
    # top up by score if we ran out of variants
    for rec in pool:
        if len(chosen) == n:
            break
        if rec not in chosen:
            chosen.append(rec)

    indistinguishable = []
    for a, b in zip(chosen, chosen[1:]):
        if not scoring.resolvable(a.score, b.score, noise_floor):
            indistinguishable.append((a.design_id, b.design_id))

    return {
        "selected": [r.design_id for r in chosen],
        "variants_covered": sorted(seen_variants),
        "n_variants": len(seen_variants),
        "min_variants_required": C.MIN_MOTIF_VARIANTS,
        "diversity_ok": len(seen_variants) >= min(C.MIN_MOTIF_VARIANTS, n),
        "noise_floor": noise_floor,
        "pairs_within_noise": indistinguishable,
        "caveat": ("pairs listed above are indistinguishable from chance; "
                   "ranking within the noise floor is superstition (handover 3.2)"),
        "records": [asdict(r) for r in chosen],
    }


def run_fieldtools_on(records: list[DesignRecord], out_root: Path | None = None) -> None:
    """Electrostatic preorganisation proxy on the finalists."""
    out_root = Path(out_root or (C.RUNS / "m5_field"))
    for rec in records:
        if not rec.backbone:
            continue
        try:
            rec.fieldtools = fieldtools(rec.backbone, out_root / rec.design_id)
        except (RunnerError, Exception) as exc:   # noqa: BLE001
            rec.fieldtools = {"error": str(exc)}
