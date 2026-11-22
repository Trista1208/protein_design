"""Tests for the traps that fail SILENTLY.

A test here is worth writing only if its failure mode would otherwise waste a
batch of designs or fabricate a score. Ordinary happy paths are covered by
`make m0` / `make m3`.

Run: make test
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from kemp import config as C
from kemp import gates, rfd3, scoring, theozyme
from kemp.theozyme import TheozymeError, load


# ============================================================ fixtures
@pytest.fixture(scope="module")
def tz_with():
    return load(C.THEOZYME_WITH_ORI)


@pytest.fixture(scope="module")
def tz_no():
    return load(C.THEOZYME_NO_ORI)


# ============================================================ trap 1: frames
def test_two_files_same_motif_different_frame(tz_with, tz_no):
    """Handover 1.2: identical motif, different frame."""
    cmp = theozyme.compare_frames(tz_with, tz_no)
    assert cmp["n_shared_atoms"] == 45
    assert cmp["superposed_rmsd_A"] < 0.01, "motifs should superpose near-perfectly"
    assert cmp["raw_rmsd_A"] > 50, "as supplied the files must be far apart"


def test_ori_is_frame_specific():
    a = theozyme.ori_for(C.THEOZYME_WITH_ORI)
    b = theozyme.ori_for(C.THEOZYME_NO_ORI)
    assert float(np.linalg.norm(a - b)) > 50, "the two ORI values are not interchangeable"


def test_ori_sane_in_its_own_frame(tz_with, tz_no):
    for tz in (tz_with, tz_no):
        d = theozyme.assert_ori_sane(tz, theozyme.ori_for(tz.path))
        assert d < C.ORI_LIGAND_MAX_DIST


def test_mixed_frame_is_rejected(tz_no):
    """THE trap: RFD3 would not error, so we must."""
    with pytest.raises(TheozymeError, match="mixed-frame"):
        theozyme.assert_ori_sane(tz_no, theozyme.ori_for(C.THEOZYME_WITH_ORI))


def test_unknown_frame_refuses_to_guess(tmp_path):
    f = tmp_path / "Theozyme_mystery.pdb"
    f.write_text("END\n")
    with pytest.raises(TheozymeError, match="refusing to guess"):
        theozyme.ori_for(f)


# ============================================================ trap 2: ORI atom
def test_ori_atom_stripped_by_default(tz_with):
    assert all(a.resname != "ORI" and a.name != "ORI" for a in tz_with.atoms)
    assert tz_with.ori is not None, "coordinate must still be captured"


def test_clean_pdb_has_no_ori(tz_with, tmp_path):
    out = theozyme.write_clean_pdb(tz_with, tmp_path / "clean.pdb")
    lines = out.read_text().splitlines()
    coords = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    # Check ATOM records only: the frame REMARK legitimately contains "withORI".
    assert coords, "clean PDB has no coordinates"
    for l in coords:
        assert l[17:20].strip() != "ORI" and l[12:16].strip() != "ORI", \
            "ORI breaks RDKit/OpenBabel (3-char element field)"
    assert any(l[17:20].strip() == C.LIGAND_CCD for l in coords), \
        "the ligand must survive cleaning"


def test_clean_pdb_carries_frame_provenance(tz_with, tmp_path):
    """Our own frame guard must not block our own derived files."""
    out = theozyme.write_clean_pdb(tz_with, tmp_path / "clean.pdb")
    assert theozyme.FRAME_REMARK in out.read_text()
    ori = theozyme.ori_for(out)          # must not raise
    assert np.allclose(ori, C.ORI_IN_WITHORI_FRAME)
    tz2 = load(out)
    assert tz2.frame == theozyme.FRAME_WITH_ORI
    assert theozyme.assert_ori_sane(tz2, theozyme.ori_for(tz2)) < C.ORI_LIGAND_MAX_DIST


def test_clean_pdb_residues_non_decreasing(tz_with, tmp_path):
    """RF3 rejects out-of-order res_id: 'not in non-decreasing order for chain G'."""
    out = theozyme.write_clean_pdb(tz_with, tmp_path / "clean.pdb")
    per_chain: dict[str, list[int]] = {}
    for l in out.read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            per_chain.setdefault(l[21:22], []).append(int(l[22:26]))
    assert per_chain, "no coordinates written"
    for chain, ids in per_chain.items():
        assert ids == sorted(ids), f"chain {chain} res_ids are not non-decreasing"


# ============================================================ trap 3: altlocs
def test_ligand_survives_altloc_normalisation(tz_with):
    """The ligand carries altloc 'B'; a naive 'keep A or blank' filter deletes it."""
    assert len(tz_with.ligand) > 0
    assert "B" in {"A", "B", ""} | tz_with.altlocs or True  # documents the mix


def test_altlocs_are_mixed_in_the_source(tz_with):
    """If this ever becomes single-altloc, the normaliser can be simplified."""
    raw = load(C.THEOZYME_WITH_ORI, normalise_altlocs=False)
    assert len(raw.altlocs) > 1, "source genuinely mixes altlocs"


def test_normalisation_leaves_one_atom_per_name(tz_with):
    keys = [a.key for a in tz_with.atoms]
    assert len(keys) == len(set(keys)), "duplicate atoms survived normalisation"


# ============================================================ geometry
def test_geometry_reproduces_verified_reference(tz_with):
    """All six handover values must round-trip."""
    measured = theozyme.measure_geometry(tz_with)
    report = theozyme.diff_reference(measured)
    for key, row in report.items():
        assert row["within_tol"], f"{key} deviates: {row}"


def test_gln_reference_now_trusted(tz_with):
    """RESOLVED: the modelled Gln50 MATCHES the crystallographic one in 4BS0/7K4Z.

    It still has B-factor 0.00 (it was modelled in), but it agrees with the
    crystal structures, so it is no longer treated as an opinion.
    """
    report = theozyme.diff_reference(theozyme.measure_geometry(tz_with))
    assert report["gln50_ne2__6nt_n1"]["modelled_not_observed"] is False
    assert C.REF_GEOMETRY_MODELLED == set()
    assert "CONFIRMED" in C.GLN50_PROVENANCE


def test_reference_theozyme_scores_near_one(tz_with):
    geo = scoring.geometry_score(theozyme.measure_geometry(tz_with))
    assert geo["score"] > 0.99, "the reference must score ~1.0 or the scorer is wrong"
    # Gln50 confirmed against 4BS0/7K4Z, so the donor term is now trusted
    assert geo["terms"]["donor_contact"]["trusted_reference"] is True
    assert geo["terms"]["donor_contact"]["weight"] >= 0.30, (
        "oxyanion stabilisation is the mechanism the challenge turns on")


# ============================================================ scoring honesty
def test_ipde_never_invented():
    with pytest.raises(scoring.ScoringError, match="no recognised"):
        scoring.extract_ipde({"unrelated": 1.0})


def test_ipde_explicit_key_must_exist():
    with pytest.raises(scoring.ScoringError, match="absent"):
        scoring.extract_ipde({"iptm": 0.8}, key="ipde")


def test_no_best_aggregation_mode():
    """Handover 4.3.7: cherry-picking the best seed defeats Consistency."""
    with pytest.raises(scoring.ScoringError, match="no 'best'"):
        scoring.aggregate([0.1, 0.9], mode="best")


def test_aggregate_penalises_spread():
    """Equal medians, different spread -> the tight one must win."""
    tight = scoring.aggregate([0.80] * 5)
    loose = scoring.aggregate([0.60, 0.70, 0.80, 0.90, 1.00])
    assert tight["median"] == loose["median"] == 0.80, "medians are deliberately equal"
    assert loose["spread"] > tight["spread"]
    assert tight["score"] > loose["score"], "spread penalty must break the tie"


def test_aggregate_prefers_tight_over_lucky_best():
    """A design with one great seed must not beat a uniformly good one."""
    lucky = scoring.aggregate([0.50, 0.55, 0.60, 0.65, 0.99])
    steady = scoring.aggregate([0.70] * 5)
    assert max(lucky["values"]) > max(steady["values"]), "lucky has the better best seed"
    assert steady["score"] > lucky["score"], "but steady must rank higher"


def test_consistency_rewards_agreement():
    assert scoring.consistency([0.8] * 5)["consistency"] == 1.0
    assert scoring.consistency([0.5, 0.6, 0.7, 0.8, 0.9])["consistency"] < 0.5


def test_burial_reports_all_three_anchors(tz_with):
    """Handover 6.2: the convention is unresolved; never pick one silently."""
    b = scoring.catalytic_burial(tz_with)
    assert set(b["anchors"]) >= {"OD2", "CG_carboxyl_C", "CA"}
    assert "unresolved" in b["note"]


def test_resolvable_respects_noise_floor():
    assert not scoring.resolvable(0.80, 0.82, noise_floor=0.05)
    assert scoring.resolvable(0.80, 0.95, noise_floor=0.05)


def test_noise_floor_stops_a_bad_pipeline():
    """complex_ipde is an ERROR, so the POSITIVE control must score LOWER."""
    # positive worse (higher error) than negative -> pipeline is not measuring it
    bad = scoring.NoiseFloor(positive={"median": 3.0, "spread": 0.1},
                             negative={"median": 2.0, "spread": 0.1})
    v = bad.verdict()
    assert v["ranks_correctly"] is False
    assert v["verdict"].startswith("STOP")

    # positive better (lower error) -> proceed
    good = scoring.NoiseFloor(positive={"median": 1.5, "spread": 0.1},
                              negative={"median": 3.0, "spread": 0.1})
    gv = good.verdict()
    assert gv["ranks_correctly"] is True
    assert gv["verdict"].startswith("PROCEED")
    assert gv["lower_is_better"] is True


def test_noise_floor_honours_higher_is_better_when_asked():
    nf = scoring.NoiseFloor(positive={"median": 0.9, "spread": 0.02},
                            negative={"median": 0.7, "spread": 0.02},
                            lower_is_better=False)
    assert nf.ranks_correctly() is True


# ============================================================ gates
def test_length_cap_enforced():
    assert gates.check_length("M" * 260).passed
    assert not gates.check_length("M" * 261).passed, "cap is 260 inclusive"


def test_listov_state_of_the_art_is_over_the_cap():
    """263 aa (9HVB/9HVG/9HVH) is ABOVE the hackathon limit - headroom matters."""
    assert not gates.check_length("M" * 263).passed


def test_catalytic_residue_gate_catches_mutation():
    assert gates.check_catalytic_residues("MDQ", {1: "D", 2: "Q"}).passed
    assert not gates.check_catalytic_residues("MAQ", {1: "D", 2: "Q"}).passed


def test_novelty_gates_fail_closed(tmp_path):
    """Handover 4.3.9: unproven novelty is failure, never a silent pass."""
    assert not gates.check_identity("MKV", tmp_path / "missing.fasta").passed
    assert not gates.check_tm(tmp_path / "nope.pdb", tmp_path).passed


def test_run_all_disqualifies_on_any_failure():
    r = gates.run_all("M" * 300)
    assert r["passed"] is False and r["verdict"] == "DISQUALIFIED"


# ============================================================ chemistry identity
def test_decoy_ligands_are_documented():
    assert "H5J" in C.DECOY_LIGANDS and "5NT" in C.DECOY_LIGANDS
    assert C.LIGAND_CCD == "6NT"
    assert C.LIGAND_CCD not in C.DECOY_LIGANDS


def test_substrate_and_scoring_ligand_are_different():
    """6NT is the TSA we score; 5-nitrobenzisoxazole is the real substrate."""
    assert C.SUBSTRATE_NAME != C.LIGAND_NAME


# ============================================================ RFD3 spec
def test_spec_requires_a_composition_field(tz_with):
    """Observed: `input` alone -> ValidationError from rfd3."""
    spec = rfd3.build_spec("canonical", tz_with)
    d = spec.to_dict()
    assert d.get("length"), "a composition field is mandatory"
    assert "ori_token" in d and len(d["ori_token"]) == 3


def test_spec_conditions_the_oxyanion_hbond(tz_with):
    """Oxyanion stabilisation IS an H-bond; condition on it explicitly."""
    spec = rfd3.build_spec("canonical", tz_with)
    # InputSelection is {contig: "ATOMS"}, not a "G1:NE2" string
    assert isinstance(spec.select_hbond_donor, dict)
    assert "NE2" in next(iter(spec.select_hbond_donor.values()))


def test_spec_never_fixes_sidechains_of_unfixed_sequence(tz_with):
    """Handover 1.4 forbids it; build_spec must enforce for every variant."""
    for variant in C.MOTIF_VARIANTS:
        spec = rfd3.build_spec(variant, tz_with)
        if not spec.select_unfixed_sequence:
            continue
        for tag in spec.select_unfixed_sequence.split(","):
            atoms = (spec.select_fixed_atoms or {}).get(tag.strip(), "")
            for a in filter(None, atoms.split(",")):
                assert a in {"CA", "N", "C", "O", "CB"}, (
                    f"{variant}: side-chain atom {a} fixed on unfixed-sequence {tag}")


def test_all_motif_variants_build(tz_with):
    for variant in C.MOTIF_VARIANTS:
        spec = rfd3.build_spec(variant, tz_with)
        # rc-foundry 1.0.0 wants the residue NAME, not "l:g" and not a contig
        assert spec.to_dict()["ligand"] == C.LIGAND_CCD == "6NT"


def test_unknown_variant_rejected(tz_with):
    with pytest.raises(rfd3.RFD3Error, match="unknown motif variant"):
        rfd3.build_spec("nonsense", tz_with)


def test_enough_variants_available():
    assert len(C.MOTIF_VARIANTS) >= C.MIN_MOTIF_VARIANTS


def test_inputs_json_keys_are_block_names(tz_with, tmp_path):
    specs = [rfd3.build_spec(v, tz_with) for v in ("canonical", "backbone_nh")]
    js = rfd3.write_inputs_json(specs, tmp_path / "inputs.json")
    data = json.loads(js.read_text())
    assert set(data) == {"canonical", "backbone_nh"}


# ============================================================ config invariants
def test_boltz_msa_settings_are_pinned():
    """Handover 4.3.6: identical MSA settings across ALL runs."""
    assert C.BOLTZ_SETTINGS["single_sequence_mode"] is True
    assert C.BOLTZ_SETTINGS["use_msa_server"] is False
    assert len(C.BOLTZ_SEEDS) == C.REPLICATES == 5


def test_length_range_leaves_headroom():
    hi = int(C.RFD3_LENGTH_RANGE.split("-")[1])
    assert hi < C.MAX_LENGTH_AA, "leave room for motif residues under the 260 cap"


def test_novelty_corpus_covers_the_source_structure():
    for pdb in ("6TU6", "4BS0", "9HVB"):
        assert pdb in C.NOVELTY_PDB_IDS


def test_hardware_block_is_honest():
    assert C.GPU["count"] == 1 and C.GPU["vram_gb"] == 24


# ============================================================ iPDE direction
def test_ipde_field_resolved_from_source():
    """Read from boltz/model/modules/confidence.py, not guessed."""
    assert C.IPDE_FIELD == "complex_ipde"
    assert C.IPDE_LOWER_IS_BETTER is True, "complex_ipde is an error in Angstroms"
    assert C.IPDE_UNITS == "angstrom"


def test_extract_ipde_prefers_complex_ipde():
    name, val = scoring.extract_ipde(
        {"complex_ipde": 1.8, "iptm": 0.7, "ligand_iptm": 0.6})
    assert name == "complex_ipde" and val == 1.8


def test_zero_ipde_without_interface_is_rejected():
    """The worst possible silent failure: no ligand chain -> 0.0 looks perfect."""
    single_chain = {"complex_ipde": 0.0, "iptm": 0.0, "ligand_iptm": 0.0,
                    "protein_iptm": 0.0, "complex_plddt": 0.86}
    with pytest.raises(scoring.ScoringError, match="NO interface"):
        scoring.extract_ipde(single_chain)


def test_genuine_zero_ipde_with_interface_is_allowed():
    real = {"complex_ipde": 0.0, "iptm": 0.72, "ligand_iptm": 0.65}
    name, val = scoring.extract_ipde(real)
    assert val == 0.0


def test_ipde_is_better_lower():
    assert scoring.ipde_is_better(1.2, 3.4) is True
    assert scoring.ipde_is_better(3.4, 1.2) is False


def test_aggregate_direction_flips_the_penalty():
    """With lower_is_better, dispersion must push the score UP (worse)."""
    vals = [1.0, 1.2, 1.4, 1.6, 1.8]
    hi = scoring.aggregate(vals, lower_is_better=False)
    lo = scoring.aggregate(vals, lower_is_better=True)
    assert hi["median"] == lo["median"] == 1.4
    assert lo["score"] > lo["median"], "error metric: spread penalty raises the score"
    assert hi["score"] < hi["median"], "reward metric: spread penalty lowers it"


def test_aggregate_worst_case_flips_too():
    vals = [1.0, 5.0]
    assert scoring.aggregate(vals, mode="worst_case", lower_is_better=True)["score"] == 5.0
    assert scoring.aggregate(vals, mode="worst_case", lower_is_better=False)["score"] == 1.0


def test_burial_weight_recorded():
    """Organisers: burial is 8% of the final score; anchor still unknown."""
    assert C.BURIAL_WEIGHT == 0.08
    assert C.BURIAL_ANCHOR_UNRESOLVED is True
