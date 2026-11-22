"""Verified constants for the Kemp eliminase challenge.

Every number here is traceable to HANDOVER_kemp_eliminase.md, which verified them
against primary sources on 21 Nov 2026. Do not edit without updating the handover.
"""

from __future__ import annotations

import os
from pathlib import Path

# --------------------------------------------------------------------------- paths
ROOT = Path(os.environ.get("KEMP_ROOT", "/home/mandalorian/Hackathons/Enzymes-BSSE2026"))
INPUTS = ROOT / "inputs"
THEOZYMES = ROOT / "theozymes-discovered"
OUTPUTS = ROOT / "outputs"
RUNS = OUTPUTS / "runs"
TOOLS = ROOT / "tools"
OPT = Path("/home/mandalorian/opt")

THEOZYME_WITH_ORI = THEOZYMES / "Theozyme_withORI.pdb"
THEOZYME_NO_ORI = THEOZYMES / "Theozyme_noORI.pdb"
#: The input the PIPELINE should scaffold from: identical geometry, but with
#: non-zero B-factors on ASP/GLN/6NT. APNovo infers motif residues from
#: atom_b_factor > 0, and the supplied Gln50 has B=0.00 on all 17 atoms, so it
#: is silently dropped - taking oxy (0.20) and orient_oxy (0.08) with it.
THEOZYME_DESIGN = THEOZYMES / "Theozyme_withORI_bfix.pdb"

# conda env interpreters (see docs/SETUP.md)
ENVS = Path("/home/mandalorian/anaconda3/envs")
PY_FOUNDRY = ENVS / "foundry/bin/python"
PY_AIZYMES = ENVS / "aizymes/bin/python"
PY_BOLTZ = ENVS / "her2_boltz2/bin/python"
PY_MM = ENVS / "mm/bin/python"
PY_QM = ENVS / "qm/bin/python"
BIN_RFD3 = ENVS / "foundry/bin/rfd3"
BIN_RF3 = ENVS / "foundry/bin/rf3"
BIN_BOLTZ = ENVS / "her2_boltz2/bin/boltz"
LIGANDMPNN_DIR = TOOLS / "LigandMPNN"
FIELDTOOLS = TOOLS / "FieldTools.py"

# --------------------------------------------------------------------------- chemistry
# Handover 4.3.5: 6-nitrobenzotriazole is the SCORING ligand (transition-state analogue).
LIGAND_CCD = "6NT"
LIGAND_FORMULA = "C6H4N4O2"
LIGAND_NAME = "6-nitrobenzotriazole"
# Decoys that must never be silently substituted (handover 4.3.5 / 1.5):
DECOY_LIGANDS = {
    "H5J": "5-nitro-1,2-benzoxazole (the real SUBSTRATE, from 7EPO) - not the TSA",
    "5NT": "5-nitrobenzotriazole (the 3NYD ligand) - wrong nitro position",
}
SUBSTRATE_NAME = "5-nitrobenzisoxazole"  # the reaction substrate (slides p.4)

# --------------------------------------------------------------------------- reference geometry
# Handover section 2, measured from the supplied theozyme; round-trips to deviation 0.0.
REF_GEOMETRY = {
    "asp127_od2__6nt_n3": 2.48,        # Angstrom - the catalytic base contact
    "asp127_od2__h3": 1.53,
    "gln50_ne2__6nt_n1": 3.33,         # CONFIRMED against 4BS0/7K4Z - see below
    "asp127_ca__gln50_ca": 10.64,      # the real scaffolding constraint
    "od2_n3_n1_angle_deg": 158.4,
    "od1__n3_nonengaging": 4.56,
}
# RESOLVED (handover open question 1): the organisers' modelled Gln50 DOES match
# the crystallographic Gln50 in 4BS0/7K4Z. The B-factor 0.00 signature means it
# was modelled in, but the model agrees with the crystal structures, so the
# 3.33 A contact is trustworthy after all and carries normal weight.
REF_GEOMETRY_MODELLED: set[str] = set()
GLN50_PROVENANCE = (
    "modelled in (B-factor 0.00, no ANISOU) but CONFIRMED to match the "
    "crystallographic Gln50 in 4BS0/7K4Z - treat as trustworthy"
)

# ORI pseudo-atom, handover 1.2. Frame-specific - NEVER mix frames.
ORI_IN_WITHORI_FRAME = (25.568, 0.048, -23.462)
ORI_IN_NOORI_FRAME = (-33.913, 30.683, -33.510)
# Assertion threshold from handover 1.2 / 4.3.2: a mixed frame lands ORI ~70.5 A away.
ORI_LIGAND_MAX_DIST = 12.0

# --------------------------------------------------------------------------- challenge rules
MAX_LENGTH_AA = 260            # slides p.7, motif residues INCLUDED
SEQ_IDENTITY_DISQUALIFY = 0.70  # >70% identity to published enzymes disqualifies
TM_SCORE_DISQUALIFY = 0.70      # TM-score >0.7 to published structures disqualifies
REPLICATES = 5                  # 5 Boltz seeds per sequence
BURIAL_RADIUS_A = 12.0          # catalytic base burial: CA atoms within 12 A
N_LAB_SLOTS = 5

#: Minimum mean complex pLDDT for a design to be scored at all.
#:
#: A design that does not fold confidently is UNMEASURABLE, not bad
#: (LESSONS 7.1): at pLDDT 0.41-0.47 the same sequence folded twice gave d_base
#: 2.56 and 13.64 A, so its geometry terms are noise and ranking on them ranks
#: noise. denovoKemp (Listov et al., Nature 2025, Methods) discards pLDDT < 90
#: for its final series; the team set 0.95 here, which is stricter.
#:
#: Calibration on this project's own data:
#:   rfd4_design_045, 5 Boltz seeds : 0.963 - 0.967   PASS
#:   teammate designs               : 0.80  - 0.90    fail at 0.95
#:   local RFD3 + LigandMPNN designs: 0.41  - 0.47    fail by a mile
#: So this gate rejects most of what the 21 Nov campaign produced. That is the
#: intended effect: those designs were never measurable.
PLDDT_MIN = 0.95

# Official metric weight for catalytic base burial, per the organisers: 8%.
# The anchor convention (carboxylate carbon / OD2 / CA) is still UNKNOWN, and
# whether the base's own CA counts is unknown. At 8% of the total the ambiguity
# is low-risk, so we report all three anchors and do not block on it.
BURIAL_WEIGHT = 0.08
BURIAL_ANCHOR_UNRESOLVED = True

# ----- Boltz-2 iPDE: RESOLVED from the installed source, not guessed ---------
# boltz/model/modules/confidence.py:451-456
#   token_interface_pair_mask = token_pair_mask * (asym_id_i != asym_id_j)
#   complex_ipde = (pde * mask).sum() / mask.sum()
# `pde` is an expected distance error in ANGSTROMS (bins 0-32), averaged over
# CROSS-CHAIN token pairs only.
IPDE_FIELD = "complex_ipde"
IPDE_LOWER_IS_BETTER = True      # it is an ERROR, unlike plddt/iptm
IPDE_UNITS = "angstrom"
# DANGER: with a single chain the interface mask sums to zero and boltz returns
# complex_ipde == 0.0, which looks like a perfect score. A design folded without
# the ligand as its own chain will therefore score spuriously well. Always check
# that an interface actually existed before trusting the number.
IPDE_NO_INTERFACE_SENTINEL = 0.0

# Handover 3.3.1: state of the art (Listov) is 263 aa, ABOVE the cap. Leave headroom.
RFD3_LENGTH_RANGE = "190-250"

# --------------------------------------------------------------------------- novelty corpus
# Handover 4.3.8: screen against a real corpus, not just 4BS0/6TU6.
NOVELTY_PDB_IDS = [
    "6TU6",  # HG3.17 Q50M + 6NT - SOURCE of the supplied theozyme
    "4BS0",  # HG3.17 E47N/N300D + 6NT
    "7K4Z", "7K4Q",  # HG3.17 / HG3 + 6NT (Kern lab)
    "3NYD",  # HG-2 + 5-nitrobenzotriazole (different ligand)
    "7EPO", "7EPN",  # M. smegmatis KSI + H5J / apo
    "9HVB", "9HVG", "9HVH",  # Listov Des27 family, 263 aa
    "1A53",  # indole-3-glycerol phosphate synthase scaffold (1A53-2 lineage)
]
NOVELTY_NOTE = (
    "KE07/KE70 (Rothlisberger 2008) and the HG series must be covered too; "
    "fetch by sequence from UniProt/PDB if no structure id applies."
)

# --------------------------------------------------------------------------- controls
# Handover 3.2: the control pair the organisers accidentally handed over.
CONTROL_POSITIVE = {
    "name": "HG3.17",
    "pdb": "4BS0",
    "note": "Gln50 present; most active Kemp eliminase of its lineage. MUST rank above the negative.",
}
CONTROL_NEGATIVE = {
    "name": "HG3.17-Q50M",
    "pdb": "6TU6",
    "note": "oxyanion-stabiliser knockout; experimentally worse. If the pipeline cannot "
            "rank the positive above this, Geometry/Consistency are not measuring "
            "oxyanion stabilisation and every later number is noise.",
}

# --------------------------------------------------------------------------- motif variants
# Handover 4.2: scaffold at least three; five correlated designs waste all five lab slots.
MOTIF_VARIANTS = {
    "canonical":      "Asp base + Gln donor, exactly as supplied",
    "crystal_gln":    "same, but the real Gln50 rotamer from 4BS0/7K4Z if it differs",
    "backbone_nh":    "Asp base, oxyanion stabilised by a backbone amide instead of a side chain",
    "alt_donor":      "Ser/Thr/Tyr hydroxyl in place of Gln",
    "base_variant":   "Glu instead of Asp (longer reach, different rotamer space)",
    # The supplied Gln BACKBONE limits d_oxy to >=3.27 A (exhaustive rotamer
    # scan: best 3.265 vs supplied 3.331). The organisers' cst target is 2.40 A,
    # so fixing that backbone locks a mediocre 0.20-weight term into every
    # design. This variant frees the donor entirely - position AND identity -
    # and conditions only the H-bond onto the ligand, letting RFD3 find a
    # scaffold where a tight contact is actually reachable.
    "oxy_free":       "Asp base fixed; oxyanion donor free to be placed anywhere",
    # Roles sanctioned by the scorer's own KE_MOTIF_ROLES that the handover
    # variant list omitted:
    "polarizer":      "canonical plus a SECOND carboxylate (D/E) polarizer role",
    "aromatic":       "canonical plus an aromatic (F/W/Y) pi-stacking residue",
}
MIN_MOTIF_VARIANTS = 3

# --------------------------------------------------------------------------- Boltz settings
# Handover 4.3.6: identical MSA settings across ALL runs or confidences are incomparable.
# De novo designs have no meaningful MSA -> single-sequence mode, fixed forever.
BOLTZ_SETTINGS = {
    "single_sequence_mode": True,
    "use_msa_server": False,
    "recycling_steps": 3,
    "diffusion_samples": 1,     # per seed; replicates come from distinct seeds
}
BOLTZ_SEEDS = [1, 2, 3, 4, 5]

# --------------------------------------------------------------------------- hardware
# Measured on this box (see docs/CONCURRENCY.md). One GPU only.
GPU = {"name": "NVIDIA GeForce RTX 3090", "vram_gb": 24, "count": 1,
       "driver": "580.178.04", "cuda_cap": "13.0"}
HOST = {"ram_gb": 92, "cores": 16}
