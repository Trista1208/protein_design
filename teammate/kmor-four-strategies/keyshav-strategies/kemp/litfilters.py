"""Literature-derived active-site filters for Kemp eliminase designs.

Sources, all supplied with the project and all for THIS reaction:

  * `filtri_backbone_v2.xlsx` - Privett et al. 2012 PNAS (HG-1/HG-2/HG-3) and
    Alexandrova et al. 2008 JACS (KE07/KE10/KE16). Hard filters H1-H5 and
    soft filters S1-S10.
  * `geometrie_kemp_risso2020_v2.xlsx` - Risso et al. 2020, Chem. Sci. 11,
    6134 (GNCA4 beta-lactamase, Kemp elimination of 5-nitrobenzisoxazole).
    Geometric filters G1-G3 derived from EVB trajectories.
  * `filters hackaton.docx` - the same material plus Bhattacharya et al. 2022.

WHY THIS MODULE EXISTS
----------------------
Our scoring reproduces the ORGANISERS' metric, which is what the competition
ranks on. These filters are different: they are what the literature says
actually correlates with catalytic ACTIVITY. Two designs can tie at 0.000 on
the official metric and differ completely here.

Three findings in those files are worth stating before any of the code:

1. **Rosetta score does NOT correlate with activity.** Risso et al. tested it
   directly on this reaction: "No correlation with activity. Good for
   stability ranking, not for catalysis." The same sheet reports that base
   pKa (PROPKA), reorganization energy and Tm likewise failed as activity
   predictors - Tm is explicitly "a gate, not a predictor".

2. **The predictive geometry is a DISTANCE and an ANGLE**, both measured on
   the donor-hydrogen-acceptor triple, not on the heavy atoms alone:
       D-A distance   WT 2.64 A, range 2.64-4.25, R^2 0.84 vs calc. dG
       D-H...A angle  WT 166.5 deg, range 129.8-167.1, R^2 0.71 vs log(kcat/KM)
   We were not measuring the angle at all.

3. **Combinations beat additivity.** Bhattacharya: "the triple mutant was 3x
   more active than predicted from the three single mutations." That is the
   same effect we measured when F42Y+H103Y did what neither single could.

CAVEATS, because these are other people's thresholds
----------------------------------------------------
The Risso sheet labels its own cutoffs "mie ipotesi dalle letture di Fig. 9:
da tarare" - hypotheses read off a figure, to be calibrated. The paper values
are EVB trajectory AVERAGES, so a single predicted model is not the same
measurement. Treat a FAIL here as a flag for attention, not as a veto, and
never let it override the organisers' own metric.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import config as C

# --------------------------------------------------------------------------- thresholds

#: Risso et al. 2020, sheet "Filtri". soft = earns a point; hard = discard.
RISSO = {
    "G1_distance": {"soft_max": 3.0, "hard_max": 4.0, "unit": "A",
                    "wt": 2.64, "range": (2.64, 4.25),
                    "what": "substrate donor carbon to catalytic Asp oxygen"},
    "G2_angle": {"soft_min": 150.0, "hard_min": 135.0, "unit": "deg",
                 "wt": 166.5, "range": (129.8, 167.1),
                 "what": "donor-hydrogen...acceptor angle"},
}

#: Privett 2012 / Alexandrova 2008, sheet "Filtri".
PRIVETT_HARD = {
    "H1_base_H_distance": {"max": 3.0, "unit": "A",
                           "what": "base to substrate acidic H"},
    "H2_base_rel_sasa": {"max": 10.0, "unit": "%",
                         "what": "catalytic base relative SASA"},
    "H3_sequestering": {"max": 0, "unit": "count",
                        "what": "charged residues within 4 A of the base"},
}
PRIVETT_SOFT = {
    "S2_apolar_fraction": {"min": 0.60, "unit": "fraction",
                           "what": "apolar fraction in the 5-6 A shell"},
    "S4_apo_holo_rmsd": {"max": 1.0, "unit": "A",
                         "what": "catalytic-site RMSD, apo vs holo"},
    "S7_donor_contacts": {"min": 6, "unit": "count",
                          "what": "donor heavy-atom contacts within 4.5 A"},
}

#: Maximum side-chain SASA for a residue, Tien et al. 2013 (theoretical).
#: Used as the denominator for relative SASA. ASP = 193 A^2.
MAX_SASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

#: Alexandrova 2008 H3: a charged side chain this close to the base competes
#: with the substrate. KE16's Lys201-Asp48 salt bridge is the worked example.
CHARGED = {"LYS": ("NZ",), "ARG": ("NH1", "NH2", "NE"),
           "ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"),
           "HIS": ("ND1", "NE2")}
APOLAR = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY", "CYS"}


@dataclass
class LitResult:
    values: dict = field(default_factory=dict)
    hard_fail: list[str] = field(default_factory=list)
    soft_pass: list[str] = field(default_factory=list)
    soft_fail: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not self.hard_fail


# --------------------------------------------------------------------------- geometry

def dh_a_angle(donor_heavy: np.ndarray, h: np.ndarray,
               acceptor: np.ndarray) -> float:
    """The Risso G2 angle: donor-H...acceptor, in degrees.

    Near-linear is good (WT 166.5). This is NOT the organisers' `angle2`,
    which is centred on the ligand acceptor and targets 120 deg - a different
    triple answering a different question. Keep them apart.
    """
    v1, v2 = donor_heavy - h, acceptor - h
    cos = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return float(np.degrees(np.arccos(max(-1.0, min(1.0, cos)))))


def place_polar_h(donor: np.ndarray, neighbours: list[np.ndarray],
                  length: float = 1.01) -> np.ndarray:
    """Idealised H on a planar donor, bisecting its neighbours and pointing out.

    Boltz output carries no hydrogens, and the Risso angle is defined through
    one. Reconstructing it geometrically is standard, but it IS a
    reconstruction: the angle below is therefore a derived quantity, flagged
    as such in the result rather than reported as a measurement.
    """
    v = np.zeros(3)
    for n in neighbours:
        d = donor - n
        v += d / np.linalg.norm(d)
    v /= np.linalg.norm(v)
    return donor + length * v


# --------------------------------------------------------------------------- the filters

def evaluate(structure, *, base_resseq: int, oxy_resseq: int,
             chain: str = "A", ligand: str = C.LIGAND_CCD) -> LitResult:
    """Apply the implementable literature filters to one folded complex.

    ``structure`` is a :class:`theozyme.Theozyme`. Filters needing MD (S1),
    an apo fold (S4) or explicit waters (S3) are not computed here; they are
    listed in ``notes`` so their absence is visible rather than silent.
    """
    from .theozyme import TheozymeError

    r = LitResult()
    prot = [a for a in structure.atoms if a.chain == chain]
    lig = [a for a in structure.atoms if a.resname == ligand]
    if not lig:
        r.hard_fail.append("no ligand in structure")
        return r

    def at(resseq: int, name: str):
        return structure.atom_at(resseq, name, chain=chain).xyz

    # ---- G1 / H1: base carboxylate to the ligand's reactive nitrogen
    try:
        od1, od2 = at(base_resseq, "OD1"), at(base_resseq, "OD2")
        n3 = structure.atom(ligand, "N3").xyz
        d_od1, d_od2 = (float(np.linalg.norm(od1 - n3)),
                        float(np.linalg.norm(od2 - n3)))
        d_base = min(d_od1, d_od2)
        near_o = od2 if d_od2 <= d_od1 else od1
        r.values["G1_distance"] = round(d_base, 3)
        g1 = RISSO["G1_distance"]
        if d_base > g1["hard_max"]:
            r.hard_fail.append(f"G1 distance {d_base:.2f} > {g1['hard_max']} A")
        (r.soft_pass if d_base <= g1["soft_max"] else r.soft_fail).append("G1")
    except TheozymeError as exc:
        r.notes.append(f"G1 unmeasurable: {exc}")
        near_o = None

    # ---- G2: the D-H...A angle, through a reconstructed H on the ligand N3
    try:
        n3 = structure.atom(ligand, "N3").xyz
        n2 = structure.atom(ligand, "N2").xyz
        c3a = structure.atom(ligand, "C3A").xyz
        h = place_polar_h(n3, [n2, c3a])
        if near_o is not None:
            ang = dh_a_angle(n3, h, near_o)
            r.values["G2_angle"] = round(ang, 2)
            r.values["G2_angle_is_reconstructed"] = True
            g2 = RISSO["G2_angle"]
            if ang < g2["hard_min"]:
                r.hard_fail.append(f"G2 angle {ang:.1f} < {g2['hard_min']} deg")
            (r.soft_pass if ang >= g2["soft_min"] else r.soft_fail).append("G2")
            r.notes.append("G2 angle uses an IDEALISED H on ligand N3 - Boltz "
                           "output has no hydrogens (Risso WT 166.5 deg)")
    except TheozymeError as exc:
        r.notes.append(f"G2 unmeasurable: {exc}")

    # ---- H3: sequestering charged residues near the base (Alexandrova KE16)
    try:
        anchor = at(base_resseq, "CG")
        seq_hits = []
        for a in prot:
            if a.resseq == base_resseq or a.resname not in CHARGED:
                continue
            if a.name not in CHARGED[a.resname]:
                continue
            d = float(np.linalg.norm(a.xyz - anchor))
            if d <= 4.0:
                seq_hits.append(f"{a.resname}{a.resseq}.{a.name}@{d:.2f}A")
        r.values["H3_sequestering"] = len(seq_hits)
        r.values["H3_detail"] = seq_hits
        if seq_hits:
            r.hard_fail.append(f"H3 sequestering: {', '.join(seq_hits)}")
    except TheozymeError as exc:
        r.notes.append(f"H3 unmeasurable: {exc}")

    # ---- S2: apolar fraction in the 5-6 A shell around the base
    try:
        anchor = at(base_resseq, "CG")
        shell: dict[int, str] = {}
        for a in prot:
            if a.resseq == base_resseq:
                continue
            if float(np.linalg.norm(a.xyz - anchor)) <= 6.0:
                shell[a.resseq] = a.resname
        if shell:
            frac = sum(1 for v in shell.values() if v in APOLAR) / len(shell)
            r.values["S2_apolar_fraction"] = round(frac, 3)
            r.values["S2_shell_size"] = len(shell)
            (r.soft_pass if frac >= PRIVETT_SOFT["S2_apolar_fraction"]["min"]
             else r.soft_fail).append("S2")
    except TheozymeError as exc:
        r.notes.append(f"S2 unmeasurable: {exc}")

    # ---- S7: donor heavy-atom contacts with the ligand, within 4.5 A
    try:
        donor_atoms = [a for a in prot if a.resseq == oxy_resseq
                       and a.element != "H"]
        lig_xyz = np.array([a.xyz for a in lig])
        n = sum(1 for a in donor_atoms
                if np.min(np.linalg.norm(lig_xyz - a.xyz, axis=1)) <= 4.5)
        r.values["S7_donor_contacts"] = n
        (r.soft_pass if n >= PRIVETT_SOFT["S7_donor_contacts"]["min"]
         else r.soft_fail).append("S7")
    except Exception as exc:                                    # noqa: BLE001
        r.notes.append(f"S7 unmeasurable: {exc}")

    # the pinner - generalised from the Tyr103 observation, measured on 50 designs
    try:
        r.values["pinner"] = gln_pinner(structure, oxy_resseq=oxy_resseq,
                                        chain=chain)
    except Exception as exc:                                    # noqa: BLE001
        r.notes.append(f"pinner unmeasurable: {exc}")

    r.notes.append("NOT computed here (need MD / an apo fold / waters): "
                   "S1 preorganization, S3 water exclusion, S4 apo-holo RMSD, "
                   "S8 pose tolerance")
    return r


def base_relative_sasa(pdb_path, base_resseq: int, chain: str = "A") -> dict:
    """H2: catalytic base relative SASA, computed in the LIGAND-FREE protein.

    Privett's H2 is about whether the base sits in a buried pocket, so the
    ligand is excluded - leaving it in would report the pocket as buried
    simply because something is sitting in it.
    """
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley

    st = PDBParser(QUIET=True).get_structure("x", str(pdb_path))
    model = st[0]
    for ch in list(model):                       # strip ligands and waters
        for res in list(ch):
            if res.id[0] != " ":
                ch.detach_child(res.id)
    ShrakeRupley().compute(model, level="R")
    try:
        res = model[chain][(" ", base_resseq, " ")]
    except KeyError:
        return {"error": f"residue {chain}{base_resseq} not found"}
    total = float(res.sasa)
    rel = 100.0 * total / MAX_SASA.get(res.get_resname(), 193.0)
    return {"resname": res.get_resname(), "sasa_A2": round(total, 2),
            "rel_sasa_pct": round(rel, 2),
            "max_sasa_ref_A2": MAX_SASA.get(res.get_resname(), 193.0),
            "passes_H2": rel <= PRIVETT_HARD["H2_base_rel_sasa"]["max"],
            "threshold_pct": PRIVETT_HARD["H2_base_rel_sasa"]["max"],
            "note": "ligand-free protein; Tien et al. 2013 max-SASA reference"}


# --------------------------------------------------------------------------- the pinner

#: Atoms that can donate or accept a hydrogen bond to the Gln46 amide.
HBOND_ATOMS = {
    "TYR": ("OH",), "TRP": ("NE1",), "SER": ("OG",), "THR": ("OG1",),
    "ASN": ("OD1", "ND2"), "GLN": ("OE1", "NE2"), "HIS": ("ND1", "NE2"),
    "LYS": ("NZ",), "ARG": ("NE", "NH1", "NH2"),
    "ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"),
}

#: MEASURED, 50 bound+measurable designs, 1 seed each:
#:
#:     nearest H-bond donor to the Gln46 amide
#:       IN BAND     (10)  median 2.34 A   range 2.29-2.72
#:       OUT OF BAND (40)  median 3.34 A   range 2.38-3.56
#:
#: All ten in-band designs have Tyr42.OH or Tyr103.OH within 2.72 A. The
#: separation is large but NOT clean - `L_B_anchor_M45F` sits at 2.38 A and is
#: still out of band - so a close pinner is necessary rather than sufficient.
#: Use it to rank and to explain, never as a hard gate.
PINNER_IN_BAND_MAX = 2.75


def gln_pinner(structure, *, oxy_resseq: int = 46, chain: str = "A") -> dict:
    """Which residue hydrogen-bonds the oxyanion Gln, and how close.

    This is the generalised form of the Tyr103 observation. The original
    reading - "the Tyr103 hydroxyl pins the rotamer" - was too narrow: in
    designs carrying F42Y it is Tyr42 that does the pinning, with Tyr103 6.7 A
    away and the dihedral still in band. What predicts the band is that SOME
    H-bond-capable atom reaches the Gln46 amide, not which residue supplies
    it.

    Returned on a 1-seed screen, so it is a cheap proxy for a quantity that
    otherwise costs five folds to establish.
    """
    import numpy as np
    tip = np.array([a.xyz for a in structure.atoms
                    if a.chain == chain and a.resseq == oxy_resseq
                    and a.name in ("OE1", "NE2")])
    if not len(tip):
        return {"distance_A": None, "residue": None,
                "note": f"no Gln{oxy_resseq} amide atoms found"}
    best_d, best_who = 1e9, None
    for a in structure.atoms:
        if a.chain != chain or a.resseq == oxy_resseq:
            continue
        if a.name not in HBOND_ATOMS.get(a.resname, ()):
            continue
        d = float(np.min(np.linalg.norm(tip - a.xyz, axis=1)))
        if d < best_d:
            best_d, best_who = d, f"{a.resname}{a.resseq}.{a.name}"
    if best_who is None:
        return {"distance_A": None, "residue": None,
                "note": "no H-bond-capable atom anywhere near the Gln amide"}
    return {"distance_A": round(best_d, 3), "residue": best_who,
            "within_in_band_range": best_d <= PINNER_IN_BAND_MAX,
            "threshold_A": PINNER_IN_BAND_MAX,
            "note": ("all 10 in-band designs measured so far have a pinner "
                     "<= 2.72 A; necessary but not sufficient")}
