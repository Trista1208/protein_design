"""The organisers' official scoring function.

Weights as supplied (they sum to exactly 1.00):

    base        0.25   drives d_base
    ipde        0.20   drives ipde
    oxy         0.20   drives d_oxy
    orient_base 0.12   drives min(angle2, dihedral12) at the base
    orient_oxy  0.08   drives min(angle2, dihedral12) at the oxyanion
    bury        0.08   drives ca_neighbors_mean
    cons        0.07   drives sd_base

Three consequences worth stating plainly, because they redirect the whole effort:

1. GEOMETRY IS 65% OF THE SCORE.
   base + oxy + orient_base + orient_oxy = 0.25+0.20+0.12+0.08 = 0.65.
   iPDE, the thing folding models are usually optimised for, is only 0.20.
   Designs should be filtered on active-site geometry FIRST and folding
   confidence second - the opposite of the usual de novo instinct.

2. OXYANION STABILISATION IS 0.28 OF THE SCORE (oxy 0.20 + orient_oxy 0.08).
   That is more than iPDE. Combined with the confirmation that Gln50 matches
   4BS0/7K4Z, the oxyanion geometry is a first-class objective, not a nicety.

3. CONSISTENCY IS sd_base, NOT the spread of iPDE.
   cons (0.07) is driven by the standard deviation of d_base ACROSS REPLICATES.
   So geometry must be measured on EVERY one of the 5 folded models, not just
   on one. A pipeline that folds 5 seeds but only measures geometry once cannot
   compute this term at all.

UNRESOLVED, and flagged rather than invented:
  * The exact atom quadruples behind `angle2` and `dihedral12`. The names imply
    an angle1/angle2 pair plus a dihedral spanning atoms 1-2 of an H-bond.
    ORIENT_ATOMS below encodes a defensible reading; swap it when confirmed.
  * Whether `min(...)` is over raw degrees or over normalised sub-scores. Taking
    a min of an angle and a dihedral in raw degrees is dimensionally odd, so we
    min over NORMALISED sub-scores, which is the conservative "both must be
    good" reading.
  * The burial anchor (handover open question 2) - still unknown, 8% weight.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from statistics import mean, stdev
from typing import Sequence

import numpy as np

from . import config as C
from .theozyme import Theozyme, TheozymeError, load

# --------------------------------------------------------------------------- weights
WEIGHTS = {
    "base": 0.25,
    "ipde": 0.20,
    "oxy": 0.20,
    "orient_base": 0.12,
    "orient_oxy": 0.08,
    "bury": 0.08,
    "cons": 0.07,
}
assert abs(sum(WEIGHTS.values()) - 1.0) < 1e-9, "official weights must sum to 1"

GEOMETRY_SHARE = WEIGHTS["base"] + WEIGHTS["oxy"] + WEIGHTS["orient_base"] + WEIGHTS["orient_oxy"]
OXYANION_SHARE = WEIGHTS["oxy"] + WEIGHTS["orient_oxy"]

# --------------------------------------------------------------------------- targets
#: Lauko scoring windows (hbonds.py:142) - the SUPPORT of each subscore.
#: Outside the window s = 0; a geometric mean then zeroes the whole design.
LAUKO = {
    "dist": (2.0, 3.6),
    "angle_sp2": (80.0, 160.0),
    "dihedral_sp2_180": (130.0, 230.0),
    "dihedral_sp2_0": (310.0, 50.0),      # wraps through 360/0
    "angle_sp3": (69.5, 149.5),
}

#: Target values measured from the CRYSTAL STRUCTURES, not from a constraint file.
#:
#: BRIEF 0.2 retired the earlier 2.40 A oxyanion target: the cst's 2.4 A is for
#: NH2O...O5 (Gln amide to the real TS OXYGEN), whereas we measure NE2...N1
#: (amide to a ring NITROGEN of the TSA). Different acceptor, different chemical
#: state, and 2.4 A is unusually short for N-H...N.
#:
#: Measured here from the two structures with a real refined Gln50:
#:            Gln50 NE2...N1      Asp127 OD2...N3
#:   4BS0 A       3.102                2.532
#:   7K4Z A       3.100                2.514
#:   supplied     3.331                2.482
#: So the crystallographic oxyanion contact is 3.10 A. Using 2.40 was the sole
#: reason `oxy` scored 0.224 on a geometrically healthy motif.
#:
#: BRIEF 0.1 also RETRACTED the "re-target 120 deg" advice: all eight
#: angle2/dihedral2 values pass both RDKit reference orderings with wide margin
#: (angles 121-125 deg, dihedrals 163-198 deg). The motif's angular geometry is
#: healthy; 120 deg stays only as the window centre, not as a correction target.
TARGETS = {
    "d_base": {"target": 2.52, "window": LAUKO["dist"],
               "provenance": "mean of 4BS0 2.532 / 7K4Z 2.514"},
    "d_oxy":  {"target": 3.10, "window": LAUKO["dist"],
               "provenance": "mean of 4BS0 3.102 / 7K4Z 3.100"},
    "angle":  {"target": 120.0, "window": LAUKO["angle_sp2"],
               "provenance": "5TS cst; theozyme measures 121-125, mid-window"},
}

#: Taper half-widths. The ramp shape between window edges is UNRESOLVED
#: (BRIEF 5.2 lists it as a question for the organisers). Fitted to the only
#: two data points in SCORING.md's worked example:
#:   angle 120 deg          -> s = 1.000   (window 80-160, half-width 40)
#:   dihedral 210 deg       -> s = 0.500   => half-width 60 about an ideal 180
#:   dihedral 301.5 deg     -> s = 0.000   (outside both sp2 bands)
ANGLE_HALF_WIDTH = 40.0
DIHEDRAL_IDEAL = 180.0
DIHEDRAL_HALF_WIDTH = 60.0

#: Because the ramp is unresolved we score BOTH models and report both, so the
#: sensitivity is visible rather than hidden in a constant.
#:   "taper" - linear decay from target to window edge (conservative)
#:   "flat"  - 1.0 anywhere inside the Lauko window (thresholds only)
RAMP_MODE = "taper"

RAMP_SHAPE = "see RAMP_MODE; both taper and flat are reported"


def _window_score(value: float, target: float, lo: float, hi: float) -> float:
    """1.0 at ``target``, linear decay to 0.0 at whichever window edge is nearer."""
    if value is None:
        return 0.0
    if not (lo <= value <= hi):
        return 0.0
    span = (hi - target) if value >= target else (target - lo)
    if span <= 0:
        return 1.0 if abs(value - target) < 1e-9 else 0.0
    return max(0.0, 1.0 - abs(value - target) / span)


def _flat_score(value: float, lo: float, hi: float) -> float:
    """Threshold-only model: 1.0 anywhere inside the window, else 0."""
    return 0.0 if value is None else (1.0 if lo <= value <= hi else 0.0)


def _dihedral_score(value: float) -> float:
    """Best score across the two sp2 dihedral bands, on [0, 360)."""
    if value is None:
        return 0.0
    v = value % 360.0
    best = 0.0
    lo, hi = LAUKO["dihedral_sp2_180"]
    if lo <= v <= hi:
        # fitted: 210 deg -> 0.500 about an ideal of 180 => half-width 60
        best = max(best, max(0.0, 1.0 - abs(v - DIHEDRAL_IDEAL) / DIHEDRAL_HALF_WIDTH))
    lo0, hi0 = LAUKO["dihedral_sp2_0"]          # 310 -> 50, wrapping
    if v >= lo0 or v <= hi0:
        d = min(abs(v - 360.0), abs(v - 0.0))
        best = max(best, max(0.0, 1.0 - d / DIHEDRAL_HALF_WIDTH))
    return best


def _dihedral_flat(value: float) -> float:
    if value is None:
        return 0.0
    v = value % 360.0
    lo, hi = LAUKO["dihedral_sp2_180"]
    lo0, hi0 = LAUKO["dihedral_sp2_0"]
    return 1.0 if (lo <= v <= hi) or (v >= lo0 or v <= hi0) else 0.0

#: INTERPRETED atom selections for the orientation terms. Swap when confirmed.
#: angle2  = angle centred on the ligand acceptor atom
#: dihedral12 = dihedral across the contact, donor-side atom -> acceptor-side atom
ORIENT_ATOMS = {
    "base": {
        "donor_heavy": ("ASP", "OD2"),
        "donor_root": ("ASP", "CG"),
        "acceptor": (C.LIGAND_CCD, "N3"),
        "acceptor_root": (C.LIGAND_CCD, "N2"),
        # measured on Theozyme_withORI.pdb (OD2...N3-N2). NOT 158.4, which is
        # the OD2-N3-N1 angle over a different atom triple.
        "angle2_ref": 124.637,
        "angle2_tol": 20.0,
        "dihedral_ref": 31.867,
        "dihedral_tol": 40.0,
    },
    "oxy": {
        "donor_heavy": ("GLN", "NE2"),
        "donor_root": ("GLN", "CD"),
        "acceptor": (C.LIGAND_CCD, "N1"),
        "acceptor_root": (C.LIGAND_CCD, "N2"),
        # measured on Theozyme_withORI.pdb (NE2...N1-N2)
        "angle2_ref": 124.829,
        "angle2_tol": 25.0,
        "dihedral_ref": 50.902,
        "dihedral_tol": 40.0,
    },
}

#: burial = cat_base_tip_ca_neighbors_mean (kemp_metrics.py:438), VERIFIED:
#:   * anchor is the TIP atom - the motif residue atom closest to the reactive
#:     ligand atom (for Asp, whichever of OD1/OD2 is nearer N3)
#:   * count protein CA within EACH of 6, 8 and 12 A and take the MEAN of the
#:     three. NOT a single 12 A count.
#:   * NO exclusion of the base's own CA.
#: Published reference threshold: cat_base_tip_ca_neighbors_mean >= 15.
CA_DIST_CUTOFFS = (6.0, 8.0, 12.0)
BURY_TARGET = {"good": 15.0, "poor": 3.0}

#: consistency: sd of d_base across replicates, in Angstroms. 0 is perfect.
CONS_TOLERANCE_A = 0.50


# --------------------------------------------------------------------------- helpers
def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    v1, v2 = a - b, c - b
    cos = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return float(np.degrees(np.arccos(max(-1.0, min(1.0, cos)))))


def _dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(b1n, v), w))
    # APNovo remaps [-180,180] -> [0,360) (hbonds.py:485-500). Taking abs()
    # instead collapses 310 deg onto 50 deg and silently scores the wrong band.
    return float(np.degrees(np.arctan2(y, x)) % 360.0)


def _circular_mean(values: Sequence[float]) -> float:
    """Mean of angles on [0, 360), via the unit-vector sum.

    A linear mean of a circular quantity is not merely imprecise, it can land
    in the opposite half of the circle: mean(350, 10) = 180, the antipode of
    the correct 0. For the sp2 dihedral bands that is the difference between
    the 310-50 band and the 230-310 dead zone, so it silently zeroes an
    orientation term that is in fact good. See combine_replicates().
    """
    if not values:
        raise ValueError("no values to average")
    s = sum(math.sin(math.radians(v)) for v in values)
    c = sum(math.cos(math.radians(v)) for v in values)
    if abs(s) < 1e-12 and abs(c) < 1e-12:
        # antipodal values cancel: the mean direction is genuinely undefined,
        # so fall back to the first value rather than inventing 0.0
        return float(values[0]) % 360.0
    return math.degrees(math.atan2(s, c)) % 360.0


def _tol_score(value: float, ref: float, tol: float) -> float:
    """1.0 at the reference, decaying linearly to 0 at 3x tolerance."""
    return max(0.0, 1.0 - abs(value - ref) / (3.0 * tol))


def _ramp(value: float, poor: float, good: float) -> float:
    """0 at or below ``poor``, 1 at or above ``good``."""
    if good == poor:
        return 0.0
    return max(0.0, min(1.0, (value - poor) / (good - poor)))


# --------------------------------------------------------------------------- raw metrics
@dataclass
class RawMetrics:
    """The organisers' quantities, measured. ``None`` means not measurable."""
    d_base: float | None = None
    d_oxy: float | None = None
    angle2_base: float | None = None
    dihedral12_base: float | None = None
    angle2_oxy: float | None = None
    dihedral12_oxy: float | None = None
    ca_neighbors_mean: float | None = None
    ipde: float | None = None
    sd_base: float | None = None
    n_replicates: int = 0
    notes: list[str] = field(default_factory=list)


@dataclass
class CatalyticSite:
    """Where the motif actually landed in the designed chain.

    RFD3 reports this in its output JSON as ``diffused_index_map``, e.g.
    ``{"G2": "A95", "G1": "A39"}`` - theozyme Asp (G2) -> design position A95,
    Gln (G1) -> A39. Without it, geometry cannot be measured on a real design:
    a 240-residue protein has ~15 aspartates and a residue-TYPE lookup is
    ambiguous.
    """
    base_resseq: int
    oxy_resseq: int
    chain: str | None = None
    ligand_resname: str = C.LIGAND_CCD

    @classmethod
    def from_diffused_index_map(cls, dim: dict[str, str], *,
                                base_tag: str, oxy_tag: str) -> "CatalyticSite":
        def parse(v: str) -> tuple[str, int]:
            return v[0], int(v[1:])
        bc, bi = parse(dim[base_tag])
        oc, oi = parse(dim[oxy_tag])
        if bc != oc:
            raise ValueError(f"motif split across chains {bc}/{oc}; unsupported")
        return cls(base_resseq=bi, oxy_resseq=oi, chain=bc)


def measure_raw(structure: Theozyme | str, *, burial_anchor: str = "OD2",
                site: CatalyticSite | None = None) -> RawMetrics:
    """Measure the per-structure geometry terms from one folded model.

    ``site`` is REQUIRED for real designs. Without it we fall back to
    residue-type lookup, which only works on the 3-residue theozyme itself.
    """
    tz = structure if isinstance(structure, Theozyme) else load(structure)
    m = RawMetrics()

    if site is None:
        def xyz(resname: str, atom: str):
            return tz.atom(resname, atom).xyz
        m.notes.append("no CatalyticSite given - residue-TYPE lookup; valid only "
                       "for the bare theozyme, NOT for designs")
    else:
        # verify the design actually kept the catalytic identities before
        # measuring geometry that would otherwise be meaningless
        got_base = tz.resname_at(site.base_resseq, chain=site.chain)
        got_oxy = tz.resname_at(site.oxy_resseq, chain=site.chain)
        if got_base not in ("ASP", "GLU"):
            m.notes.append(f"catalytic base at {site.chain}{site.base_resseq} is "
                           f"{got_base}, not ASP/GLU - motif was not preserved")
        if got_oxy not in ("GLN", "ASN", "SER", "THR", "TYR"):
            m.notes.append(f"oxyanion donor at {site.chain}{site.oxy_resseq} is "
                           f"{got_oxy} - unexpected donor identity")

        def xyz(resname: str, atom: str):
            if resname == C.LIGAND_CCD:
                return tz.atom(resname, atom).xyz      # ligand is unique
            rs = site.base_resseq if resname in ("ASP", "GLU") else site.oxy_resseq
            return tz.atom_at(rs, atom, chain=site.chain).xyz

    # distances
    try:
        m.d_base = float(np.linalg.norm(xyz("ASP", "OD2") - xyz(C.LIGAND_CCD, "N3")))
    except TheozymeError as exc:
        m.notes.append(f"d_base unmeasurable: {exc}")
    try:
        m.d_oxy = float(np.linalg.norm(xyz("GLN", "NE2") - xyz(C.LIGAND_CCD, "N1")))
    except TheozymeError as exc:
        m.notes.append(f"d_oxy unmeasurable: {exc}")

    # orientation terms  (loop var is `which`, NOT `site` - that would shadow
    # the CatalyticSite parameter used by xyz())
    for which, cfg in ORIENT_ATOMS.items():
        try:
            dh = xyz(*cfg["donor_heavy"])
            dr = xyz(*cfg["donor_root"])
            ac = xyz(*cfg["acceptor"])
            ar = xyz(*cfg["acceptor_root"])
            ang = _angle(dh, ac, ar)          # centred on the ligand acceptor
            dih = _dihedral(dr, dh, ac, ar)
            if which == "base":
                m.angle2_base, m.dihedral12_base = ang, dih
            else:
                m.angle2_oxy, m.dihedral12_oxy = ang, dih
        except TheozymeError as exc:
            m.notes.append(f"orient_{which} unmeasurable: {exc}")

    # burial
    # burial: TIP atom -> mean CA count over (6, 8, 12) A, no own-CA exclusion
    try:
        reactive = xyz(C.LIGAND_CCD, "N3")
        tip_candidates = ("OD1", "OD2") if site is None else ("OD1", "OD2", "OE1", "OE2")
        tip, tip_name, tip_d = None, None, 1e9
        for nm in tip_candidates:
            try:
                cand = xyz("ASP", nm)
            except TheozymeError:
                continue
            dd = float(np.linalg.norm(cand - reactive))
            if dd < tip_d:
                tip, tip_name, tip_d = cand, nm, dd
        if tip is None:
            raise TheozymeError("no carboxylate oxygen found for the tip atom")
        counts = []
        cas = [a.xyz for a in tz.atoms if a.name == "CA"]
        for cutoff in CA_DIST_CUTOFFS:
            counts.append(sum(1 for c in cas
                              if float(np.linalg.norm(c - tip)) <= cutoff))
        m.ca_neighbors_mean = float(sum(counts) / len(counts))
        m.notes.append(
            f"burial: tip={tip_name} ({tip_d:.2f} A from N3), CA counts "
            f"{counts} at {CA_DIST_CUTOFFS} A, mean={m.ca_neighbors_mean:.2f}; "
            "own CA NOT excluded (matches kemp_metrics.py:438)")
    except Exception as exc:                  # noqa: BLE001
        m.notes.append(f"burial unmeasurable: {exc}")
    return m


def combine_replicates(per_replicate: Sequence[RawMetrics],
                       ipde_values: Sequence[float] | None = None) -> RawMetrics:
    """Fold N replicate measurements into the scored quantities.

    ``sd_base`` is the standard deviation of d_base across replicates - the
    official Consistency term. This is why geometry must be measured on every
    replicate, not once.
    """
    if not per_replicate:
        raise ValueError("no replicate measurements supplied")

    out = RawMetrics(n_replicates=len(per_replicate))
    for fld in ("d_base", "d_oxy", "angle2_base", "angle2_oxy",
                "ca_neighbors_mean"):
        vals = [getattr(r, fld) for r in per_replicate if getattr(r, fld) is not None]
        setattr(out, fld, float(mean(vals)) if vals else None)

    # Dihedrals are CIRCULAR and must not be averaged linearly.
    #
    # Measured on rfd4_design_045, 5 Boltz seeds, base dihedral:
    #     354.8  350.2  349.6  83.4  355.6
    # Four of five sit near 0/360 - deep inside the 310-50 sp2 band and close
    # to ideal. A linear mean returns 298.7 deg, which is in the 230-310 DEAD
    # ZONE and scores 0.000; the circular mean returns 6.7 deg and scores
    # 0.889. So the linear mean was zeroing `orient_base` (weight 0.12) on a
    # design whose base dihedral is actually fine, purely as an averaging
    # artefact. Averaging 350 and 10 must give 0, never 180.
    for fld in ("dihedral12_base", "dihedral12_oxy"):
        vals = [getattr(r, fld) for r in per_replicate if getattr(r, fld) is not None]
        setattr(out, fld, _circular_mean(vals) if vals else None)

    bases = [r.d_base for r in per_replicate if r.d_base is not None]
    if len(bases) >= 2:
        out.sd_base = float(stdev(bases))
    elif len(bases) == 1:
        # Leave sd_base as None. With one seed the SD is 0.0 by construction,
        # which would hand every prefiltered design a free 0.07. The organisers
        # note ICC 0.47-0.59 and say that argues for MORE SEEDS, not a lower
        # weight - so a 1-seed run must not score this term at all.
        out.notes.append("sd_base NOT scored: one replicate only. The Consistency "
                         "term requires >=2 seeds; a single-seed SD of 0.0 is an "
                         "artefact, not a measurement.")

    if ipde_values:
        finite = [float(v) for v in ipde_values if v is not None and not math.isnan(float(v))]
        if finite:
            out.ipde = float(mean(finite))
    for r in per_replicate:
        out.notes.extend(r.notes)
    return out


# --------------------------------------------------------------------------- score
def score(raw: RawMetrics, *, ipde_scale: float = 8.0,
          strict: bool = True) -> dict:
    """Apply the official weights. Returns subscores, weighted score and gaps.

    ``strict`` (handover 4.3.9, fail closed) refuses to emit a total when any
    weighted term is unmeasurable, instead of quietly scoring a partial design
    against a full-weight scale.
    """
    subs: dict[str, float | None] = {}
    detail: dict[str, dict] = {}

    # --- base (0.25)
    if raw.d_base is not None:
        t = TARGETS["d_base"]
        subs["base"] = _window_score(raw.d_base, t["target"], *t["window"])
        detail["base"] = {"d_base_A": round(raw.d_base, 3), **t}
    else:
        subs["base"] = None

    # --- oxy (0.20)
    if raw.d_oxy is not None:
        t = TARGETS["d_oxy"]
        subs["oxy"] = _window_score(raw.d_oxy, t["target"], *t["window"])
        detail["oxy"] = {"d_oxy_A": round(raw.d_oxy, 3), **t}
    else:
        subs["oxy"] = None

    # --- ipde (0.20) - complex_ipde is an ERROR in Angstroms, lower is better
    if raw.ipde is not None:
        subs["ipde"] = max(0.0, 1.0 - raw.ipde / ipde_scale)
        detail["ipde"] = {"ipde_A": round(raw.ipde, 4), "scale_A": ipde_scale,
                          "lower_is_better": True,
                          "field": C.IPDE_FIELD}
    else:
        subs["ipde"] = None

    # --- orientation (0.12 / 0.08): min over NORMALISED subscores
    for site, key in (("base", "orient_base"), ("oxy", "orient_oxy")):
        cfg = ORIENT_ATOMS[site]
        ang = raw.angle2_base if site == "base" else raw.angle2_oxy
        dih = raw.dihedral12_base if site == "base" else raw.dihedral12_oxy
        if ang is None or dih is None or cfg["angle2_ref"] is None:
            subs[key] = None
            continue
        at = TARGETS["angle"]
        s_ang = _window_score(ang, at["target"], *at["window"])
        s_dih = _dihedral_score(dih)
        subs[key] = min(s_ang, s_dih)
        # emitted separately so a low orientation score can be attributed to the
        # angle or to the torsion (SCORING.md emits s_orient_*_ang / _dih)
        detail[key] = {"angle2_deg": round(ang, 2), "dihedral2_deg": round(dih, 2),
                       f"s_{key}_ang": round(s_ang, 4),
                       f"s_{key}_dih": round(s_dih, 4),
                       "rule": "min(angle2, dihedral2) - both must be good",
                       "angle_target": at["target"], "angle_window": at["window"],
                       "ramp": RAMP_SHAPE}

    # --- bury (0.08)
    if raw.ca_neighbors_mean is not None:
        subs["bury"] = _ramp(raw.ca_neighbors_mean, BURY_TARGET["poor"], BURY_TARGET["good"])
        detail["bury"] = {"ca_neighbors_mean": round(raw.ca_neighbors_mean, 2),
                          **BURY_TARGET, "cutoffs_A": CA_DIST_CUTOFFS,
                          "note": "weight 0.08 - meet it and stop optimising it"}
    else:
        subs["bury"] = None

    # --- cons (0.07) - sd of d_base across replicates, lower is better
    if raw.sd_base is not None:
        subs["cons"] = max(0.0, 1.0 - raw.sd_base / CONS_TOLERANCE_A)
        detail["cons"] = {"sd_base_A": round(raw.sd_base, 4),
                          "tolerance_A": CONS_TOLERANCE_A,
                          "n_replicates": raw.n_replicates,
                          "lower_is_better": True}
    else:
        subs["cons"] = None

    missing = [k for k, v in subs.items() if v is None]
    covered = sum(WEIGHTS[k] for k, v in subs.items() if v is not None)

    # ---------------------------------------------------------------- the scalar
    # S = prod_i s_i ** w_i  =  exp(sum_i w_i * ln s_i)
    # A WEIGHTED GEOMETRIC MEAN. There is no compensation between terms: a
    # design that is excellent on six components and 0.000 on the seventh
    # scores exactly 0.000. The filter must therefore be "non-zero everywhere"
    # BEFORE it is ever "high on average" - the opposite of the usual
    # optimise-the-weighted-sum instinct.
    present = {k: v for k, v in subs.items() if v is not None}
    zeroed = [k for k, v in present.items() if v <= 0.0]
    if zeroed:
        total = 0.0
    else:
        total = math.exp(sum(WEIGHTS[k] * math.log(v) for k, v in present.items()))

    # arithmetic mean kept ONLY as a diagnostic: it shows how much the geometric
    # mean is punishing a single weak term, and gives a usable ranking signal
    # among designs that all score a hard 0.
    arithmetic = sum(WEIGHTS[k] * v for k, v in present.items())
    # shadow score: same geometric form with a small floor, so designs that are
    # "bad but not impossible" remain rankable internally. NEVER report this as
    # the official score.
    floored = {k: max(v, 1e-3) for k, v in present.items()}
    shadow = math.exp(sum(WEIGHTS[k] * math.log(v) for k, v in floored.items()))

    # ---------------------------------------------------------------- triage
    # ln S = sum_i w_i ln s_i. Rank in LOG SPACE: any arithmetic combination
    # reorders candidates relative to the official score.
    # The limiting term is argmin_i (w_i * ln s_i) - the component dragging the
    # design down hardest. Marginal value is dS/ds_i = S * w_i / s_i, so the
    # weakest high-weight term is simultaneously the biggest risk and the
    # biggest available gain.
    contributions, marginal = {}, {}
    for k, v in present.items():
        contributions[k] = (WEIGHTS[k] * math.log(max(v, 1e-12)))
        marginal[k] = WEIGHTS[k] / max(v, 1e-12)
    ln_s = sum(contributions.values()) if not zeroed else float("-inf")
    limiting = min(contributions, key=contributions.get) if contributions else None

    result = {
        "weights": dict(WEIGHTS),
        "ln_S": (None if ln_s == float("-inf") else round(ln_s, 6)),
        "log_contributions": {k: round(v, 4) for k, v in contributions.items()},
        "limiting_term": limiting,
        "limiting_detail": (
            f"{limiting}: s={present.get(limiting):.4f} w={WEIGHTS.get(limiting)} "
            f"contributes {contributions.get(limiting, 0):.4f} to ln S"
            if limiting else None),
        "marginal_value": {k: round(v, 3) for k, v in
                           sorted(marginal.items(), key=lambda kv: -kv[1])},
        "triage": ("fix the limiting term first; dS/ds_i = S*w_i/s_i means the "
                   "weakest high-weight term is the biggest available gain"),
        "subscores": {k: (round(v, 4) if v is not None else None) for k, v in subs.items()},
        "detail": detail,
        "aggregation": "weighted geometric mean: S = prod s_i ** w_i",
        "geometric_total": round(total, 6),
        "arithmetic_total_diagnostic": round(arithmetic, 4),
        "shadow_total_internal_only": round(shadow, 6),
        "zeroed_terms": zeroed,
        "weight_covered": round(covered, 4),
        "missing_terms": missing,
        "geometry_share": GEOMETRY_SHARE,
        "oxyanion_share": OXYANION_SHARE,
        "notes": raw.notes,
    }
    if zeroed:
        result["zero_warning"] = (
            f"terms {zeroed} are 0.000, so the geometric mean is 0.000 regardless "
            "of the other six. Fix these before optimising anything else.")
    # A geometric mean over a subset must have its weights renormalised, or the
    # total is systematically inflated (missing weight behaves like s_i == 1).
    if present and covered > 0:
        if zeroed:
            renorm = 0.0
        else:
            renorm = math.exp(sum((WEIGHTS[k] / covered) * math.log(v)
                                  for k, v in present.items()))
        result["prefilter_score"] = round(renorm, 6)
        result["prefilter_note"] = (
            f"geometric mean over the {len(present)} measurable terms with weights "
            f"renormalised from {covered:.2f} to 1.00. NOT the official score.")

    if missing and strict:
        result["official_score"] = None
        result["verdict"] = (
            f"INCOMPLETE - {missing} unmeasurable, covering only "
            f"{covered:.2f} of 1.00 weight. Refusing to report a total that "
            "would flatter a partially-measured design (fail closed). "
            f"Use prefilter_score={result.get('prefilter_score')} for ranking only."
        )
    else:
        result["official_score"] = round(total, 6)
        result["verdict"] = "COMPLETE"
    return result


def score_replicates(per_replicate: Sequence[RawMetrics], *,
                     ipde_values: Sequence[float] | None = None,
                     ipde_scale: float = 8.0) -> dict:
    """Score EACH replicate, then aggregate the SCORES. Use this, not score().

    WHY THIS EXISTS - a false positive it was built to prevent
    ------------------------------------------------------------
    Averaging geometry and scoring once is only valid when the geometry is
    unimodal. Measured on `ins079_r103_42_F42Y_H103F`:

        per-seed oxyanion dihedral  101.1  96.4  100.0  246.7  244.7
        per-seed s_dih                0.0   0.0    0.0    0.0    0.0
        circular mean 138.8 deg  ->  s_dih 0.3138

    ZERO of five seeds has the dihedral in an allowed band. The design is
    bimodal: it sits in the 50-130 dead zone or the 230-310 dead zone, never
    between them. The circular mean of those two clusters lands at 138.8,
    which IS inside the 130-230 band, and the design was reported at an
    official score of 0.624 - entirely an artefact of the averaging step.

    Scored per replicate, `orient_oxy` is 0.000 on every seed and the design
    scores 0.000, which is the truth.

    This is the handover's own rule (4.3.7, echoed in LESSONS 5.2): rank on
    median plus dispersion or worst case, never on a single flattering
    aggregate. It also removes any need to reason about circular means.

    AGGREGATION, and why `min` for the total
    -----------------------------------------
    The official scalar is a geometric mean, so a term that is zero on ANY
    replicate makes that replicate 0.000. Reporting the median total would
    hide a design that works one time in five. We therefore report:

        per_replicate   every replicate's full score
        median_total    the middle replicate
        worst_total     the worst - the conservative number to quote
        n_nonzero       how many replicates scored at all

    `cons` and `ipde` are whole-batch quantities and are injected into every
    replicate before scoring, since they cannot be computed from one model.
    """
    if not per_replicate:
        raise ValueError("no replicate measurements supplied")

    combined = combine_replicates(per_replicate, ipde_values=ipde_values)
    results = []
    for r in per_replicate:
        r2 = RawMetrics(**{f: getattr(r, f) for f in
                           ("d_base", "d_oxy", "angle2_base", "dihedral12_base",
                            "angle2_oxy", "dihedral12_oxy", "ca_neighbors_mean")})
        # whole-batch terms, identical across replicates by construction
        r2.ipde = combined.ipde
        r2.sd_base = combined.sd_base
        r2.n_replicates = combined.n_replicates
        results.append(score(r2, ipde_scale=ipde_scale, strict=False))

    totals = [s.get("geometric_total", 0.0) or 0.0 for s in results]
    nonzero = [t for t in totals if t > 0]
    ordered = sorted(totals)
    median_total = ordered[len(ordered) // 2]

    # per-term: the worst replicate, because a geometric mean is governed by
    # its weakest component on its weakest replicate
    terms = list(WEIGHTS)
    worst_subs = {}
    for t in terms:
        vals = [s["subscores"].get(t) for s in results
                if s["subscores"].get(t) is not None]
        worst_subs[t] = (min(vals) if vals else None)

    zeroed_any = sorted({t for s in results
                         for t in s.get("zeroed_terms", [])})
    return {
        "method": "score each replicate, then aggregate the scores",
        "n_replicates": len(results),
        "per_replicate_total": [round(t, 6) for t in totals],
        "worst_total": round(min(totals), 6),
        "median_total": round(median_total, 6),
        "best_total": round(max(totals), 6),
        "n_nonzero_replicates": len(nonzero),
        "worst_subscores": {k: (round(v, 4) if v is not None else None)
                            for k, v in worst_subs.items()},
        "terms_zeroed_on_any_replicate": zeroed_any,
        "per_replicate": results,
        "averaged_geometry_total": score(combined, strict=False).get(
            "geometric_total"),
        "warning": (
            "`averaged_geometry_total` averages the geometry and scores once. "
            "It is reported ONLY for comparison and can be badly optimistic on "
            "a bimodal design - it gave 0.624 where every replicate scored "
            "0.000. Quote worst_total, or median_total with the dispersion."),
    }


def rank(scored: Sequence[dict]) -> list[int]:
    """Indices ordered best-first on official_score. Entries with None sort last."""
    idx = list(range(len(scored)))
    idx.sort(key=lambda i: (scored[i].get("official_score") is None,
                            -(scored[i].get("official_score") or 0.0)))
    return idx


# --------------------------------------------------------------------------- backbone gate
#: Thresholds for RFD3's self-reported metrics. These are FREE (already in the
#: output JSON) and directly predictive of the geometry terms: a motif placed
#: 1.2 A off cannot satisfy a 2.48 A contact, so there is no point spending
#: ~80 s of Boltz on it. Tuned from the first smoke design, which had
#: insertion_rmsd 1.22, n_chainbreaks 2, alanine_content 0.45 and scored base=0.0.
#: RECALIBRATED against real 200-timestep output (the originals were tuned on a
#: single 50-timestep design and rejected 2 of 3 backbones):
#:   insertion_rmsd   0.47-0.63 observed  -> 1.00 is comfortable, KEEP (this is
#:                    the one that actually predicts the geometry terms)
#:   max_ca_deviation 0.23-0.47 observed  -> 1.50 comfortable, KEEP
#:   n_chainbreaks    0 observed          -> KEEP at 0, a broken chain is fatal
#:   clashes          0-1 observed        -> 0 was punitive; one inter-residue
#:                    clash in a 200-residue diffusion output is unremarkable
#:                    and gets repacked by LigandMPNN. Relaxed to 2.
#:   alanine_content  0.27-0.44 observed  -> 0.35 was coin-flipping away half the
#:                    sweep. LigandMPNN REPLACES the sequence wholesale, so the
#:                    backbone's alanine fraction says little about the final
#:                    design. Kept only as a sanity check against a degenerate
#:                    poly-Ala output.
BACKBONE_GATE = {
    "insertion_rmsd_max": 1.00,
    "max_ca_deviation_max": 1.50,
    "n_chainbreaks_max": 0,
    "clashes_max": 2,
    "alanine_content_max": 0.70,
    "min_residues": 150,
}


def gate_backbone(metrics: dict, *, thresholds: dict | None = None) -> dict:
    """Reject hopeless backbones before paying for sequence design and folding.

    Because the official scalar is a geometric mean, a design with a misplaced
    motif scores a hard 0.000 no matter how good everything else is. Catching
    that here costs nothing.
    """
    t = {**BACKBONE_GATE, **(thresholds or {})}
    checks: dict[str, dict] = {}

    def chk(name, value, ok, limit):
        checks[name] = {"value": value, "limit": limit, "passed": bool(ok)}

    ins = metrics.get("insertion_rmsd")
    # nan is not a pass: it means the metric could not be computed
    ins_ok = ins is not None and ins == ins and ins <= t["insertion_rmsd_max"]
    chk("insertion_rmsd", ins, ins_ok, t["insertion_rmsd_max"])

    mcd = metrics.get("max_ca_deviation")
    chk("max_ca_deviation", mcd,
        mcd is not None and mcd == mcd and mcd <= t["max_ca_deviation_max"],
        t["max_ca_deviation_max"])

    cb = metrics.get("n_chainbreaks")
    chk("n_chainbreaks", cb, cb is not None and cb <= t["n_chainbreaks_max"],
        t["n_chainbreaks_max"])

    clash = ((metrics.get("n_clashing.interresidue_clashes_w_backbone") or 0)
             + (metrics.get("n_clashing.interresidue_clashes_w_sidechain") or 0))
    chk("clashes", clash, clash <= t["clashes_max"], t["clashes_max"])

    ala = metrics.get("alanine_content")
    chk("alanine_content", ala,
        ala is None or ala <= t["alanine_content_max"], t["alanine_content_max"])

    nres = metrics.get("num_residues")
    chk("num_residues", nres, nres is None or nres >= t["min_residues"],
        t["min_residues"])

    failed = [k for k, v in checks.items() if not v["passed"]]
    return {
        "passed": not failed,
        "failed": failed,
        "checks": checks,
        "rationale": ("the official scalar is a geometric mean, so one bad "
                      "geometry term zeroes the design; reject it before "
                      "spending GPU time on folding"),
    }
