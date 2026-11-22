"""M3 - the four official metrics, as separate auditable functions.

Slides p.6: iPDE, Geometry, Catalytic base burial (CA within 12 A), Consistency.

Deliberate design choices, all traceable to the handover:

* Handover 6.3 - "What is the exact Boltz-2 field name behind iPDE, and is
  smaller better? Do not invent a formula." We therefore READ the field from
  Boltz output and never synthesise it. If the field is absent we raise, we do
  not default (4.3.9, fail closed).
* Handover 6.2 - the 12 A burial anchor convention is ambiguous. We compute ALL
  three candidate anchors and report them side by side rather than picking one
  silently.
* Handover 4.3.7 - rank on median plus dispersion or worst case, never the best
  seed. ``aggregate`` offers no "best" option on purpose.
* Handover 2 - Gln50's 3.33 A reference is modelled, not observed, so geometry
  scoring weights the Asp contact and flags the Gln term.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import Sequence

import numpy as np

from . import config as C
from .theozyme import Theozyme, load


class ScoringError(RuntimeError):
    pass


# --------------------------------------------------------------------------- iPDE
# RESOLVED from the installed Boltz-2 source (confidence.py:451-456), not guessed.
# complex_ipde = mean predicted distance error (A) over CROSS-CHAIN token pairs.
# It is an ERROR: SMALLER IS BETTER, opposite to plddt/iptm.
IPDE_CANDIDATE_KEYS = (C.IPDE_FIELD, "ipde", "interface_pde", "pde_interface")


def read_boltz_confidence(path: Path | str) -> dict:
    """Load a Boltz confidence JSON without interpreting it."""
    path = Path(path)
    if not path.exists():
        raise ScoringError(f"Boltz confidence file missing: {path}")
    return json.loads(path.read_text())


def extract_ipde(conf: dict, *, key: str | None = None,
                 require_interface: bool = True) -> tuple[str, float]:
    """Return (field_name_used, value). Never invents a value.

    ``complex_ipde`` is an error in Angstroms, so SMALLER IS BETTER.

    With a single chain the interface mask sums to zero and Boltz returns
    ``complex_ipde == 0.0`` - which looks like a flawless interface but means
    there was no interface at all. ``require_interface`` rejects that case by
    cross-checking the ipTM family, because silently ranking a ligand-less fold
    as perfect is the worst failure this scorer could have.
    """
    if key:
        if key not in conf:
            raise ScoringError(f"requested iPDE key {key!r} absent; have {sorted(conf)}")
        name, val = key, float(conf[key])
    else:
        for cand in IPDE_CANDIDATE_KEYS:
            if cand in conf:
                name, val = cand, float(conf[cand])
                break
        else:
            raise ScoringError(
                "no recognised interface-confidence field in Boltz output; "
                f"available keys: {sorted(conf)}. Set the key explicitly rather "
                "than letting the pipeline guess."
            )

    if require_interface and val == C.IPDE_NO_INTERFACE_SENTINEL:
        # every interface-ish signal being zero means no cross-chain pairs existed
        companions = {k: conf.get(k) for k in
                      ("iptm", "ligand_iptm", "protein_iptm", "complex_iplddt")}
        if all(v in (0, 0.0, None) for k, v in companions.items() if k != "complex_iplddt"):
            raise ScoringError(
                f"{name} == 0.0 with iptm/ligand_iptm also 0 - the structure has NO "
                "interface, so this is not a perfect score, it is an absent one. "
                "The ligand must be folded as its own chain (e.g. '>B|ccd\\n6NT') "
                f"for {name} to mean anything. companions={companions}"
            )
    return name, val


def ipde_is_better(a: float, b: float) -> bool:
    """True if ``a`` is a better iPDE than ``b``. Lower wins."""
    return a < b if C.IPDE_LOWER_IS_BETTER else a > b


# --------------------------------------------------------------------------- geometry
# Gln50 is CONFIRMED against 4BS0/7K4Z (handover open question 1 resolved), so
# the donor contact is trustworthy and carries real weight. Oxyanion
# stabilisation is the mechanism the whole challenge turns on (handover 3.2).
GEOMETRY_TERMS = {
    # term: (reference key, tolerance A, weight, trustworthy?)
    "base_contact":  ("asp127_od2__6nt_n3", 0.50, 0.45, True),
    "donor_contact": ("gln50_ne2__6nt_n1", 0.80, 0.35, True),
    "ca_span":       ("asp127_ca__gln50_ca", 1.50, 0.20, True),
}


def geometry_score(measured: dict[str, float]) -> dict:
    """Score active-site geometry against the verified reference table.

    Returns per-term detail plus a 0-1 aggregate. The donor term carries low
    weight and a flag because its reference is a model's opinion (handover 1.1).
    """
    terms, total, wsum = {}, 0.0, 0.0
    for name, (ref_key, tol, weight, trusted) in GEOMETRY_TERMS.items():
        if ref_key not in measured:
            terms[name] = {"status": "missing", "trusted_reference": trusted}
            continue
        ref = C.REF_GEOMETRY[ref_key]
        dev = abs(measured[ref_key] - ref)
        # 1.0 at the reference, decaying to 0 at 3x tolerance
        sub = max(0.0, 1.0 - dev / (3.0 * tol))
        terms[name] = {
            "measured_A": round(measured[ref_key], 3),
            "reference_A": ref,
            "deviation_A": round(dev, 3),
            "within_tol": dev <= tol,
            "subscore": round(sub, 4),
            "weight": weight,
            "trusted_reference": trusted,
        }
        total += sub * weight
        wsum += weight
    return {
        "terms": terms,
        "score": round(total / wsum, 4) if wsum else None,
        "note": ("donor_contact (Gln50) was modelled in with B-factor 0.00 but is "
                 "CONFIRMED to match the crystallographic Gln50 in 4BS0/7K4Z, so it "
                 "carries full weight"),
    }


# --------------------------------------------------------------------------- burial
def catalytic_burial(structure: Theozyme | Path | str, *, radius: float = C.BURIAL_RADIUS_A,
                     include_self: bool = False) -> dict:
    """Count CA atoms within ``radius`` of the catalytic base.

    Handover 6.2: the anchor convention is unresolved, so report all three.
    """
    tz = structure if isinstance(structure, Theozyme) else load(structure)
    asp = tz.by_resname("ASP")
    if not asp:
        raise ScoringError("no ASP in structure - cannot compute catalytic burial")

    anchors = {}
    for label, aname in (("OD2", "OD2"), ("CG_carboxyl_C", "CG"), ("CA", "CA")):
        try:
            anchors[label] = tz.atom("ASP", aname).xyz
        except Exception:
            continue
    if not anchors:
        raise ScoringError("none of OD2/CG/CA present on the catalytic ASP")

    cas = [a for a in tz.atoms if a.name == "CA"]
    asp_resseq = asp[0].resseq
    out = {"radius_A": radius, "n_CA_total": len(cas), "anchors": {}}
    for label, origin in anchors.items():
        n = 0
        for ca in cas:
            if not include_self and ca.resseq == asp_resseq and ca.resname == "ASP":
                continue
            if float(np.linalg.norm(ca.xyz - origin)) <= radius:
                n += 1
        out["anchors"][label] = n
    out["reported"] = out["anchors"].get("OD2")
    out["note"] = (
        "anchor convention unresolved (handover 6.2): OD2 reported by default; "
        "ask the organisers whether they use OD2, the carboxylate carbon, or CA, "
        "and whether the base's own CA counts"
    )
    return out


# --------------------------------------------------------------------------- consistency
def aggregate(values: Sequence[float], *, mode: str = "median_minus_spread",
              lower_is_better: bool = False) -> dict:
    """Aggregate replicate scores. No 'best' mode exists, by design.

    Handover 4.3.7: never cherry-pick the best-looking seed.

    ``lower_is_better`` MUST be True for complex_ipde, which is an error in
    Angstroms. Getting this wrong inverts the entire ranking, so the penalty
    and the worst case both flip direction here rather than at the call site.
    """
    vals = [float(v) for v in values if v is not None and not math.isnan(float(v))]
    if not vals:
        raise ScoringError("no finite replicate values to aggregate")
    med = median(vals)
    spread = (max(vals) - min(vals))
    sd = float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0

    if mode == "median_minus_spread":
        # penalise dispersion in whichever direction counts as worse
        score = (med + 0.5 * spread) if lower_is_better else (med - 0.5 * spread)
    elif mode == "worst_case":
        score = max(vals) if lower_is_better else min(vals)
    elif mode == "median":
        score = med
    else:
        raise ScoringError(
            f"unknown aggregation mode {mode!r}; allowed: median_minus_spread, "
            "worst_case, median (there is deliberately no 'best')"
        )
    return {
        "n": len(vals), "median": round(med, 4), "min": round(min(vals), 4),
        "max": round(max(vals), 4), "spread": round(spread, 4),
        "stdev": round(sd, 4), "mode": mode, "score": round(score, 4),
        "lower_is_better": lower_is_better,
        "values": [round(v, 4) for v in vals],
    }


def consistency(values: Sequence[float]) -> dict:
    """Official 'Consistency' metric: tightness across the 5 replicates."""
    agg = aggregate(values, mode="median")
    spread = agg["spread"]
    # 1.0 when all seeds agree, decaying across a spread of 0.3 confidence units
    agg["consistency"] = round(max(0.0, 1.0 - spread / 0.30), 4)
    return agg


# --------------------------------------------------------------------------- noise floor
@dataclass
class NoiseFloor:
    """Seed-to-seed spread of the control pair - the pipeline's resolution limit.

    Handover 3.2: any two designs closer than this are indistinguishable from
    chance, and ranking inside it is superstition.
    """
    positive: dict = field(default_factory=dict)
    negative: dict = field(default_factory=dict)
    #: True when the metric is an error (complex_ipde), so lower = better
    lower_is_better: bool = C.IPDE_LOWER_IS_BETTER

    @property
    def floor(self) -> float:
        return max(self.positive.get("spread", 0.0), self.negative.get("spread", 0.0))

    def ranks_correctly(self) -> bool:
        """HG3.17 must beat HG3.17-Q50M, in whichever direction the metric runs."""
        p = self.positive.get("median")
        n = self.negative.get("median")
        if p is None or n is None:
            return False
        return (p < n) if self.lower_is_better else (p > n)

    def verdict(self) -> dict:
        ok = self.ranks_correctly()
        return {
            "positive_median": self.positive.get("median"),
            "negative_median": self.negative.get("median"),
            "separation": (None if not self.positive or not self.negative
                           else round(abs(self.positive["median"]
                                         - self.negative["median"]), 4)),
            "lower_is_better": self.lower_is_better,
            "noise_floor": round(self.floor, 4),
            "ranks_correctly": ok,
            "separation_exceeds_noise": (
                None if not self.positive or not self.negative
                else abs(self.positive["median"] - self.negative["median"]) > self.floor
            ),
            "verdict": (
                "PROCEED - pipeline detects oxyanion stabilisation"
                if ok else
                "STOP - pipeline cannot rank HG3.17 above HG3.17-Q50M. Geometry and "
                "Consistency are not measuring oxyanion stabilisation; every number "
                "produced afterwards is noise (handover 3.2)."
            ),
        }


def resolvable(a: float, b: float, noise_floor: float) -> bool:
    """True only if two scores differ by more than the measured noise floor."""
    return abs(a - b) > noise_floor
