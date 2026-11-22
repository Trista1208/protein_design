"""Dense internal ranking objective, degeneracy alarm, and protocol stamp.

Three instruments that LESSONS.md says were missing, and why each exists.

----------------------------------------------------------------------------
1. A DENSE OBJECTIVE  (LESSONS 1.1, 1.3, 4.3)
----------------------------------------------------------------------------
The official score is a weighted GEOMETRIC mean, so one zero sub-score makes
the total exactly 0.000. Measured yesterday: 5 designs scored, 0 non-zero.
A reward that is 0 for ~100% of samples carries no gradient - it cannot rank,
select, or inform any search.

Worse, the zero is not even monotone in the underlying geometry. The dihedral
sub-score is 0 across the WHOLE dead zone (50-130 and 230-310 deg), so:

    dihedral 237 deg  ->  s_dih = 0.000   (7 deg from the band edge)
    dihedral 290 deg  ->  s_dih = 0.000   (20 deg from the band edge)

Those two designs are not equally bad, and the official score cannot say so.
`rfd4_design_045` sits at 237.2 deg and misses by 3.4 deg; under the official
scalar it is indistinguishable from a design that is nowhere near.

So for INTERNAL ranking we score a continuous BADNESS per official term -
defined everywhere, including deep in the dead zones - then percentile-rank
each term across the batch and combine with the official weights. The result
orders a batch in which every official score is 0.000.

This is a surrogate for SEARCH ONLY. `official.score()` remains the only thing
that may be reported. Never quote a number from this module as a score.

----------------------------------------------------------------------------
2. THE DEGENERACY ALARM  (LESSONS 2.5, 4.4)
----------------------------------------------------------------------------
RFD3 reported d_base 2.486 / d_oxy 3.338 / CA-CA 10.606 IDENTICAL TO THREE
DECIMALS across designs whose motifs landed at 167/114, 128/138 and 58/13 deg.
It was a rigid copy of the input, carrying zero information about the design.

Any metric that does not vary between designs is measuring your input.
`audit()` fails the run on all-missing, >50% missing, <=1 unique value, or
std < 1e-8. LESSONS 4.4: "Would have fired on day 1 here."

----------------------------------------------------------------------------
3. THE PROTOCOL STAMP  (LESSONS 3.1, 3.2)
----------------------------------------------------------------------------
A crystallographically perfect Kemp eliminase scores 0.000 folded WITHOUT an
MSA and 0.775 WITH one. The fold protocol therefore dominates the score, and
nothing in designs.csv or the ledger recorded which protocol produced a number.
That is how the controls (natural, MSA) and the designs (de novo, no MSA)
silently became incomparable.

`protocol_stamp()` is attached to every scored record so two numbers can never
again be compared without the comparison being checkable.
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence

from . import config as C
from .official import LAUKO, TARGETS, WEIGHTS, RAMP_MODE

# --------------------------------------------------------------------------- 1. dense badness

#: Direction of every dense feature, stated EXPLICITLY rather than inferred.
#:
#: LESSONS 2.6: "I proposed three cheap proxies today and verified the sign of
#: none of them. The one I did test inverted." A surrogate whose sign is wrong
#: ranks the batch backwards and looks perfectly healthy while doing it, so the
#: direction is data here, never a convention assumed at the call site.
LOWER_IS_BETTER = {
    "base": True,          # |d_base - 2.52| A
    "oxy": True,           # |d_oxy - 3.10| A
    "ipde": True,          # complex_ipde, an error in A
    "orient_base": True,   # max(angle deviation, dihedral band gap) in deg
    "orient_oxy": True,    # ditto
    "bury": False,         # ca_neighbors_mean, more is more buried
    "cons": True,          # sd_base across replicates, in A
}
assert set(LOWER_IS_BETTER) == set(WEIGHTS), "dense terms must mirror official terms"


def dihedral_band_gap(value: float | None) -> float | None:
    """Degrees of rotation needed to reach the nearest sp2 band; 0.0 if in band.

    This is the term that makes the dead zones rankable. The official
    sub-score is a flat 0 across 50-130 and 230-310 deg; this is the distance
    to the nearest edge of an allowed band, so 237 deg (gap 7) is correctly
    ordered ahead of 290 deg (gap 20) instead of tying at zero.

    Returns ``None`` for an unmeasurable dihedral - an ABSENT quantity is not
    a good one (LESSONS 1.5), and the caller must not let it default to 0.0.
    """
    if value is None:
        return None
    v = float(value) % 360.0
    lo180, hi180 = LAUKO["dihedral_sp2_180"]        # 130, 230
    lo0, hi0 = LAUKO["dihedral_sp2_0"]              # 310, 50 - wraps through 0
    in_180 = lo180 <= v <= hi180
    in_0 = (v >= lo0) or (v <= hi0)
    if in_180 or in_0:
        return 0.0
    if hi0 < v < lo180:                             # dead zone 50-130
        return min(v - hi0, lo180 - v)
    return min(v - hi180, lo0 - v)                  # dead zone 230-310


def angle_deviation(value: float | None) -> float | None:
    """Absolute deviation from the 120 deg angle target, in degrees.

    Unlike the official angle sub-score this stays informative outside the
    80-160 deg window instead of saturating at 0.
    """
    if value is None:
        return None
    return abs(float(value) - TARGETS["angle"]["target"])


def _orient_badness(angle: float | None, dihedral: float | None) -> float | None:
    """Worse of the two orientation defects, in degrees.

    The official term is ``min(angle_score, dihedral_score)`` - BOTH must be
    good (LESSONS 1.3). The dense mirror of a min over goodness is a MAX over
    badness. Both inputs are already in degrees, so they are comparable without
    a scale factor that would need its own justification.
    """
    da, dg = angle_deviation(angle), dihedral_band_gap(dihedral)
    if da is None or dg is None:
        return None
    return max(da, dg)


def badness(raw) -> dict[str, float | None]:
    """Continuous per-term defect for one design. Lower is better except ``bury``.

    ``raw`` is an :class:`official.RawMetrics`. Every value is defined wherever
    the underlying measurement exists, including where the official sub-score
    is a flat 0, which is the entire point of this module.
    """
    out: dict[str, float | None] = {}
    out["base"] = (None if raw.d_base is None
                   else abs(raw.d_base - TARGETS["d_base"]["target"]))
    out["oxy"] = (None if raw.d_oxy is None
                  else abs(raw.d_oxy - TARGETS["d_oxy"]["target"]))
    out["ipde"] = raw.ipde
    out["orient_base"] = _orient_badness(raw.angle2_base, raw.dihedral12_base)
    out["orient_oxy"] = _orient_badness(raw.angle2_oxy, raw.dihedral12_oxy)
    out["bury"] = raw.ca_neighbors_mean
    out["cons"] = raw.sd_base
    return out


def diagnostics(raw) -> dict[str, float | None]:
    """The individual dense quantities, for the per-objective band tables.

    Reported alongside every design because LESSONS 7.3 is explicit: "Report
    every term, every design, every replicate - especially the dihedrals. Six
    of eight designs died on a term no table showed."
    """
    return {
        "dev_d_base_A": (None if raw.d_base is None
                         else round(abs(raw.d_base - TARGETS["d_base"]["target"]), 4)),
        "dev_d_oxy_A": (None if raw.d_oxy is None
                        else round(abs(raw.d_oxy - TARGETS["d_oxy"]["target"]), 4)),
        "dev_angle_base_deg": _r(angle_deviation(raw.angle2_base)),
        "dev_angle_oxy_deg": _r(angle_deviation(raw.angle2_oxy)),
        "gap_dihedral_base_deg": _r(dihedral_band_gap(raw.dihedral12_base)),
        "gap_dihedral_oxy_deg": _r(dihedral_band_gap(raw.dihedral12_oxy)),
        "dihedral_base_deg": _r(raw.dihedral12_base, 2),
        "dihedral_oxy_deg": _r(raw.dihedral12_oxy, 2),
    }


def _r(v, nd: int = 3):
    return None if v is None else round(float(v), nd)


# --------------------------------------------------------------------------- 2. percentile rank

def percentile(values: Sequence[float | None], *, lower_is_better: bool) -> list[float | None]:
    """Percentile rank in [0, 1], 1.0 = best. Ties share the average rank.

    ``None`` in, ``None`` out: a missing measurement must not be imputed to the
    middle of the batch, which would flatter an unmeasurable design (the same
    failure mode as ``official.score(strict=False)`` in LESSONS 1.5).
    """
    present = [(i, float(v)) for i, v in enumerate(values) if v is not None]
    out: list[float | None] = [None] * len(values)
    n = len(present)
    if n == 0:
        return out
    if n == 1:
        out[present[0][0]] = 0.5     # a single design has no rank within a batch
        return out
    present.sort(key=lambda iv: iv[1], reverse=not lower_is_better)
    # average ranks over ties so that a degenerate column maps to a constant
    # 0.5 rather than to an arbitrary order - the alarm in audit() then fires
    # on it instead of it quietly driving the ranking.
    i = 0
    while i < n:
        j = i
        while j + 1 < n and present[j + 1][1] == present[i][1]:
            j += 1
        avg_rank = (i + j) / 2.0
        for k in range(i, j + 1):
            out[present[k][0]] = 1.0 - avg_rank / (n - 1)
        i = j + 1
    return out


def rank_batch(badnesses: Sequence[dict[str, float | None]],
               *, weights: dict[str, float] | None = None) -> list[dict]:
    """Percentile-rank every term across the batch and combine with the weights.

    LESSONS 4.3: percentile-rank scalarisation, NEVER a product, for internal
    ranking - "Gives a usable ordering even when every official score is
    0.000." A missing term contributes 0 to its own axis only, so it costs that
    design exactly that term's weight instead of nullifying the whole row.

    Returns one dict per input design, in input order:
      ``dense_rank``      weighted sum of percentiles, in [0, 1], higher better
      ``term_percentile`` the per-term percentiles behind it
      ``weight_present``  fraction of official weight actually measured
      ``worst_term``      biggest weighted shortfall among the terms that VARY
      ``constant_terms``  terms identical across the batch - see below

    ``worst_term`` IS NOT ``official.limiting_term``, AND IT IS NOT THE ZERO.
    It answers "which term is this design worst at RELATIVE TO THIS BATCH",
    which is the right question for choosing the next move in a search. It is
    the wrong question for diagnosing a zero, because a percentile is blind to
    whether the absolute value is acceptable.

    Worked example, the real `rfd4_design_045` in a 4-design batch:
      ``worst_term`` = ipde           - mid-pack in a batch, weight 0.20
      ``limiting_term`` = orient_oxy  - the hard 0.000 that zeroes the design
    Its ipde sub-score is 0.914 and needs no work at all. Optimise against
    ``worst_term`` and you would tune a healthy term while the zero stands.

    So: read ``official.limiting_term`` to find the zero, and ``worst_term``
    only to break ties among designs that have no zero left. Both columns are
    published side by side for exactly this reason.
    """
    w = dict(weights or WEIGHTS)
    terms = list(LOWER_IS_BETTER)
    cols = {t: percentile([b.get(t) for b in badnesses],
                          lower_is_better=LOWER_IS_BETTER[t]) for t in terms}

    # A term identical across the whole batch percentile-ranks to a flat 0.5,
    # which would then generate a phantom shortfall of w*0.5 for EVERY design -
    # and because `base` carries the largest weight, the highest-weight constant
    # term would always be named as the thing to fix. It cannot be: a term that
    # does not vary is not what separates these designs (LESSONS 2.5). So
    # `worst_term` is chosen only among terms that actually vary, and the
    # constant ones are listed explicitly for audit() to rule on.
    varying, constant = [], []
    for t in terms:
        vals = {b.get(t) for b in badnesses if b.get(t) is not None}
        (varying if len(vals) > 1 else constant).append(t)

    out = []
    for i in range(len(badnesses)):
        pcts = {t: cols[t][i] for t in terms}
        total = sum(w[t] * (pcts[t] or 0.0) for t in terms)
        present = sum(w[t] for t in terms if pcts[t] is not None)
        shortfall = {t: w[t] * (1.0 - (pcts[t] or 0.0)) for t in terms}
        pool = {t: shortfall[t] for t in varying} or shortfall
        out.append({
            "dense_rank": round(total, 6),
            "term_percentile": {t: _r(pcts[t], 4) for t in terms},
            "weight_present": round(present, 4),
            "worst_term": max(pool, key=pool.get) if pool else None,
            "weighted_shortfall": {t: round(v, 4) for t, v in
                                   sorted(shortfall.items(), key=lambda kv: -kv[1])},
            "constant_terms": list(constant),
            "note": "INTERNAL SEARCH SURROGATE - not a score, never report it",
        })
    return out


# --------------------------------------------------------------------------- 3. degeneracy alarm

#: A column is degenerate when it cannot distinguish designs. Thresholds from
#: LESSONS 4.4 / ESM/audit_esm_pipeline_run.py:428-452.
DEGENERACY = {
    "max_missing_frac": 0.50,
    "min_unique": 2,
    "min_std": 1e-8,
    "min_rows": 3,          # below this a batch has no distribution to judge
}


def audit(rows: Sequence[dict], columns: Iterable[str] | None = None,
          *, thresholds: dict | None = None) -> dict:
    """Fail the run when a score column carries no information.

    This exists because RFD3 emitted d_base/d_oxy/CA-CA identical to three
    decimals across structurally different designs and the campaign ranked on
    them anyway (LESSONS 2.5). Identical numbers across different designs are
    inherited, not earned.

    ``rows`` are flat dicts (one per design). Returns a verdict; call
    :func:`assert_ok` to turn a failure into an exception at a launch gate.
    """
    t = {**DEGENERACY, **(thresholds or {})}
    n = len(rows)
    if columns is None:
        columns = sorted({k for r in rows for k, v in r.items()
                          if isinstance(v, (int, float)) and not isinstance(v, bool)})
    cols = list(columns)

    if n < t["min_rows"]:
        return {"passed": True, "skipped": True, "n_rows": n,
                "failures": [], "columns": {},
                "reason": f"only {n} row(s); need >= {t['min_rows']} to judge a "
                          "distribution. NOT a pass, just not yet judgeable."}

    report: dict[str, dict] = {}
    failures: list[str] = []
    for col in cols:
        vals = [r.get(col) for r in rows]
        present = [float(v) for v in vals
                   if isinstance(v, (int, float)) and not isinstance(v, bool)
                   and not math.isnan(float(v))]
        missing_frac = 1.0 - len(present) / n
        uniq = len(set(present))
        if len(present) >= 2:
            mu = sum(present) / len(present)
            std = math.sqrt(sum((x - mu) ** 2 for x in present) / (len(present) - 1))
        else:
            std = 0.0
        why = []
        if not present:
            why.append("all values missing or NaN")
        else:
            if missing_frac > t["max_missing_frac"]:
                why.append(f"{missing_frac:.0%} missing > {t['max_missing_frac']:.0%}")
            if uniq < t["min_unique"]:
                why.append(f"{uniq} unique value(s) - constant column, so it is "
                           "measuring your input, not your designs")
            elif std < t["min_std"]:
                why.append(f"std {std:.2e} < {t['min_std']:.0e}")
        report[col] = {"n_present": len(present), "missing_frac": round(missing_frac, 4),
                       "n_unique": uniq, "std": std, "degenerate": bool(why),
                       "reasons": why}
        if why:
            failures.append(col)
    return {
        "passed": not failures,
        "skipped": False,
        "n_rows": n,
        "failures": failures,
        "columns": report,
        "rationale": ("a column that does not vary between designs is measuring "
                      "the input, not the designs (LESSONS 2.5); ranking on it "
                      "silently orders by noise"),
    }


class DegenerateScoresError(RuntimeError):
    """Raised when a score column carries no information about the designs."""


def assert_ok(rows: Sequence[dict], columns: Iterable[str] | None = None,
              **kw) -> dict:
    """:func:`audit`, but raise on failure."""
    rep = audit(rows, columns, **kw)
    if not rep["passed"]:
        detail = "; ".join(
            f"{c}: {', '.join(rep['columns'][c]['reasons'])}" for c in rep["failures"])
        raise DegenerateScoresError(
            f"{len(rep['failures'])} degenerate score column(s) over "
            f"{rep['n_rows']} designs -> {detail}")
    return rep


#: Designs to score before the campaign checkpoint can judge anything.
#:
#: NOT a launch gate: at launch there are zero rows, so a pre-flight check is
#: theatre. The real failure (LESSONS 2.5) is only visible once several designs
#: have been scored and their geometry turns out identical. At the measured
#: ~34 s/design (500 designs = 4h41m, LESSONS 4.8) this checkpoint costs about
#: six minutes out of a full campaign - cheap insurance against waking up to
#: 500 copies of the same number.
CAMPAIGN_CHECKPOINT_AFTER = 10


def campaign_checkpoint(badnesses: Sequence[dict], *,
                        after: int = CAMPAIGN_CHECKPOINT_AFTER) -> dict:
    """Judge an in-flight campaign, separating fatal degeneracy from advisory.

    The distinction matters, and conflating the two is how a useful alarm turns
    into one people switch off:

    * **FATAL** - a term is constant, or has effectively zero variance, across
      scored designs. That is the LESSONS 2.5 signature: RFD3 reported d_base
      2.486 / d_oxy 3.338 identical to three decimals across designs whose
      motifs landed at 167/114, 128/138 and 58/13 deg. The pipeline is copying
      its input and the rest of the campaign cannot learn anything. Stop.
    * **ADVISORY** - a term is often missing. Designs legitimately fail gates
      and fail to fold, so a high missing fraction early in a sweep is
      expected, not broken. Report it; do not kill the run for it.

    Returns a verdict; the caller decides whether to raise.
    """
    n = len(badnesses)
    if n < after:
        return {"judged": False, "n_scored": n, "needs": after,
                "fatal": [], "advisory": [], "passed": True,
                "reason": f"{n}/{after} scored designs; not yet judgeable. "
                          "NOT a pass."}

    rep = audit([dict(b) for b in badnesses], columns=list(LOWER_IS_BETTER),
                thresholds={"min_rows": after})
    fatal, advisory = [], []
    for col in rep["failures"]:
        reasons = rep["columns"][col]["reasons"]
        if any(("unique value" in r) or ("std" in r) or ("all values missing" in r)
               for r in reasons):
            fatal.append(col)
        else:
            advisory.append(col)
    return {
        "judged": True, "n_scored": n, "needs": after,
        "fatal": fatal, "advisory": advisory,
        "passed": not fatal,
        "columns": rep["columns"],
        "reason": ("constant or zero-variance term(s) across scored designs: "
                   "the pipeline is measuring its input, not its designs"
                   if fatal else "no term is degenerate"),
    }


# --------------------------------------------------------------------------- 4. protocol stamp

def protocol_stamp(*, replicates: int | None = None,
                   extra: dict | None = None) -> dict:
    """Everything that changes a score without changing the design.

    Attached to every scored record. LESSONS 3.1: the SAME sequence and seed
    scored 0.000 without an MSA and 0.775 with one, so a score without its
    protocol is not a measurement. LESSONS 3.2: "Record the MSA setting
    alongside every score. This is how the control and the designs silently
    became incomparable."
    """
    b = dict(C.BOLTZ_SETTINGS)
    msa = ("none:single_sequence" if b.get("single_sequence_mode")
           else ("msa_server" if b.get("use_msa_server") else "msa:unspecified"))
    stamp = {
        "msa_mode": msa,
        "single_sequence_mode": bool(b.get("single_sequence_mode")),
        "use_msa_server": bool(b.get("use_msa_server")),
        "recycling_steps": b.get("recycling_steps"),
        "diffusion_samples": b.get("diffusion_samples"),
        "seeds": list(C.BOLTZ_SEEDS),
        "replicates": replicates,
        "ligand_ccd": C.LIGAND_CCD,
        "ramp_mode": RAMP_MODE,
        "targets": {k: v["target"] for k, v in TARGETS.items()},
        "comparability_warning": (
            "scores are comparable ONLY across records with an identical "
            "msa_mode. De novo designs have no MSA; natural controls do. "
            "Folding HG3.17 single-sequence scores 0.000 (LESSONS 3.1)."),
    }
    if extra:
        stamp.update(extra)
    return stamp


def comparable(a: dict, b: dict) -> bool:
    """True when two protocol stamps may have their scores compared directly."""
    keys = ("msa_mode", "ramp_mode", "ligand_ccd", "recycling_steps")
    return all(a.get(k) == b.get(k) for k in keys)
