"""Tests for the dense objective, the degeneracy alarm and the protocol stamp.

The band-gap cases use the REAL measured dihedrals from LESSONS.md so that a
regression shows up as a familiar number going wrong:

  rfd4_design_045, Boltz refold, 5 seeds : 233.4 - 239.8 deg, mean 237.2
  rfd4_design_045, OpenFold3             : 86.1 deg
  theozyme / crystal                     : 163 - 198 deg (in band)
"""

from __future__ import annotations

import math

import pytest

from kemp import objective as O
from kemp.official import LAUKO, RawMetrics, WEIGHTS


# ----------------------------------------------------------------- band gap
def test_in_band_dihedrals_have_zero_gap():
    for v in (130.0, 180.0, 198.0, 230.0, 310.0, 350.0, 0.0, 50.0):
        assert O.dihedral_band_gap(v) == 0.0, v


def test_d045_misses_the_band_by_the_measured_amount():
    # LESSONS 6: mean 237.2 deg against a band edge of 230 -> misses by 3.4 deg
    # at the closest seed (233.4) and 7.2 at the mean.
    assert O.dihedral_band_gap(233.4) == pytest.approx(3.4, abs=1e-9)
    assert O.dihedral_band_gap(237.2) == pytest.approx(7.2, abs=1e-9)
    assert O.dihedral_band_gap(239.8) == pytest.approx(9.8, abs=1e-9)


def test_dead_zone_designs_are_ordered_where_the_official_score_ties():
    # Both score a flat 0.000 officially. The whole point of this module is
    # that 237 is nearer the band than 290 and must rank ahead of it.
    near, far = O.dihedral_band_gap(237.0), O.dihedral_band_gap(290.0)
    assert near == pytest.approx(7.0)
    assert far == pytest.approx(20.0)
    assert near < far


def test_openfold3_and_boltz_land_in_different_dead_zones():
    # LESSONS 3.8: OF3 chose 86.1 deg, Boltz ~237 deg. Both out of band, in
    # DIFFERENT dead zones - so both must report a non-zero gap.
    assert O.dihedral_band_gap(86.1) == pytest.approx(36.1, abs=1e-9)
    assert O.dihedral_band_gap(237.0) > 0


def test_band_gap_wraps_and_is_periodic():
    for v in (237.0, 86.1, 180.0, 300.0):
        g = O.dihedral_band_gap(v)
        assert O.dihedral_band_gap(v + 360.0) == pytest.approx(g, abs=1e-9)
        assert O.dihedral_band_gap(v - 360.0) == pytest.approx(g, abs=1e-9)


def test_band_gap_peaks_at_the_middle_of_each_dead_zone():
    # dead zones are (50,130) and (230,310); midpoints 90 and 270, half-width 40
    assert O.dihedral_band_gap(90.0) == pytest.approx(40.0)
    assert O.dihedral_band_gap(270.0) == pytest.approx(40.0)
    for v in (60.0, 120.0, 240.0, 300.0):
        assert 0 < O.dihedral_band_gap(v) < 40.0


def test_absent_dihedral_is_none_not_zero():
    # LESSONS 1.5: an absent quantity must not look like a perfect one.
    assert O.dihedral_band_gap(None) is None
    assert O.angle_deviation(None) is None
    assert O._orient_badness(None, 180.0) is None
    assert O._orient_badness(120.0, None) is None


# ----------------------------------------------------------------- badness
def test_orient_badness_takes_the_worse_of_angle_and_dihedral():
    # official term is min(angle_score, dihedral_score) - both must be good,
    # so the dense mirror is a MAX over defects.
    assert O._orient_badness(120.0, 180.0) == 0.0          # both perfect
    assert O._orient_badness(120.0, 237.0) == pytest.approx(7.0)   # dihedral limits
    assert O._orient_badness(100.0, 180.0) == pytest.approx(20.0)  # angle limits
    assert O._orient_badness(100.0, 237.0) == pytest.approx(20.0)  # angle is worse


def test_badness_mirrors_the_official_terms_exactly():
    assert set(O.badness(RawMetrics())) == set(WEIGHTS)


def test_badness_is_zero_at_the_crystal_targets():
    raw = RawMetrics(d_base=2.52, d_oxy=3.10, angle2_base=120.0,
                     dihedral12_base=180.0, angle2_oxy=120.0,
                     dihedral12_oxy=180.0)
    b = O.badness(raw)
    for term in ("base", "oxy", "orient_base", "orient_oxy"):
        assert b[term] == pytest.approx(0.0), term


def test_badness_is_defined_where_the_official_subscore_is_flat_zero():
    # d_oxy 10.216 A is far outside the 2.0-3.6 window, so oxy scores 0.000.
    # The dense term must still distinguish it from 12.698 A (LESSONS 3.1).
    worse = O.badness(RawMetrics(d_oxy=12.698))["oxy"]
    bad = O.badness(RawMetrics(d_oxy=10.216))["oxy"]
    assert bad < worse


# ----------------------------------------------------------------- percentile
def test_percentile_best_is_one_worst_is_zero():
    p = O.percentile([1.0, 2.0, 3.0], lower_is_better=True)
    assert p == [1.0, 0.5, 0.0]
    p = O.percentile([1.0, 2.0, 3.0], lower_is_better=False)
    assert p == [0.0, 0.5, 1.0]


def test_percentile_preserves_none():
    p = O.percentile([1.0, None, 3.0], lower_is_better=True)
    assert p[1] is None
    assert p[0] == 1.0 and p[2] == 0.0


def test_percentile_of_a_constant_column_is_a_flat_half():
    # a degenerate column must not impose an arbitrary order; audit() is what
    # catches it, and it can only do so if the ranking has not already used it.
    assert O.percentile([5.0] * 4, lower_is_better=True) == [0.5] * 4


def test_single_design_gets_no_rank_within_a_batch():
    assert O.percentile([2.0], lower_is_better=True) == [0.5]


def test_empty_and_all_missing_columns():
    assert O.percentile([], lower_is_better=True) == []
    assert O.percentile([None, None], lower_is_better=True) == [None, None]


def test_bury_is_the_one_higher_is_better_term():
    assert O.LOWER_IS_BETTER["bury"] is False
    assert all(v for k, v in O.LOWER_IS_BETTER.items() if k != "bury")


# ----------------------------------------------------------------- rank_batch
def _raw(**kw) -> RawMetrics:
    base = dict(d_base=2.52, d_oxy=3.10, angle2_base=120.0, dihedral12_base=180.0,
                angle2_oxy=120.0, dihedral12_oxy=180.0, ca_neighbors_mean=20.0,
                ipde=0.5, sd_base=0.07)
    base.update(kw)
    return RawMetrics(**base)


def test_rank_batch_orders_a_batch_whose_official_scores_are_all_zero():
    from kemp import official
    # Three designs with out-of-band oxyanion dihedrals. All score a hard
    # 0.000 officially, and the dense rank must still order them by how far
    # they are from the band: 237 (gap 7) > 290 (gap 20) > 270 (gap 40, the
    # dead-zone centre and therefore the worst of the three).
    dihedrals = [237.0, 290.0, 270.0]
    raws = [_raw(dihedral12_oxy=d) for d in dihedrals]
    for r in raws:
        assert official.score(r, strict=False)["geometric_total"] == 0.0
    ranked = O.rank_batch([O.badness(r) for r in raws])
    got = [r["dense_rank"] for r in ranked]
    assert got[0] > got[1] > got[2]
    # the ordering is exactly the band-gap ordering
    by_rank = sorted(range(3), key=lambda i: -got[i])
    by_gap = sorted(range(3), key=lambda i: O.dihedral_band_gap(dihedrals[i]))
    assert by_rank == by_gap


def test_rank_batch_names_the_term_to_fix():
    raws = [_raw(dihedral12_oxy=290.0), _raw(), _raw()]
    ranked = O.rank_batch([O.badness(r) for r in raws])
    assert ranked[0]["worst_term"] == "orient_oxy"


def test_worst_term_is_not_the_zeroing_term():
    """The documented trap: worst_term is batch-relative, not a zero diagnosis.

    Real numbers for rfd4_design_045 in a 4-design batch. Its ipde sub-score is
    0.914 and needs no work, but ipde carries weight 0.20 against orient_oxy's
    0.08, so the percentile shortfall names ipde. The zero must be read from
    official.limiting_term instead.
    """
    from kemp import official
    d045 = official.RawMetrics(
        d_base=2.858, d_oxy=3.103, angle2_base=120.0, dihedral12_base=218.2,
        angle2_oxy=120.0, dihedral12_oxy=237.2, ca_neighbors_mean=20.0,
        ipde=0.688, sd_base=0.0759, n_replicates=5)
    others = [
        official.RawMetrics(d_base=2.532, d_oxy=3.102, angle2_base=120.0,
                            dihedral12_base=180.0, angle2_oxy=120.0,
                            dihedral12_oxy=180.0, ca_neighbors_mean=20.0,
                            ipde=0.53, sd_base=0.068, n_replicates=5),
        official.RawMetrics(d_base=3.037, d_oxy=10.216, angle2_base=120.0,
                            dihedral12_base=180.0, angle2_oxy=120.0,
                            dihedral12_oxy=180.0, ca_neighbors_mean=20.33,
                            ipde=0.655, sd_base=0.30, n_replicates=5),
        official.RawMetrics(d_base=2.157, d_oxy=12.698, angle2_base=120.0,
                            dihedral12_base=180.0, angle2_oxy=120.0,
                            dihedral12_oxy=180.0, ca_neighbors_mean=20.0,
                            ipde=3.00, sd_base=0.30, n_replicates=5),
    ]
    batch = [d045, *others]
    ranked = O.rank_batch([O.badness(r) for r in batch])
    scored = official.score(d045)

    # the official score cannot distinguish d045 from the two other zeros
    assert scored["official_score"] == 0.0
    # but the dense rank can, and puts it top of the zeros
    zeros = [1 if official.score(r)["official_score"] == 0.0 else 0 for r in batch]
    dense = [r["dense_rank"] for r in ranked]
    assert zeros[0] == 1
    assert dense[0] == max(d for d, z in zip(dense, zeros) if z)

    # and the two columns disagree, by design
    assert scored["limiting_term"] == "orient_oxy"      # the actual zero
    assert ranked[0]["worst_term"] == "ipde"            # batch-relative only
    assert scored["subscores"]["ipde"] > 0.9            # ipde needs no work


def test_worst_term_ignores_constant_columns():
    # Every term except orient_oxy is identical across this batch, so `base`
    # (weight 0.25) sits at a flat 0.5 percentile and would otherwise be named
    # as the thing to fix. A term that does not vary cannot be the
    # differentiator (LESSONS 2.5).
    raws = [_raw(dihedral12_oxy=290.0), _raw(), _raw()]
    ranked = O.rank_batch([O.badness(r) for r in raws])
    assert "base" in ranked[0]["constant_terms"]
    assert "orient_oxy" not in ranked[0]["constant_terms"]
    assert ranked[0]["worst_term"] == "orient_oxy"
    # the phantom shortfall is still visible for audit, just not acted on
    assert ranked[0]["weighted_shortfall"]["base"] > 0


def test_a_missing_term_costs_only_its_own_weight():
    good, missing = _raw(), _raw(ca_neighbors_mean=None)
    ranked = O.rank_batch([O.badness(good), O.badness(missing), O.badness(_raw())])
    assert ranked[1]["term_percentile"]["bury"] is None
    assert ranked[1]["weight_present"] == pytest.approx(1.0 - WEIGHTS["bury"])
    # it still has a usable rank rather than being nullified
    assert ranked[1]["dense_rank"] > 0


def test_dense_rank_is_bounded():
    raws = [_raw(), _raw(d_base=9.0), _raw(dihedral12_oxy=270.0)]
    for r in O.rank_batch([O.badness(x) for x in raws]):
        assert 0.0 <= r["dense_rank"] <= 1.0


def test_rank_batch_is_order_preserving_in_input():
    raws = [_raw(d_base=2.52), _raw(d_base=4.0)]
    ranked = O.rank_batch([O.badness(r) for r in raws])
    assert len(ranked) == 2
    assert ranked[0]["term_percentile"]["base"] == 1.0
    assert ranked[1]["term_percentile"]["base"] == 0.0


# ----------------------------------------------------------------- audit
def test_audit_fires_on_the_rfd3_inherited_geometry():
    # LESSONS 2.5: d_base 2.486 / d_oxy 3.338 identical to three decimals
    # across designs whose motifs landed at 167/114, 128/138 and 58/13 deg.
    rows = [{"d_base_A": 2.486, "d_oxy_A": 3.338, "dihedral_oxy_deg": d}
            for d in (114.0, 138.0, 13.0)]
    rep = O.audit(rows)
    assert not rep["passed"]
    assert set(rep["failures"]) == {"d_base_A", "d_oxy_A"}
    assert "constant column" in " ".join(rep["columns"]["d_base_A"]["reasons"])
    # the column that DOES vary must not be flagged
    assert "dihedral_oxy_deg" not in rep["failures"]


def test_audit_passes_a_healthy_batch():
    rows = [{"score": s, "d_base_A": d}
            for s, d in zip((0.1, 0.4, 0.7, 0.2), (2.5, 2.9, 3.1, 2.7))]
    assert O.audit(rows)["passed"]


def test_audit_fires_on_all_missing_and_on_mostly_missing():
    rows = [{"a": None, "b": 1.0}, {"a": None, "b": None},
            {"a": None, "b": None}, {"a": None, "b": 2.0}]
    rep = O.audit(rows, columns=["a", "b"])
    assert not rep["passed"]
    assert "a" in rep["failures"]            # all missing
    assert rep["columns"]["a"]["n_present"] == 0


def test_audit_fires_on_nan_only_column():
    rows = [{"a": float("nan")} for _ in range(4)]
    rep = O.audit(rows, columns=["a"])
    assert not rep["passed"]
    assert rep["columns"]["a"]["n_present"] == 0


def test_audit_fires_on_tiny_std():
    rows = [{"a": 1.0}, {"a": 1.0 + 1e-12}, {"a": 1.0 - 1e-12}, {"a": 1.0}]
    rep = O.audit(rows, columns=["a"])
    assert not rep["passed"]


def test_audit_skips_a_batch_too_small_to_judge_and_says_so():
    rep = O.audit([{"a": 1.0}, {"a": 1.0}])
    assert rep["skipped"] is True
    assert rep["passed"] is True
    assert "NOT a pass" in rep["reason"]


def test_audit_ignores_booleans_when_choosing_columns():
    rows = [{"eligible": True, "score": s} for s in (0.1, 0.2, 0.3)]
    rep = O.audit(rows)
    assert "eligible" not in rep["columns"]


def test_assert_ok_raises_with_the_offending_column_named():
    rows = [{"d_base_A": 2.486} for _ in range(4)]
    with pytest.raises(O.DegenerateScoresError, match="d_base_A"):
        O.assert_ok(rows, columns=["d_base_A"])


def test_assert_ok_returns_the_report_when_clean():
    rows = [{"a": float(i)} for i in range(4)]
    assert O.assert_ok(rows, columns=["a"])["passed"]


# ----------------------------------------------------------------- campaign checkpoint
def _bad(**kw) -> dict:
    base = dict(base=0.5, oxy=0.4, ipde=0.6, orient_base=10.0,
                orient_oxy=7.0, bury=20.0, cons=0.07)
    base.update(kw)
    return base


def test_checkpoint_refuses_to_judge_before_enough_designs():
    cp = O.campaign_checkpoint([_bad() for _ in range(4)], after=10)
    assert cp["judged"] is False
    assert cp["passed"] is True              # not a failure...
    assert "NOT a pass" in cp["reason"]      # ...but explicitly not a pass
    assert cp["n_scored"] == 4 and cp["needs"] == 10


def test_checkpoint_aborts_on_the_rfd3_inherited_geometry():
    """The LESSONS 2.5 failure: identical geometry across different designs.

    RFD3 reported d_base 2.486 / d_oxy 3.338 to three decimals across designs
    whose motifs landed at 167/114, 128/138 and 58/13 deg. Here `base` and
    `oxy` are copied from the input while everything else varies.
    """
    rows = [_bad(base=0.034, oxy=0.238,          # <- inherited from the input
                 ipde=0.5 + 0.01 * i, orient_base=float(i),
                 orient_oxy=float(20 - i), bury=15.0 + i, cons=0.01 * i)
            for i in range(10)]
    cp = O.campaign_checkpoint(rows, after=10)
    assert cp["judged"] is True
    assert cp["passed"] is False
    # exactly the two copied terms, and none of the five that vary
    assert set(cp["fatal"]) == {"base", "oxy"}
    assert "measuring its input" in cp["reason"]


def test_checkpoint_passes_a_varying_campaign():
    rows = [_bad(base=0.1 * i, oxy=0.05 * i, ipde=0.5 + 0.01 * i,
                 orient_base=float(i), orient_oxy=float(20 - i),
                 bury=15.0 + i, cons=0.01 * i) for i in range(10)]
    cp = O.campaign_checkpoint(rows, after=10)
    assert cp["passed"] is True and cp["judged"] is True
    assert cp["fatal"] == []


def test_a_mostly_missing_term_is_advisory_not_fatal():
    """Designs legitimately fail to fold; that must not kill a whole campaign."""
    rows = [_bad(cons=None if i % 3 else 0.01 * i, base=0.1 * i,
                 oxy=0.05 * i, ipde=0.5 + 0.01 * i, orient_base=float(i),
                 orient_oxy=float(20 - i), bury=15.0 + i) for i in range(10)]
    cp = O.campaign_checkpoint(rows, after=10)
    assert "cons" in cp["advisory"]
    assert "cons" not in cp["fatal"]
    assert cp["passed"] is True              # the sweep continues


def test_an_entirely_missing_term_is_fatal():
    rows = [_bad(oxy=None, base=0.1 * i, ipde=0.5 + 0.01 * i,
                 orient_base=float(i), orient_oxy=float(20 - i),
                 bury=15.0 + i, cons=0.01 * i) for i in range(10)]
    cp = O.campaign_checkpoint(rows, after=10)
    assert "oxy" in cp["fatal"]
    assert cp["passed"] is False


def test_checkpoint_can_be_disabled_by_a_zero_threshold():
    # `after=0` means every batch is judgeable; the sweep disables it by never
    # calling this (checkpoint_done starts True). Guard the arithmetic anyway.
    cp = O.campaign_checkpoint([_bad()], after=0)
    assert cp["judged"] is True


# ----------------------------------------------------------------- protocol
def test_protocol_stamp_records_the_msa_mode():
    s = O.protocol_stamp(replicates=5)
    assert s["msa_mode"] == "none:single_sequence"
    assert s["single_sequence_mode"] is True
    assert s["use_msa_server"] is False
    assert s["replicates"] == 5
    assert s["ligand_ccd"]
    assert "0.000" in s["comparability_warning"]


def test_protocol_stamp_carries_the_ramp_mode_and_targets():
    s = O.protocol_stamp()
    assert s["ramp_mode"] in ("taper", "flat")
    assert s["targets"]["d_base"] == 2.52
    assert s["targets"]["d_oxy"] == 3.10


def test_comparable_rejects_a_cross_protocol_comparison():
    design = O.protocol_stamp()
    control = {**O.protocol_stamp(), "msa_mode": "msa_server"}
    assert O.comparable(design, O.protocol_stamp())
    assert not O.comparable(design, control)


def test_comparable_rejects_a_ramp_mode_mismatch():
    a = O.protocol_stamp()
    assert not O.comparable(a, {**a, "ramp_mode": "flat"})
