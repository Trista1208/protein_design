"""Publish results to team-sharing/kmor/four-strategies, live.

    python scripts/live_publish_team.py [--push] [--loop SECONDS]

Collects every design we have measured - screen, campaign, finalists - into
one table, writes plots and structures, and commits to the team branch.
Designs are published whether they are good or bad: a 0.000 with a named
limiting term tells a teammate more than a missing row.

COLUMNS, and why each is there
------------------------------
Official (the competition metric, kemp/official.py):
    official_score, the seven sub-scores, limiting_term, sd_base
Measurement validity (docs/good_practices.md §2 - without these a number
cannot be interpreted):
    bound_ligand, mean_plddt, msa_mode, n_seeds
Literature activity filters (kemp/litfilters.py, from the supplied sheets):
    G1 distance, G2 angle, H3 charged-near-base, H2 base SASA, S2 apolar
Stability, as a GATE not a predictor - Risso et al. report that Tm and
Rosetta score do not correlate with activity on this reaction:
    rosetta_total_score, plddt

CHART RULES (dataviz skill): one panel per metric, single series, no dual
axes, thin marks, recessive grid, sequential single hue for magnitude, status
colour reserved for pass/fail and always paired with a text label.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
import json
import shutil
import subprocess
import time
from datetime import datetime
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C                                    # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TEAM = Path("/home/mandalorian/Hackathons/team-sharing")
DEST = TEAM / "keyshav-strategies/results"
BRANCH = "kmor/four-strategies"
IN_WINDOW = (2.0, 3.6)

#: dataviz: single sequential hue for magnitude; status reserved, with labels
INK = "#1c1c1c"
MUTED = "#6b6b6b"
GRID = "#e4e4e4"
SEQ = "#3b6fb6"
GOOD = "#2e7d4f"
BAD = "#b3493f"


def collect() -> list[dict]:
    rows: dict[str, dict] = {}

    def add(did, **kw):
        rows.setdefault(did, {"design_id": did}).update(
            {k: v for k, v in kw.items() if v is not None})

    # 1-seed screens
    # every screen, or the published table silently omits most of the run:
    # this read 107 designs while 298 had been folded, because library/ and
    # library2/ were never listed here
    for src, stage in ((ROOT / "outputs/insurance/screen.json", "screen"),
                       (ROOT / "outputs/library/screen.json", "library1"),
                       (ROOT / "outputs/library2/screen.json", "library2"),
                       (ROOT / "outputs/campaign/campaign.json", "campaign"),
                       (ROOT / "outputs/gen2/screen.json", "gen2")):
        if not src.exists():
            continue
        try:
            data = json.loads(src.read_text())
        except Exception:                                       # noqa: BLE001
            continue
        for r in data:
            if r.get("status") not in (None, "ok"):
                continue
            db, do = r.get("d_base"), r.get("d_oxy")
            add(r["design_id"], stage=stage, n_seeds=1,
                sequence=r.get("sequence"),
                mutations="+".join(r.get("mutations", []) or
                                   r.get("edits", []) or []),
                d_base=db, d_oxy=do, dihedral_oxy=r.get("dih_oxy"),
                gap_deg=(r.get("gap") if r.get("gap") is not None
                         else r.get("gap_dih_oxy")),
                mean_plddt=r.get("plddt"), ipde=r.get("ipde"),
                bound_ligand=(db is not None and do is not None
                              and IN_WINDOW[0] <= db <= IN_WINDOW[1]
                              and IN_WINDOW[0] <= do <= IN_WINDOW[1]),
                pinner_A=((r.get("pinner") or {}).get("distance_A")),
                pinner_residue=((r.get("pinner") or {}).get("residue")),
                **{f"lit_{k}": v for k, v in
                   (r.get("lit", {}).get("values", {}) or {}).items()
                   if not isinstance(v, (list, dict))})

    # 5-seed official scores override the 1-seed row
    fin = ROOT / "outputs/finalists"
    if fin.exists():
        for vd in fin.glob("*/verdict.json"):
            try:
                v = json.loads(vd.read_text())
            except Exception:                                   # noqa: BLE001
                continue
            s = v.get("subscores", {})
            add(vd.parent.name, stage="finalist_5seed",
                n_seeds=len(v.get("replicates", [])) or C.REPLICATES,
                sequence=v.get("sequence"),
                official_score=v.get("official_score"),
                official_best=v.get("official_score_best"),
                n_nonzero_replicates=v.get("n_nonzero_replicates"),
                per_replicate_total=str(v.get("per_replicate_total") or ""),
                limiting_term=v.get("limiting_term"),
                sd_base=v.get("sd_base"), mean_plddt=v.get("mean_plddt"),
                msa_mode=(v.get("protocol") or {}).get("msa_mode"),
                **{f"s_{k}": val for k, val in s.items()})
    return sorted(rows.values(),
                  key=lambda r: (-(r.get("official_score") or 0),
                                 r.get("gap_deg") if r.get("gap_deg") is not None
                                 else 999))


COLS = ["design_id", "stage", "n_seeds", "official_best",
        "n_nonzero_replicates", "official_score", "per_replicate_total",
        "limiting_term",
        "s_base", "s_oxy", "s_ipde", "s_orient_base", "s_orient_oxy",
        "s_bury", "s_cons", "sd_base", "bound_ligand", "mean_plddt",
        "d_base", "d_oxy", "dihedral_oxy", "gap_deg", "ipde",
        "lit_G1_distance", "lit_G2_angle", "lit_H3_sequestering",
        "lit_S2_apolar_fraction", "lit_S7_donor_contacts", "pinner_A", "pinner_residue",
        "msa_mode", "mutations", "sequence"]


def write_table(rows, dest: Path):
    dest.mkdir(parents=True, exist_ok=True)
    with (dest / "designs.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_markdown(rows, dest: Path, meta: dict):
    """Landing page: plots first, best designs, then what is failing and why.

    GitHub renders this file on the directory page, so the plots must be
    EMBEDDED, not merely present as files. And the headline number must be
    `best` with its reproducibility - the median is 0.000 for a design that
    scores on one seed in five, which hides the only real result we have.
    """
    scored = [r for r in rows if r.get("n_seeds", 0) >= 5]
    hits = sorted((r for r in scored if (r.get("official_best") or 0) > 0),
                  key=lambda r: -(r["official_best"]))
    zeros = [r for r in scored if (r.get("official_best") or 0) <= 0]
    bound = [r for r in rows if r.get("bound_ligand")]
    best = hits[0]["official_best"] if hits else 0.0

    L = ["# Kemp eliminase designs — live results", "",
         f"_Auto-updated every 10 min. Last: {meta['updated']}_", "",
         "| | |", "|---|---|",
         f"| **Best official score** | **{best:.4f}**"
         + (f" ({hits[0]['n_nonzero_replicates']}/5 seeds)" if hits else "")
         + " |",
         f"| Designs measured | {len(rows)} |",
         f"| With a bound ligand | {len(bound)} |",
         f"| Fully scored (5 seeds) | {len(scored)} |",
         f"| Scoring on ≥1 replicate | **{len(hits)}** |",
         "| HG3.17 crystal (yardstick) | 0.959 |",
         "| HG3.17 folded with MSA | 0.834 |", ""]

    plots = sorted((dest / "plots").glob("*.png")) if (dest / "plots").exists() else []
    wanted = [p for p in plots
              if p.name in ("official_score.png", "dihedral_gap.png",
                            "distances.png", "progress.png")]
    if wanted:
        L += ["## Plots", ""]
        for p in wanted:
            L += [f"### {p.stem.replace('_', ' ').title()}", "",
                  f"![{p.stem}](plots/{p.name})", ""]

    L += ["## How to read a score here", "",
          "The official metric is a weighted **geometric** mean: one zero "
          "sub-score makes the whole design 0.000. Every design is scored on "
          "**5 Boltz-2 seeds, each scored separately**, then aggregated — "
          "never by averaging the geometry and scoring once. Averaging a "
          "bimodal rotamer produced a phantom 0.624 here for a design where "
          "0/5 replicates scored.", "",
          "So a result is quoted as **best achievable + reproducibility**: "
          "\"0.743 on 1/5 seeds\", never a bare average.", "",
          "`bound_ligand` is the first filter — geometry is meaningless "
          "unless both `d_base` and `d_oxy` sit inside 2.0–3.6 Å. "
          "`mean_plddt` below 0.95 means *unmeasurable*, not bad.", ""]

    if hits:
        L += ["## ✅ Designs that score", "",
              "| design | best | seeds | per-replicate | limiting | pLDDT |",
              "|---|---|---|---|---|---|"]
        for r in hits[:20]:
            L.append(f"| `{r['design_id']}` | **{r['official_best']:.4f}** | "
                     f"{r.get('n_nonzero_replicates')}/5 | "
                     f"`{r.get('per_replicate_total','')}` | "
                     f"{r.get('limiting_term','—')} | "
                     f"{(r.get('mean_plddt') or 0):.3f} |")
        L.append("")

    gdir = dest / "gifs"
    gifs = sorted(gdir.glob("*.gif")) if gdir.exists() else []
    if gifs:
        L += ["### Rotating 3D views", "",
              "Cartoon protein, **6NT** analogue in salmon sticks, catalytic "
              "**Asp102 / Gln46** in blue, the **Tyr42 / Tyr103** pinners in "
              "green, and the two scored contacts as dashed lines. Rendered "
              "from the seed that scored.", ""]
        for g in gifs[:6]:
            did = g.stem
            r = next((x for x in hits if x["design_id"] == did), None)
            cap = (f"{did} — best {r['official_best']:.4f}, "
                   f"{r.get('n_nonzero_replicates')}/5 seeds" if r else did)
            L += [f"**{cap}**", "", f"![{did}](gifs/{g.name})", ""]

    if zeros:
        cnt = Counter(r.get("limiting_term") or "?" for r in zeros)
        L += ["## ❌ Designs that score 0.000 — and which term killed them", "",
              "| limiting term | designs |", "|---|---|"]
        for k, n in cnt.most_common():
            L.append(f"| `{k}` | {n} |")
        L += ["", "| design | limiting | base | oxy | ipde | or_base | "
              "or_oxy | bury | cons | pLDDT |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for r in zeros[:20]:
            g = lambda k: ("—" if r.get(k) is None else f"{r[k]:.3f}")  # noqa: E731
            L.append(f"| `{r['design_id']}` | {r.get('limiting_term','—')} | "
                     f"{g('s_base')} | {g('s_oxy')} | {g('s_ipde')} | "
                     f"{g('s_orient_base')} | {g('s_orient_oxy')} | "
                     f"{g('s_bury')} | {g('s_cons')} | {g('mean_plddt')} |")
        L.append("")

    lit = [r for r in bound if r.get("lit_G1_distance") is not None]
    if lit:
        L += ["## Literature activity filters", "",
              "From the supplied sheets (Risso 2020; Privett 2012; "
              "Alexandrova 2008), all for this reaction. Reference: G1 "
              "2.64 Å, G2 166.5°, H3 = 0 charged residues within 4 Å of the "
              "base, S2 ≥ 0.60.", "",
              "| design | G1 Å | G2 ° | H3 | S2 |", "|---|---|---|---|---|"]
        for r in lit[:15]:
            h3 = r.get("lit_H3_sequestering")
            L.append(f"| `{r['design_id']}` | {r.get('lit_G1_distance')} | "
                     f"{r.get('lit_G2_angle')} | "
                     f"{'0 ✅' if h3 == 0 else f'{h3} ⚠️'} | "
                     f"{r.get('lit_S2_apolar_fraction')} |")
        L.append("")

    L += biophysics_section()

    L += ["## Correction — if you saw an earlier commit", "",
          "Commit `b93b0608` reported **best S = 0.624098**. That is "
          "**retracted**: it averaged the geometry across seeds and scored "
          "once, and the design's dihedral is bimodal (~99° / ~245°, never "
          "between). The mean of those clusters lands inside the allowed "
          "band that no replicate occupies. Per replicate, 0/5 score.", "",
          "## Files", "",
          "- `designs.csv` — every measured design, all metrics",
          "- `structures/` — folded PDBs (all 5 seeds for designs that score)",
          "- `plots/` — the panels above",
          "- `status.json` — machine-readable summary", ""]
    (dest / "README.md").write_text("\n".join(L))


def scoring_section() -> list[str]:
    """How a number on this page was produced. Read from kemp/official.py so
    the page cannot drift from the code that produced the scores."""
    from kemp import official as O
    t = O.TARGETS
    d180 = O.LAUKO["dihedral_sp2_180"]
    d0 = O.LAUKO["dihedral_sp2_0"]
    return [
        "## How each score is computed", "",
        "The official metric is a **weighted geometric mean** of seven "
        "subscores:", "",
        "```",
        "S = PROD_i  s_i ** w_i        equivalently  S = exp( SUM_i w_i * ln s_i )",
        "```", "",
        "**One zero kills the design.** If any `s_i` is 0 then S = 0 no matter "
        "how good the other six are. That single property shaped the whole "
        "pipeline: a batch of all-zero designs is unrankable, so we optimise "
        "against a *dense* surrogate (degrees outside the allowed dihedral "
        "band) and score with the real metric.", "",
        "| term | measured quantity | formula | target / window |",
        "|---|---|---|---|",
        f"| `base` | `d_base` — Asp102 OD2 ⋯ ligand N3 | window score | "
        f"target **{t['d_base']['target']} Å**, window "
        f"{t['d_base']['window'][0]}–{t['d_base']['window'][1]} Å |",
        f"| `ipde` | Boltz-2 `complex_ipde` (an *error*) | "
        f"`max(0, 1 − iPDE/8.0)` | lower is better |",
        f"| `oxy` | `d_oxy` — Gln46 NE2 ⋯ ligand N1 | window score | "
        f"target **{t['d_oxy']['target']} Å**, window "
        f"{t['d_oxy']['window'][0]}–{t['d_oxy']['window'][1]} Å |",
        "| `orient_base` | angle + dihedral at the base | "
        "`min(s_angle, s_dihedral)` | — |",
        "| `orient_oxy` | angle + dihedral at the oxyanion | "
        "`min(s_angle, s_dihedral)` | — |",
        f"| `bury` | mean Cα neighbours at "
        f"{'/'.join(str(int(c)) for c in O.CA_DIST_CUTOFFS)} Å | linear ramp | "
        f"0 at ≤{int(O.BURY_TARGET['poor'])}, 1 at ≥{int(O.BURY_TARGET['good'])} |",
        f"| `cons` | SD of `d_base` across the 5 seeds | "
        f"`max(0, 1 − σ/{O.CONS_TOLERANCE_A})` | lower is better |",
        "",
        "### Window score (distances and angles)", "",
        "```",
        "s = 0                        if v outside [lo, hi]",
        "s = 1 − |v − t| / span       otherwise",
        "span = (hi − t) if v ≥ t else (t − lo)",
        "```", "",
        f"1.0 at the target, linear decay to 0 at whichever window edge is "
        f"nearer. Angles use target **{t['angle']['target']}°**, window "
        f"**{t['angle']['window'][0]}–{t['angle']['window'][1]}°**.", "",
        "### Dihedral score — the term that decided everything", "",
        f"Two allowed sp2 bands, **{int(d180[0])}–{int(d180[1])}°** and "
        f"**{int(d0[0])}–{int(d0[1])}°** (wrapping through 0). Ideal "
        f"{int(O.DIHEDRAL_IDEAL)}° (or 0°), half-width "
        f"{int(O.DIHEDRAL_HALF_WIDTH)}°:", "",
        "```",
        "s_dih = max(0, 1 − |Δ| / 60)     Δ = deviation from the nearer band's ideal",
        "```", "",
        f"Everything in **{int(d180[0])-80}–{int(d180[0])}°** and "
        f"**{int(d180[1])}–{int(d0[0])}°** scores **exactly 0**. The team's "
        f"`rfd4_design_045` sat at **86.1°** — dead in the first gap — which "
        f"is why it scored 0.000 on all five seeds despite good distances. "
        f"Three mutations moved it into band.", "",
        "`orient_* = min(angle, dihedral)`, never an average: per "
        "SUBMISSION.md *\"the angle and the in-plane terms do not average — "
        "the worse of the two decides.\"* A perfect angle cannot rescue a bad "
        "torsion.", "",
        "### How a design on this page got its number", "",
        "1. Fold the sequence with the 6NT ligand in Boltz-2 — **5 seeds**.",
        "2. Measure every quantity **from each folded structure**, never from "
        "stored numbers.",
        "3. Gate: both distances inside 2.0–3.6 Å, mean pLDDT ≥ 0.95.",
        "4. **Score each replicate separately**, then aggregate the *scores*.",
        "5. Report median, best, worst and **reproducibility = n/5 non-zero**.",
        "", "---", ""]


def biophysics_section() -> list[str]:
    """Biochemistry / biophysics / thermodynamics, inline on the landing page.

    Read from published/biophysics.json, which biophysics_table.py regenerates
    earlier in this same publish, so these numbers can never be stale relative
    to the scores above them. Returns [] if that file is missing rather than
    inventing a section - an absent table is better than a wrong one.
    """
    f = ROOT / "published/biophysics.json"
    if not f.exists():
        return []
    try:
        rows = json.loads(f.read_text())
    except Exception:                                            # noqa: BLE001
        return []
    if not rows:
        return []

    def g(r, k, n=3, sign=False):
        v = r.get(k)
        if v is None:
            return "–"
        return f"{v:+.{n}f}" if sign else f"{v:.{n}f}"

    ref = next((r for r in rows if r["design_id"].startswith("REF_")), None)
    designs = [r for r in rows if not r["design_id"].startswith("REF_")]

    L = scoring_section()

    L += ["## Stability — biochemistry, biophysics, thermodynamics", "",
         "Every design here is folded, not measured in a tube. These columns "
         "ask a different question from the score above: *is this protein "
         "buildable at all?*", ""]
    if ref:
        L += [f"The **HG3.17 (4BS0) crystal structure is carried through as "
              f"the control arm** — a real enzyme that has been expressed, "
              f"purified, crystallised and kinetically characterised. It "
              f"relaxes to **{g(ref,'relax_per_res')} REU/residue**, so every "
              f"design below that beats it is better packed than the enzyme "
              f"this reaction is actually run with.", ""]
    # per-design verdict: the competition's own gates, then buildability
    if any(r.get("stability_verdict") for r in designs):
        n_clean = sum(1 for r in designs if r.get("stability_verdict") == "CLEAN")
        n_flag = sum(1 for r in designs if r.get("stability_verdict") == "FLAGGED")
        n_elig = sum(1 for r in designs if r.get("eligibility") == "NOVEL")
        L += ["### Verdict per design", "",
              "Eligibility is the competition's rule set. Stability is **not** "
              "a competition criterion — it is a buildability check and was "
              "never used to rank anything.", "",
              "A design is CLEAN when all of these hold: ThermoMPNN "
              "ΔΔG ≤ 1.5 kcal/mol (~3× the ~0.5 experimental error), "
              "instability index < 40, Rosetta relax per residue better than "
              "the HG3.17 control, MD drift ≤ 0.5 Å/ns, Rg swell ≤ 5%. "
              "FLAGGED names the term that failed rather than hiding it.", "",
              "| design | repro | official | eligible | TM | identity | "
              "base/oxy | stability | why flagged |",
              "|---|---|---|---|---|---|---|---|---|"]
        for r in designs:
            tm, idn = r.get("tm_max"), r.get("identity_max")
            sv = r.get("stability_verdict")
            mark = {"CLEAN": "✅ CLEAN", "FLAGGED": "⚠️ FLAGGED"}.get(sv, "–")
            el = r.get("eligibility")
            L.append(f"| `{r['design_id']}` | {r.get('repro','–')} | "
                     f"{g(r,'official_median',4)} | "
                     f"{'✅ NOVEL' if el == 'NOVEL' else (el or '–')} | "
                     f"{'–' if tm is None else f'{tm:.3f}'} | "
                     f"{'–' if idn is None else f'{idn:.2f}'} | "
                     f"{r.get('base_idx','–')}/{r.get('oxy_idx','–')} | "
                     f"{mark} | {r.get('stability_reasons') or '–'} |")
        L += ["", f"**{n_elig} eligible, {n_clean} clean, {n_flag} flagged.** "
              f"A flagged design is still submitted: SUBMISSION.md ranks each "
              f"design individually against other teams' entries, so a weaker "
              f"one cannot drag the others down, and an occupied slot beats an "
              f"empty one.", ""]

    # thermal MD, when the run has reached these designs
    md = [r for r in designs + ([ref] if ref else [])
          if r and r.get("md_rmsd_A") is not None]
    if md:
        t = md[0]
        L += ["### Thermal — molecular dynamics", "",
              f"{t['md_ns']:g} ns implicit-solvent (GBn2) MD at "
              f"{t['md_temp_K']:g} K, protein only, identical settings for "
              f"every run including the control. A folded protein's RMSD "
              f"plateaus (drift ≈ 0); an unfolding one keeps climbing, and "
              f"its radius of gyration swells.", "",
              "| design | RMSD Å | max Å | drift Å/ns | Rg start → end Å | swell % |",
              "|---|---|---|---|---|---|"]
        for r in md:
            rg0, rg1 = r.get("md_rg_start_A"), r.get("md_rg_end_A")
            sw = f"{100*(rg1-rg0)/rg0:+.1f}" if rg0 and rg1 else "–"
            tag = (f"**`{r['design_id']}`** (control)"
                   if r["design_id"].startswith("REF_") else f"`{r['design_id']}`")
            L.append(f"| {tag} | {r['md_rmsd_A']:.2f} | "
                     f"{g(r,'md_rmsd_max_A',2)} | {r['md_drift_A_per_ns']:+.2f} | "
                     f"{rg0:.2f} → {rg1:.2f} | {sw} |")
        L += ["", "**Not a melting temperature.** Two nanoseconds of implicit "
              "solvent cannot unfold a stable protein, so a flat trace is weak "
              "evidence of stability while a rising one is strong evidence of "
              "a problem. Only the comparison between identically-run "
              "simulations means anything.", ""]

    L += ["### Thermodynamics", "",
          "`relax/res` — Rosetta relax total_score per residue, protein only; "
          "**lower is better packed**.",
          "`ΔΔG` — ThermoMPNN, summed over each design's mutations relative to "
          "the parent shell. **Positive = destabilising.** The sign was "
          "verified empirically, not assumed: a self-mutation scores exactly "
          "0.000, →Pro +2.06, →Gly +1.36, hydrophobic→charged +2.06. "
          "Experimental error on a measured ΔΔG is ~0.5 kcal/mol, so anything "
          "under that is indistinguishable from zero.", "",
          "| design | repro | official | relax/res | ΔΔG kcal/mol | mutations |",
          "|---|---|---|---|---|---|"]
    for r in designs:
        L.append(f"| `{r['design_id']}` | {r.get('repro','–')} | "
                 f"{g(r,'official_median',4)} | **{g(r,'relax_per_res')}** | "
                 f"{g(r,'thermompnn_ddg_sum',2,sign=True)} | "
                 f"{r.get('mutations','–')} |")
    if ref:
        L.append(f"| **`{ref['design_id']}`** (control) | – | – | "
                 f"**{g(ref,'relax_per_res')}** | – | reference enzyme |")

    L += ["", "### Biophysics — fold and surface", "",
          "`core` — buried hydrophobics (relative SASA < 0.15) ÷ length. "
          "`exp_hyd` / `exp_chg` — residues with relative SASA > 0.40. High "
          "exposed charge favours soluble expression.", "",
          "| design | core | exp_hyd % | exp_chg % | SASA Å² | pLDDT | iPDE |",
          "|---|---|---|---|---|---|---|"]
    for r in designs + ([ref] if ref else []):
        tag = f"**`{r['design_id']}`** (control)" if r is ref else f"`{r['design_id']}`"
        L.append(f"| {tag} | {g(r,'core_frac')} | {g(r,'exposed_hyd_pct',1)} | "
                 f"{g(r,'exposed_chg_pct',1)} | {g(r,'sasa_total',0)} | "
                 f"{g(r,'plddt')} | {g(r,'ipde',2)} |")

    L += ["", "### Biochemistry — sequence", "",
          "Instability index < 40 is the classical Guruprasad threshold for a "
          "protein predicted stable *in vivo*. GRAVY below zero is "
          "hydrophilic.", "",
          "| design | len | pI | net charge | GRAVY | instability | arom % | Cys |",
          "|---|---|---|---|---|---|---|---|"]
    for r in designs + ([ref] if ref else []):
        tag = f"**`{r['design_id']}`** (control)" if r is ref else f"`{r['design_id']}`"
        nc = r.get("net_charge")
        L.append(f"| {tag} | {r.get('length','–')} | {g(r,'pI',2)} | "
                 f"{'–' if nc is None else f'{nc:+d}'} | {g(r,'gravy',2)} | "
                 f"{g(r,'instability_index',1)} | {g(r,'aromatic_pct',1)} | "
                 f"{r.get('n_cys','–')} |")

    L += ["", "**What this does and does not say.** It says nothing here "
          "indicates a folding, packing or expression problem. It does **not** "
          "say these proteins are stable: Rosetta score is a packing proxy and "
          "ThermoMPNN is a sequence–structure predictor, and neither is a "
          "melting temperature. The supplied Risso 2020 sheet is explicit that "
          "Tm is *a gate, not a predictor* of activity — so none of this was "
          "used to rank designs, only to check that the ranked ones are "
          "buildable.", "",
          "Full tables, including every column: [`BIOPHYSICS.md`](BIOPHYSICS.md) "
          "· [`biophysics.csv`](biophysics.csv)", ""]
    return L


def make_plots(rows, dest: Path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:                                           # noqa: BLE001
        return []
    dest.mkdir(parents=True, exist_ok=True)
    made = []
    bound = [r for r in rows if r.get("bound_ligand")
             and r.get("gap_deg") is not None]

    def style(ax):
        ax.spines[["top", "right"]].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=8, length=3)
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.set_axisbelow(True)

    # one panel per metric, single series each - never a shared axis
    if bound:
        fig, ax = plt.subplots(figsize=(7, 3.2), dpi=160)
        g = sorted(r["gap_deg"] for r in bound)
        ax.bar(range(len(g)), g, color=SEQ, width=0.8)
        ax.axhline(0, color=GOOD, lw=1.2)
        ax.set_title("Oxyanion dihedral: degrees outside the allowed band",
                     color=INK, fontsize=10, loc="left")
        ax.set_ylabel("gap (°)  — 0 is in band", color=MUTED, fontsize=8)
        ax.set_xlabel("bound-ligand designs, sorted", color=MUTED, fontsize=8)
        style(ax)
        fig.tight_layout(); fig.savefig(dest / "dihedral_gap.png"); plt.close(fig)
        made.append("dihedral_gap.png")

        fig, ax = plt.subplots(figsize=(7, 3.2), dpi=160)
        xs = [r["d_base"] for r in bound if r.get("d_base")]
        ys = [r["d_oxy"] for r in bound if r.get("d_oxy")]
        ax.scatter(xs, ys, s=26, color=SEQ, alpha=0.85, edgecolor="white", lw=0.6)
        ax.scatter([2.532], [3.102], s=70, color=GOOD, marker="*",
                   edgecolor="white", lw=0.8, label="4BS0 crystal")
        ax.axvspan(2.0, 3.6, color=GRID, alpha=0.35, zorder=0)
        ax.set_title("Catalytic distances (bound-ligand designs only)",
                     color=INK, fontsize=10, loc="left")
        ax.set_xlabel("d_base  Å (target 2.52)", color=MUTED, fontsize=8)
        ax.set_ylabel("d_oxy  Å (target 3.10)", color=MUTED, fontsize=8)
        ax.legend(frameon=False, fontsize=8, labelcolor=MUTED)
        style(ax)
        fig.tight_layout(); fig.savefig(dest / "distances.png"); plt.close(fig)
        made.append("distances.png")

    scored = [r for r in rows if r.get("n_seeds", 0) >= 5]
    if scored:
        fig, ax = plt.subplots(figsize=(7, 3.2), dpi=160)
        # BEST achievable per design, not the median: a design scoring on one
        # seed in five has median 0.000, which hides the only real result.
        s = sorted((r.get("official_best") or 0.0) for r in scored)[::-1]
        cols = [GOOD if x > 0 else BAD for x in s]
        ax.bar(range(len(s)), s, color=cols, width=0.8)
        ax.axhline(0.959, color=MUTED, lw=1, ls="--")
        ax.text(0, 0.965, "HG3.17 crystal 0.959", color=MUTED, fontsize=7)
        ax.set_title("Best official score per design (5 seeds, scored per replicate)",
                     color=INK, fontsize=10, loc="left")
        ax.set_ylabel("weighted geometric mean", color=MUTED, fontsize=8)
        style(ax)
        fig.tight_layout(); fig.savefig(dest / "official_score.png"); plt.close(fig)
        made.append("official_score.png")
    return made


def stage_structures(rows, dest: Path, limit: int = 140):
    """One PDB per measured design; ALL five seeds for anything that scores.

    475 PDBs (95 designs x 5 seeds) is ~57 MB, too much for a branch people
    have to clone. A design that scores on at least one replicate is worth
    its full set - the per-seed spread IS the result there - and everything
    else needs only the structure the metrics were measured on.
    """
    dest.mkdir(parents=True, exist_ok=True)
    n = 0
    for r in rows:
        did = r["design_id"]
        # a scoring design gets every seed, regardless of the limit
        if (r.get("official_score_best") or 0) > 0:
            src = ROOT / "outputs/finalists" / did
            for sd in sorted(src.glob("seed*/folded.pdb")):
                shutil.copy2(sd, dest / f"{did}__{sd.parent.name}.pdb")
                n += 1
            continue
        if n >= limit:
            continue
        for cand in (ROOT / "outputs/finalists" / did / "seed1" / "folded.pdb",
                     ROOT / "outputs/insurance/folds" / did / "folded.pdb"):
            if cand.exists():
                shutil.copy2(cand, dest / f"{did}.pdb")
                n += 1
                break
        else:
            for c in (ROOT / "outputs/campaign").rglob(f"{did}/folded.pdb"):
                shutil.copy2(c, dest / f"{did}.pdb")
                n += 1
                break
    return n


def publish(push: bool) -> dict:
    rows = collect()
    meta = {"updated": datetime.now().strftime("%F %T"), "n_designs": len(rows)}
    write_table(rows, DEST)
    write_markdown(rows, DEST, meta)
    plots = make_plots(rows, DEST / "plots")
    nstruct = stage_structures(rows, DEST / "structures")
    # biochemistry / biophysics / thermodynamics for every scored design,
    # carrying the HG3.17 crystal through as the control arm. Regenerated on
    # every publish so the branch never shows stale stability numbers beside
    # fresh scores. Rosetta relax is NOT re-run here (minutes per structure);
    # scripts/stability_top20.py caches it and this step reads that cache.
    for step in ("scripts/make_designs_fasta.py", "scripts/biophysics_table.py",
                 "scripts/archive_results.py"):
        try:
            subprocess.run([sys.executable, str(ROOT / step)],
                           capture_output=True, text=True, timeout=1200,
                           env={**__import__("os").environ,
                                "KEMP_ROOT": str(ROOT)})
        except Exception as exc:                                # noqa: BLE001
            print(f"  {step}: {type(exc).__name__}")
    for name in ("BIOPHYSICS.md", "biophysics.csv", "biophysics.json",
                 "ARCHIVE.md", "archive.json"):
        src_f = ROOT / "published" / name
        if src_f.exists():
            shutil.copy2(src_f, DEST / name)
    for name in ("designs.fasta", "designs_manifest.json"):
        src_f = ROOT / "published/submission" / name
        if src_f.exists():
            (DEST / "submission").mkdir(parents=True, exist_ok=True)
            shutil.copy2(src_f, DEST / "submission" / name)

    # rotating 3D GIFs for anything scoring on >=1 replicate, including
    # partial passes - run_monitor.md 4: a 1/5 pass is a real result
    ngif = 0
    try:
        subprocess.run([sys.executable, str(ROOT / "scripts/make_gifs.py"),
                        "--max", "10", "--frames", "30"],
                       capture_output=True, text=True, timeout=1800,
                       env={**__import__("os").environ, "KEMP_ROOT": str(ROOT)})
        src = ROOT / "outputs/gifs"
        if src.exists():
            (DEST / "gifs").mkdir(parents=True, exist_ok=True)
            for g in src.glob("*.gif"):
                shutil.copy2(g, DEST / "gifs" / g.name)
                ngif += 1
            if (src / "gifs.json").exists():
                shutil.copy2(src / "gifs.json", DEST / "gifs" / "gifs.json")
    except Exception as exc:                                    # noqa: BLE001
        print(f"  gif step: {type(exc).__name__}")
    scored = [r for r in rows if r.get("official_score") is not None]
    nonzero = [r for r in scored if (r["official_score"] or 0) > 0]
    status = {**meta, "n_bound": sum(1 for r in rows if r.get("bound_ligand")),
              "n_scored_5seed": len(scored), "n_nonzero": len(nonzero),
              # BEST achievable, not the median. A design scoring on 1 of 5
              # seeds has median 0.000, so the commit message read
              # "best S=0" while the README correctly showed 0.7435.
              "best_official": max((r.get("official_best") or 0
                                    for r in scored), default=0),
              "best_reproducibility": next(
                  (f"{r.get('n_nonzero_replicates')}/5" for r in
                   sorted(scored, key=lambda x: -(x.get("official_best") or 0))
                   if (r.get("official_best") or 0) > 0), None),
              "n_structures": nstruct, "n_gifs": ngif, "plots": plots}
    (DEST / "status.json").write_text(json.dumps(status, indent=2))

    msg = (f"live: {len(rows)} designs, {status['n_bound']} bound, "
           f"{len(nonzero)} non-zero, best S={status['best_official']} "
           f"({meta['updated']})")
    subprocess.run(["git", "add", "-A", "keyshav-strategies/results"],
                   cwd=TEAM, check=False)
    c = subprocess.run(["git", "commit", "-m", msg], cwd=TEAM,
                       capture_output=True, text=True)
    if c.returncode == 0 and push:
        p = subprocess.run(["git", "push", "origin", BRANCH], cwd=TEAM,
                           capture_output=True, text=True, timeout=180)
        status["pushed"] = p.returncode == 0
        if p.returncode != 0:
            status["push_error"] = p.stderr[-200:]
    print(f"  {msg}  pushed={status.get('pushed')}")
    return status


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--push", action="store_true")
    ap.add_argument("--loop", type=int, default=0,
                    help="republish every N seconds")
    ap.add_argument("--until", default=None, help="stop time HH:MM")
    a = ap.parse_args()
    stop = None
    if a.until:
        hh, mm = (int(x) for x in a.until.split(":"))
        stop = datetime.now().replace(hour=hh, minute=mm, second=0).timestamp()
    while True:
        try:
            publish(a.push)
        except Exception as exc:                                # noqa: BLE001
            print(f"  publish error: {type(exc).__name__}: {exc}")
        if not a.loop or (stop and time.time() > stop):
            return 0
        time.sleep(a.loop)


if __name__ == "__main__":
    raise SystemExit(main())
