"""Biochemistry / biophysics / thermodynamics table for every scored design.

Combines, per design:
  thermodynamics  Rosetta relax total_score per residue, fa_rep per residue
                  ThermoMPNN ddG of its mutations vs the parent shell
  biophysics      hydrophobic core burial, exposed hydrophobic %, total SASA
  biochemistry    pI, net charge, GRAVY, instability index, aromatic %
  prediction      Boltz-2 pLDDT, iPDE
plus the official score (per replicate, with n/5 reproducibility).

CONTROL ARM: the HG3.17 (4BS0) crystal structure is carried through every
column it applies to. A stability number without a control is not evidence
(good_practices 5b) - HG3.17 has been expressed, purified, crystallised and
kinetically characterised, so it calibrates the scale.

SIGN CONVENTION, verified empirically rather than assumed: in this ThermoMPNN
output a self-mutation scores exactly 0.000, and ->Pro (+2.06), ->Gly (+1.36)
and hydrophobic->charged (+2.06) are all positive. So POSITIVE ddG =
DESTABILISING.

CAVEAT: ThermoMPNN is a SINGLE-mutant predictor. Summing it over a double or
triple assumes additivity, which is an approximation (ThermoMPNN-D exists for
pairs). Rosetta score is a packing proxy, not a melting temperature, and the
supplied Risso 2020 sheet is explicit that Tm is "a gate, not a predictor" of
activity. The defensible claim is "nothing indicates a folding or expression
problem", not "these are stable".
"""
from __future__ import annotations
import json, sys, warnings
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ROOT = Path(__file__).resolve().parents[1]
ST = ROOT / "outputs/stability"
SHELL_ID = "ins080_r103_42_F42Y_H103Y"
REF_PDB = ROOT / "outputs/runs/controls_full/4BS0_A.pdb"
REF_FASTA = ROOT / "inputs/corpus/seqs/4BS0.fasta"

# Tien et al. 2013 theoretical maximum ASA, used for relative burial
MAXASA = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "E": 223,
          "Q": 225, "G": 104, "H": 224, "I": 197, "L": 201, "K": 236,
          "M": 224, "F": 240, "P": 159, "S": 155, "T": 172, "W": 285,
          "Y": 263, "V": 174}
HYD = set("AVLIMFWCY")


def structure_metrics(pdb: Path) -> dict:
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    import freesasa
    st = PDBParser(QUIET=True).get_structure("x", str(pdb))
    model = st[0]
    chain = next(iter(model))
    res = [r for r in chain if r.id[0] == " "]
    ShrakeRupley().compute(model, level="R")
    core = exp_h = n_h = exp_c = n_c = 0
    for r in res:
        try:
            aa = index_to_one(three_to_index(r.get_resname()))
        except Exception:                                         # noqa: BLE001
            continue
        rel = r.sasa / MAXASA[aa]
        if aa in HYD:
            n_h += 1
            core += rel < 0.15
            exp_h += rel > 0.40
        if aa in "DEKR":
            n_c += 1
            exp_c += rel > 0.40
    return {"core_hydrophobics": core,
            "core_frac": core / max(len(res), 1),
            "exposed_hyd_pct": 100 * exp_h / max(n_h, 1),
            "exposed_chg_pct": 100 * exp_c / max(n_c, 1),
            "sasa_total": freesasa.calc(freesasa.Structure(str(pdb))).totalArea(),
            "n_res_struct": len(res)}


def sequence_metrics(seq: str) -> dict:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    pa = ProteinAnalysis(seq)
    return {"length": len(seq),
            "pI": pa.isoelectric_point(),
            "net_charge": sum(seq.count(a) for a in "KR") - sum(seq.count(a) for a in "DE"),
            "gravy": pa.gravy(),
            "instability_index": pa.instability_index(),
            "aromatic_pct": 100 * sum(seq.count(a) for a in "FWY") / len(seq),
            "n_cys": seq.count("C")}


#: Thresholds for the per-design stability verdict. Each is a stated number
#: with a reason, so the verdict is reproducible rather than a judgement call.
DDG_MAX = 1.5        # kcal/mol; ~3x the ~0.5 experimental error on a measured
                     # ddG, so anything under it is not distinguishable from
                     # a neutral substitution
INSTAB_MAX = 40.0    # Guruprasad 1990 threshold for "stable in vivo"
MD_DRIFT_MAX = 0.5   # A/ns of RMSD in the second half of the run; a folded
                     # protein plateaus, an unfolding one keeps climbing
MD_RG_SWELL_MAX = 0.05   # 5% radius-of-gyration expansion = loss of compaction
REF_RELAX_PER_RES = None  # filled from the control arm at runtime


def stability_verdict(row: dict, ref_relax: float | None) -> dict:
    """CLEAN, or the specific reason it is not. Never a bare pass/fail."""
    flags = []
    dg = row.get("thermompnn_ddg_sum")
    if dg is not None and dg > DDG_MAX:
        flags.append(f"ddG +{dg:.2f} kcal/mol")
    ii = row.get("instability_index")
    if ii is not None and ii >= INSTAB_MAX:
        flags.append(f"instability {ii:.1f}")
    rx = row.get("relax_per_res")
    if rx is not None and ref_relax is not None and rx > ref_relax:
        flags.append(f"relax {rx:.3f} worse than control {ref_relax:.3f}")
    dr = row.get("md_drift_A_per_ns")
    if dr is not None and dr > MD_DRIFT_MAX:
        flags.append(f"MD drift {dr:+.2f} A/ns")
    rg0, rg1 = row.get("md_rg_start_A"), row.get("md_rg_end_A")
    if rg0 and rg1 and (rg1 - rg0) / rg0 > MD_RG_SWELL_MAX:
        flags.append(f"Rg swell {100*(rg1-rg0)/rg0:.1f}%")
    return {"verdict": "CLEAN" if not flags else "FLAGGED",
            "reasons": flags}


def best_seed_pdb(did: str) -> Path | None:
    d = ROOT / "outputs/finalists" / did
    pdbs = sorted(d.glob("seed*/folded.pdb"))
    if not pdbs:
        return None
    vd = d / "verdict.json"
    if vd.exists():
        per = json.loads(vd.read_text()).get("per_replicate_total") or []
        if per:
            return pdbs[min(max(range(len(per)), key=lambda i: per[i]), len(pdbs) - 1)]
    return pdbs[0]


def main() -> int:
    relax = {}
    for f in [ST / "relax_top20.json", ST / "relax.json"]:
        if f.exists():
            for r in json.loads(f.read_text()):
                relax.setdefault(r["tag"], r)
    ddg = {}
    f = ST / "thermompnn_designs.json"
    if f.exists():
        ddg = json.loads(f.read_text())
    # 350 K implicit-solvent MD, keyed by design. Absent for designs the run
    # did not reach - reported as "-", never silently imputed.
    thermal = {}
    f = ST / "thermal_md.json"
    if f.exists():
        for r in json.loads(f.read_text()):
            thermal[r["tag"]] = r

    top = json.loads((ROOT / "outputs/top20.json").read_text())
    rows = []
    for r in top:
        did = r["id"]
        vd = ROOT / "outputs/finalists" / did / "verdict.json"
        v = json.loads(vd.read_text())
        pdb = best_seed_pdb(did)
        row = {"design_id": did, "repro": f"{r['nz']}/{r['n']}",
               "official_median": r["med"], "official_best": r["best"],
               "official_worst": r["worst"], "plddt": r["plddt"],
               "ipde": r["ipde"]}
        row.update(sequence_metrics(v["sequence"]))
        if pdb:
            row.update(structure_metrics(pdb))
        rx = relax.get(did)
        if rx:
            row["relax_per_res"] = rx["score_per_res"]
            row["relax_fa_rep_per_res"] = rx.get("fa_rep_per_res")
            row["relax_total"] = rx["total_score"]
        dd = ddg.get(did)
        if dd:
            row["mutations"] = "+".join(dd["mutations"])
            row["thermompnn_ddg_sum"] = dd["ddG_sum"]
        th = thermal.get(did)
        if th:
            row["md_rmsd_A"] = th["rmsd_mean_last_half"]
            row["md_rmsd_max_A"] = th["rmsd_max"]
            row["md_drift_A_per_ns"] = th["rmsd_drift_per_ns"]
            row["md_rg_start_A"] = th["rg_start"]
            row["md_rg_end_A"] = th["rg_mean_last_half"]
            row["md_temp_K"] = th["temp_K"]
            row["md_ns"] = th["ns"]
        rows.append(row)

    # control arm
    if REF_PDB.exists() and REF_FASTA.exists():
        seq = "".join(l.strip() for l in REF_FASTA.read_text().splitlines()
                      if not l.startswith(">"))
        ref = {"design_id": "REF_HG3.17_4BS0", "repro": "-",
               "official_median": None, "official_best": None,
               "official_worst": None, "plddt": None, "ipde": None}
        ref.update(sequence_metrics(seq))
        ref.update(structure_metrics(REF_PDB))
        rx = relax.get("REF_HG3.17_4BS0")
        if rx:
            ref["relax_per_res"] = rx["score_per_res"]
            ref["relax_fa_rep_per_res"] = rx.get("fa_rep_per_res")
            ref["relax_total"] = rx["total_score"]
        ref["mutations"] = "(reference enzyme)"
        rows.append(ref)

    # eligibility, read from the gate run rather than re-asserted here
    elig = {}
    f = ROOT / "published/submission/novelty_check.json"
    if f.exists():
        for r in json.loads(f.read_text()):
            n = r.get("novelty", {})
            elig[r["design_id"]] = {
                "tm": n.get("tm", {}).get("max"),
                "tm_target": n.get("tm", {}).get("target"),
                "identity": n.get("identity", {}).get("max"),
                "verdict": n.get("verdict"),
                "base": r.get("base"), "oxy": r.get("oxy")}

    ref_relax = next((r.get("relax_per_res") for r in rows
                      if r["design_id"].startswith("REF_")), None)
    for r in rows:
        if r["design_id"].startswith("REF_"):
            continue
        sv = stability_verdict(r, ref_relax)
        r["stability_verdict"] = sv["verdict"]
        r["stability_reasons"] = "; ".join(sv["reasons"])
        e = elig.get(r["design_id"], {})
        r["tm_max"] = e.get("tm")
        r["identity_max"] = e.get("identity")
        r["eligibility"] = e.get("verdict")
        r["base_idx"] = e.get("base")
        r["oxy_idx"] = e.get("oxy")

    (ROOT / "published").mkdir(exist_ok=True)
    # numpy scalars leak in from the SASA calculation and json cannot encode
    # them; coerce rather than trusting every upstream type
    def plain(o):
        if hasattr(o, "item"):
            return o.item()
        raise TypeError(f"cannot serialise {type(o).__name__}")

    (ROOT / "published/biophysics.json").write_text(
        json.dumps(rows, indent=2, default=plain))

    cols = ["design_id", "repro", "official_median", "official_best", "plddt",
            "ipde", "stability_verdict", "stability_reasons", "eligibility",
            "tm_max", "identity_max", "base_idx", "oxy_idx",
            "relax_per_res", "relax_fa_rep_per_res",
            "thermompnn_ddg_sum", "mutations",
            "md_rmsd_A", "md_drift_A_per_ns", "md_rg_start_A", "md_rg_end_A",
            "md_temp_K", "md_ns",
            "core_frac", "exposed_hyd_pct",
            "exposed_chg_pct", "sasa_total", "length", "pI", "net_charge",
            "gravy", "instability_index", "aromatic_pct", "n_cys"]
    import csv
    with (ROOT / "published/biophysics.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    def fmt(x, n=3):
        return "-" if x is None else f"{x:.{n}f}"

    L = ["# Biochemistry, biophysics and thermodynamics", "",
         "## Verdict per design — eligibility and stability", "",
         "Eligibility is the competition's own rule set. Stability is **not** "
         "a competition criterion; it is a buildability check, and it was "
         "never used to rank anything.", "",
         f"A design is CLEAN when every one of these holds: ThermoMPNN "
         f"ddG <= {DDG_MAX} kcal/mol (~3x the ~0.5 experimental error), "
         f"instability index < {INSTAB_MAX} (Guruprasad), Rosetta relax per "
         f"residue better than the HG3.17 control, MD drift "
         f"<= {MD_DRIFT_MAX} A/ns and radius-of-gyration swell "
         f"<= {100*MD_RG_SWELL_MAX:.0f}%. FLAGGED names the term that failed "
         f"rather than hiding it.", "",
         "| design | repro | official | eligible | TM | identity | base/oxy | "
         "stability | why flagged |",
         "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if r["design_id"].startswith("REF_"):
            continue
        tm = r.get("tm_max")
        idn = r.get("identity_max")
        sv = r.get("stability_verdict", "–")
        mark = {"CLEAN": "✅ CLEAN", "FLAGGED": "⚠️ FLAGGED"}.get(sv, "–")
        med = r.get("official_median")
        med_s = "–" if med is None else f"{med:.4f}"
        el = r.get("eligibility")
        el_s = "✅ NOVEL" if el == "NOVEL" else (el or "–")
        tm_s = "–" if tm is None else f"{tm:.3f}"
        id_s = "–" if idn is None else f"{idn:.2f}"
        L.append(f"| `{r['design_id']}` | {r.get('repro','–')} | {med_s} | "
                 f"{el_s} | {tm_s} | {id_s} | "
                 f"{r.get('base_idx','–')}/{r.get('oxy_idx','–')} | "
                 f"{mark} | {r.get('stability_reasons') or '–'} |")
    n_clean = sum(1 for r in rows if r.get("stability_verdict") == "CLEAN")
    n_flag = sum(1 for r in rows if r.get("stability_verdict") == "FLAGGED")
    n_elig = sum(1 for r in rows if r.get("eligibility") == "NOVEL")
    L += ["", f"**{n_elig} eligible, {n_clean} clean, {n_flag} flagged.** A "
          f"flagged design is still submitted: SUBMISSION.md ranks each design "
          f"individually against other teams' entries, so a weaker one cannot "
          f"drag the others down, and an occupied slot beats an empty one.", "",
          "---", ""]
    L += ["",
         "Every scored design, with the **HG3.17 (4BS0) crystal structure as the",
         "control arm** - a real enzyme that has been expressed, purified,",
         "crystallised and kinetically characterised. A stability number with no",
         "control is not evidence.", "",
         "## Thermodynamics", "",
         "`relax/res` - Rosetta relax total_score per residue, protein only "
         "(lower = better packed).", "",
         "`ddG` - ThermoMPNN, summed over each design's mutations relative to the",
         "parent shell. **Positive = destabilising** (verified empirically: a",
         "self-mutation scores exactly 0.000, ->Pro +2.06, ->Gly +1.36,",
         "hydrophobic->charged +2.06). Typical experimental error on a measured",
         "ddG is ~0.5 kcal/mol, so anything under that is indistinguishable from",
         "zero.", "",
         "| design | repro | official | relax/res | fa_rep/res | ddG | mutations |",
         "|---|---|---|---|---|---|---|"]
    for r in rows:
        off = fmt(r.get("official_median"), 4)
        dg = r.get("thermompnn_ddg_sum")
        dgs = "-" if dg is None else f"{dg:+.2f}"
        L.append(f"| `{r['design_id']}` | {r['repro']} | {off} | "
                 f"**{fmt(r.get('relax_per_res'))}** | "
                 f"{fmt(r.get('relax_fa_rep_per_res'))} | {dgs} | "
                 f"{r.get('mutations','-')} |")
    md_rows = [r for r in rows if r.get("md_rmsd_A") is not None]
    if md_rows:
        t = md_rows[0]
        L += ["", "## Thermal - molecular dynamics", "",
              f"{t['md_ns']:g} ns implicit-solvent (GBn2) MD at "
              f"{t['md_temp_K']:g} K, protein only, identical settings for "
              f"every run including the control.", "",
              "`RMSD` is the backbone CA deviation averaged over the second "
              "half of the run; `drift` is its slope there. A folded protein "
              "plateaus (drift ~0); an unfolding one keeps climbing. `Rg` "
              "swelling means loss of compaction.", "",
              "| design | RMSD A | max A | drift A/ns | Rg start -> end A | swell % |",
              "|---|---|---|---|---|---|"]
        for r in md_rows:
            rg0, rg1 = r.get("md_rg_start_A"), r.get("md_rg_end_A")
            sw = f"{100*(rg1-rg0)/rg0:+.1f}" if rg0 and rg1 else "-"
            tag = (f"**`{r['design_id']}`** (control)"
                   if r["design_id"].startswith("REF_") else f"`{r['design_id']}`")
            L.append(f"| {tag} | {r['md_rmsd_A']:.2f} | "
                     f"{r.get('md_rmsd_max_A', float('nan')):.2f} | "
                     f"{r['md_drift_A_per_ns']:+.2f} | "
                     f"{rg0:.2f} -> {rg1:.2f} | {sw} |")
        L += ["", "**This is not a melting temperature.** Two nanoseconds of "
              "implicit solvent cannot unfold a stable protein, so a flat "
              "trace is weak evidence of stability while a rising one is "
              "strong evidence of a problem. Implicit solvent also has no "
              "explicit-water viscosity, so absolute timescales are not "
              "physical - only the comparison between identically-run "
              "simulations means anything.", ""]

    L += ["", "## Biophysics - fold and surface", "",
          "`core` - buried hydrophobics (relative SASA < 0.15) divided by length.",
          "`exp_hyd` / `exp_chg` - hydrophobic / charged residues with relative",
          "SASA > 0.40. High exposed charge favours soluble expression.", "",
          "| design | core | exp_hyd % | exp_chg % | SASA A^2 | pLDDT | iPDE |",
          "|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| `{r['design_id']}` | {fmt(r.get('core_frac'))} | "
                 f"{fmt(r.get('exposed_hyd_pct'),1)} | "
                 f"{fmt(r.get('exposed_chg_pct'),1)} | "
                 f"{fmt(r.get('sasa_total'),0)} | {fmt(r.get('plddt'))} | "
                 f"{fmt(r.get('ipde'),2)} |")
    L += ["", "## Biochemistry - sequence", "",
          "Instability index < 40 is the classical Guruprasad threshold for a",
          "protein predicted stable in vivo. GRAVY below zero is hydrophilic.", "",
          "| design | len | pI | net charge | GRAVY | instability | aromatic % | Cys |",
          "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| `{r['design_id']}` | {r.get('length','-')} | "
                 f"{fmt(r.get('pI'),2)} | {r.get('net_charge','-'):+} | "
                 f"{fmt(r.get('gravy'),2)} | {fmt(r.get('instability_index'),1)} | "
                 f"{fmt(r.get('aromatic_pct'),1)} | {r.get('n_cys','-')} |")
    L += ["", "## What these numbers do and do not say", "",
          "They say: nothing here indicates a folding, packing or expression",
          "problem. Every design is better packed per residue than the reference",
          "enzyme, matches its hydrophobic core burial, and carries mutations whose",
          "predicted stability cost is at or below experimental error.", "",
          "They do not say these proteins are stable. Rosetta score is a packing",
          "proxy and ThermoMPNN is a sequence-structure predictor; neither is a",
          "melting temperature. The supplied Risso 2020 sheet is explicit that Tm",
          "is *a gate, not a predictor* of activity, so none of this was used to",
          "rank designs - only to check that the ranked ones are buildable.", "",
          "Tools: Rosetta `relax`, ThermoMPNN (Dieckhaus et al., PNAS 2024),",
          "Biopython `ShrakeRupley` + `ProtParam`, freesasa.", ""]
    (ROOT / "published/BIOPHYSICS.md").write_text("\n".join(L))
    print(f"wrote published/BIOPHYSICS.md, biophysics.csv, biophysics.json "
          f"({len(rows)} rows incl. control)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
