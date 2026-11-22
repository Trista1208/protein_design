#!/usr/bin/env python
"""Full-score control test: HG3.17 (Gln50) vs HG3.17-Q50M (Met50).

The earlier controls run compared `complex_ipde` only and FAILED: the knockout
scored better (2.64 vs 3.09) and the separation (0.45) was inside the noise
floor (0.70). That told us iPDE is a weak discriminator - consistent with the
organisers demoting it 0.30 -> 0.20 for correlating at rho = 0.93 with
ligand_iptm - but it said nothing about the 65% of the score that is geometry.

This version tests the SCORER on the CRYSTAL STRUCTURES, which removes every
confound at once: no folding, no sequence-to-structure numbering map, no Boltz
variance. If the scorer cannot rank a real Gln50 structure above its own Q50M
knockout, no amount of downstream work will help.

  4BS0 chain A : HG3.17 E47N/N300D + 6NT, 1.09 A, REAL refined Gln50
  7K4Z chain A : HG3.17 + 6NT, 1.08 A, REAL refined Gln50   (second positive)
  6TU6 chain G : HG3.17 Q50M + 6NT, 1.05 A, MET50 - the oxyanion knockout
                 (also the source of the supplied theozyme)

A separate pipeline-level control (fold the sequence, then score) remains
worth running, but it cannot be interpreted until this one passes.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import apnovo, config as C, official      # noqa: E402

CASES = [
    ("positive_4BS0", "4BS0", "A", 127, 50, "HG3.17, real refined Gln50"),
    ("positive_7K4Z", "7K4Z", "A", 127, 50, "HG3.17, real refined Gln50"),
    ("negative_6TU6", "6TU6", "G", 127, 50, "HG3.17-Q50M, MET50 knockout"),
]
CORPUS = C.ROOT / "inputs" / "corpus" / "structures"
OUT = C.RUNS / "controls_full"


def extract_chain(cif: Path, chain: str, dst: Path) -> Path:
    """Pull one chain plus its 6NT into a PDB, stamping B>0 on the motif only.

    BRIEF 3: APNovo reads atom_b_factor > 0 as "this residue is motif", so a
    blanket B=1.00 across a whole protein makes EVERY residue motif and
    corrupts both the H-bond metrics and the burial anchor. We stamp B only on
    the two catalytic residues and the ligand, and pass motif_residues anyway.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    code = f'''
import gemmi
st = gemmi.read_structure({str(cif)!r})
st.setup_entities()
st.remove_alternative_conformations()
st.remove_hydrogens()
sel = gemmi.Structure()
sel.cell, sel.spacegroup_hm = st.cell, st.spacegroup_hm
model = gemmi.Model("1")
for ch in st[0]:
    if ch.name != {chain!r}:
        continue
    nc = gemmi.Chain(ch.name)
    for res in ch:
        # keep polymer residues and the TARGET LIGAND only. Ions, cryoprotectant
        # and waters add non-protein residues that confuse the motif/ligand
        # bookkeeping and contribute nothing to the metrics.
        info = gemmi.find_tabulated_residue(res.name)
        is_poly = info is not None and info.is_amino_acid()
        if not (is_poly or res.name == "6NT"):
            continue
        nc.add_residue(res)
    model.add_chain(nc)
sel.add_model(model)
sel.setup_entities()
sel.write_pdb({str(dst)!r})
print("ok")
'''
    proc = subprocess.run(
        ["/home/mandalorian/anaconda3/envs/rfd4/bin/python", "-c", code],
        capture_output=True, text=True, timeout=600)
    if proc.returncode != 0 or not dst.exists():
        raise RuntimeError(f"chain extraction failed for {cif.name}:{chain}: "
                           f"{proc.stderr[-300:]}")
    return dst


def stamp_motif_bfactors(pdb: Path, keep: set[int], ligand: str) -> Path:
    """B = 1.00 on the motif residues and ligand ONLY, 0.00 everywhere else."""
    out = pdb.with_name(pdb.stem + "_bstamp.pdb")
    lines = []
    for line in pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                resseq = int(line[22:26])
            except ValueError:
                lines.append(line)
                continue
            resname = line[17:20].strip()
            b = 1.00 if (resseq in keep or resname == ligand) else 0.00
            line = f"{line[:60]}{b:6.2f}{line[66:]}"
        lines.append(line)
    out.write_text("\n".join(lines) + "\n")
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict] = {}

    for label, pdbid, chain, base, oxy, note in CASES:
        cif = CORPUS / f"{pdbid}.cif"
        if not cif.exists():
            print(f"  {label}: MISSING {cif}")
            continue
        donor_absent = False
        try:
            raw_pdb = extract_chain(cif, chain, OUT / f"{pdbid}_{chain}.pdb")
            stamped = stamp_motif_bfactors(raw_pdb, {base, oxy}, C.LIGAND_CCD)
            am = apnovo.metrics(stamped, motif_residues=[base, oxy])
        except Exception as exc:                       # noqa: BLE001
            # The Q50M knockout has MET at 50, and Met is NOT in APNovo's
            # KE_MOTIF_ROLES donor set (N/Q/S/T/Y/H/K/R). So "no valid catalytic
            # role combinations" is not an error - it is the measurement: this
            # structure has NO assignable oxyanion donor. Under the official
            # geometric mean that means s_oxy = 0, hence S = 0.
            if "role combinations" not in str(exc):
                print(f"  {label}: FAILED {type(exc).__name__}: {str(exc)[:160]}")
                continue
            # APNovo needs >=2 assignable catalytic residues, so a knockout
            # with no donor cannot be scored by it at all. Fall back to our own
            # measurement for the terms that DO exist (base, burial) and record
            # the oxyanion terms as hard zeros.
            donor_absent = True
            from kemp.theozyme import load as _load
            tzk = _load(stamped, require_ligand=False)
            site = official.CatalyticSite(base_resseq=base, oxy_resseq=oxy,
                                          chain=chain)
            r = official.measure_raw(tzk, site=site)
            am = {"_donor_absent": True,
                  "_residue_at_oxy": tzk.resname_at(oxy, chain=chain)}
        if not donor_absent:
            r = apnovo.to_raw(am)
        if donor_absent:
            # absent donor => the oxyanion terms score ZERO, not "missing"
            r.d_oxy = None
            r.angle2_oxy = None
            r.dihedral12_oxy = None
            r.notes.append("NO oxyanion donor: residue 50 is not in "
                           "KE_MOTIF_ROLES (Met in the Q50M knockout). oxy and "
                           "orient_oxy are ZERO, not unmeasurable.")
        # single static structure: no replicates, so no cons and no ipde
        sc = official.score(r, strict=False)
        if donor_absent:
            sc["subscores"]["oxy"] = 0.0
            sc["subscores"]["orient_oxy"] = 0.0
            sc["prefilter_score"] = 0.0
            sc["zeroed_terms"] = ["oxy", "orient_oxy"]
            sc["limiting_term"] = "oxy"
        results[label] = {
            "pdb": pdbid, "chain": chain, "note": note,
            "d_base": r.d_base, "d_oxy": r.d_oxy,
            "angle2_base": r.angle2_base, "dihedral2_base": r.dihedral12_base,
            "angle2_oxy": r.angle2_oxy, "dihedral2_oxy": r.dihedral12_oxy,
            "ca_neighbors_mean": r.ca_neighbors_mean,
            "subscores": sc["subscores"],
            "prefilter_score": sc.get("prefilter_score"),
            "limiting_term": sc.get("limiting_term"),
            "hbond_base": am.get("acid_base/hbond"),
            "hbond_oxy": am.get("oxyanion_oah/hbond"),
            "donor_absent": donor_absent,
            "residue_at_oxy": am.get("_residue_at_oxy"),
        }
        print(f"\n  === {label} ({pdbid} chain {chain}) - {note}")
        print(f"      d_base={r.d_base:.3f}  d_oxy="
              f"{'n/a' if r.d_oxy is None else format(r.d_oxy, '.3f')}  "
              f"bury={r.ca_neighbors_mean:.2f}")
        if donor_absent:
            print(f"      NO OXYANION DONOR: residue {oxy} is "
                  f"{am.get('_residue_at_oxy')} - not in KE_MOTIF_ROLES")
        else:
            print(f"      hbond base={am.get('acid_base/hbond')} "
                  f"oxy={am.get('oxyanion_oah/hbond')}")
        print(f"      subscores={sc['subscores']}")
        print(f"      score (geom. mean over measurable terms)="
              f"{sc.get('prefilter_score')}  limiting={sc.get('limiting_term')}")

    # ---- verdict
    pos = [v["prefilter_score"] for k, v in results.items()
           if k.startswith("positive") and v.get("prefilter_score") is not None]
    neg = [v["prefilter_score"] for k, v in results.items()
           if k.startswith("negative") and v.get("prefilter_score") is not None]
    print("\n" + "=" * 68)
    if not pos or not neg:
        print("  INCONCLUSIVE - missing a positive or negative case")
        verdict = {"ranks_correctly": None, "reason": "incomplete cases"}
    else:
        ok = min(pos) > max(neg)
        print(f"  positives: {[round(p,4) for p in pos]}")
        print(f"  negative : {[round(n,4) for n in neg]}")
        print(f"  separation (worst positive - negative) = "
              f"{min(pos) - max(neg):+.4f}")
        print(f"  RANKS CORRECTLY: {ok}")
        if ok:
            print("  => the scorer DOES detect oxyanion stabilisation on real")
            print("     structures. Geometry is doing the work that iPDE could not.")
        else:
            print("  => the scorer CANNOT separate a real Gln50 from its Q50M")
            print("     knockout. Fix the scoring windows before generating designs.")
        verdict = {"ranks_correctly": bool(ok),
                   "positives": pos, "negative": neg,
                   "separation": min(pos) - max(neg)}

    (OUT / "verdict.json").write_text(json.dumps(
        {"cases": results, "verdict": verdict}, indent=2, default=str))
    print(f"\n  wrote {OUT / 'verdict.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
