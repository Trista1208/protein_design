"""Measure chain integrity, catalytic strain and catalytic distances on every
RFD3 output we have. Prints the distributions the STRAIN_GATE thresholds are
calibrated from. Read-only."""
import glob, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import strain

def site_for(run_dir: Path):
    for js in sorted(run_dir.rglob("*_model_*.json")):
        try:
            d = json.loads(js.read_text())
        except json.JSONDecodeError:
            continue
        dim = d.get("diffused_index_map")
        if isinstance(dim, dict):
            return dim
    return None

rows = []
for cif in sorted(glob.glob("outputs/**/rfd3/*model_0.cif*", recursive=True)):
    cif = Path(cif)
    dim = site_for(cif.parent)
    if not dim:
        print(f"{cif}: no diffused_index_map"); continue
    # map tags -> designed positions, e.g. {"A127": "A131", "A50": "A100"}
    vals = {k: v for k, v in dim.items()}
    nums = sorted(int(v[1:]) for v in vals.values())
    atoms = strain.load_atoms(cif)
    res = strain.residues(atoms)
    base = oxy = None
    for v in vals.values():
        ch, seq = v[0], int(v[1:])
        rn = strain._resname(atoms, (ch, seq))
        if rn == "ASP": base = seq
        elif rn == "GLN": oxy = seq
    if base is None or oxy is None:
        print(f"{cif}: motif not ASP/GLN at {vals} -> {[ (v, strain._resname(atoms,(v[0],int(v[1:])))) for v in vals.values()]}")
        continue
    v = strain.validate(cif, base_resseq=base, oxy_resseq=oxy)
    rows.append((str(cif), base, oxy, v))
    print(f"\n=== {cif}")
    print(f"  motif: base={base} oxy={oxy}  map={vals}")
    i, s, d = v.integrity, v.strain, v.distances
    print(f"  CHAIN    bonds={i['n_peptide_bonds']} cuts={i['n_cuts']} nonphysical={i['n_nonphysical']} "
          f"distorted={i['n_distorted']}")
    print(f"           C-N mean={i['cn_mean_A']} min={i['cn_min_A']} max={i['cn_max_A']}")
    print(f"           max |dev|: at motif junction={i['junction_max_dev_A']}  everywhere else={i['nonjunction_max_dev_A']}")
    for b in i['junction_bonds']:
        print(f"           junction {b['between'][0]}-{b['between'][1]:6} C-N={b['C_N_A']:>6}  dev {b['dev']:+.3f}"
              f"{'  <-- CUT' if b['cut'] else ('  <-- NON-PHYSICAL' if b['nonphysical'] else '')}")
    for b in i['cuts']:
        if not b['junction']: print(f"           CUT (not at motif) {b['between']} C-N={b['C_N_A']}")
    for role, rv in s["residues"].items():
        if not rv.get("present"):
            print(f"  {role.upper():8} MISSING"); continue
        print(f"  {role.upper():8} {rv['resname']}{rv['resseq']} atoms={rv['n_atoms']} "
              f"strained={rv['n_strained']} omega_planar={rv['omega_planar']} dih={rv['dihedrals']}")
        wb, wa = rv["worst_bond"], rv["worst_angle"]
        if wb: print(f"           worst bond : {wb['atoms']:12} {wb['value']:>7} vs {wb['ideal']:>6}  dev {wb['dev']:+.3f}")
        if wa: print(f"           worst angle: {wa['atoms']:12} {wa['value']:>7} vs {wa['ideal']:>6}  dev {wa['dev']:+.2f}")
        for t in rv["terms"]:
            if t["strained"]:
                print(f"           STRAINED {t['kind']:5} {t['atoms']:12} {t['value']:>7} vs {t['ideal']:>6} dev {t['dev']:+.3f}")
    print(f"  LIGAND   present={d['ligand_present']} atoms={d['ligand_atoms']}")
    for k in ("d_base", "d_oxy", "ca_ca", "od2_n3_n1_angle_deg"):
        tt = d[k]
        print(f"  {k:20} {str(tt['value']):>9}  ref {tt['reference']:>7}  dev {str(tt['deviation']):>8}  "
              f"{'ok' if tt['within_tol'] else ('UNMEASURABLE' if not tt['measurable'] else 'OFF')}")
    print(f"  VERDICT  {v.summary()}")

print("\n\n================ SUMMARY ================")
for path, base, oxy, v in rows:
    print(f"{'PASS' if v.passed else 'FAIL':5} {Path(path).parts[1]:28} {v.summary()}")
