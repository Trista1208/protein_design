import shutil
from pathlib import Path
import numpy as np
import pandas as pd
from atomworks.io.parser import parse
from rfproteina.utils.standard_annotations import SA_CRDSRC

OUT_DIR = Path("/cluster/scratch/gomest/hackaton/enzyme-design-hackathon-2026/rfd4/enzdes_tutorial/runs/scaffold_n25")
CSV_PATH = OUT_DIR / "inference_metrics" / "per_example_metrics.csv"
BEST_DIR = OUT_DIR / "best"
TOP_K = 5
RMSD_COL = "design_unindexed_motif_allatom_rmsd.rmsd"

df = pd.read_csv(CSV_PATH)
if "ame.pass" in df.columns:
    ranked = df.sort_values(["ame.pass", RMSD_COL], ascending=[False, True]).reset_index(drop=True)
    passed = ranked[ranked["ame.pass"] == 1]
    if passed.empty:
        raise RuntimeError(f"0 of {len(ranked)} designs passed ame.pass; inspect {CSV_PATH}")
else:
    ranked = df.sort_values(RMSD_COL, ascending=True).reset_index(drop=True)
    passed = ranked

BEST_DIR.mkdir(parents=True, exist_ok=True)
ranked.to_csv(BEST_DIR / "selection.csv", index=False)

copied_files = []
for rank_idx, row in passed.head(TOP_K).iterrows():
    src = Path(row["designed_structure_path"])
    dst = BEST_DIR / f"rank{rank_idx + 1}_{src.name}"
    shutil.copy2(src, dst)
    copied_files.append(dst)

best_row = passed.iloc[0]
best_cif = copied_files[0]
print(f"Selected: {best_row['example_id']} ({len(passed)}/{len(ranked)} candidates, "
      f"motif RMSD = {best_row[RMSD_COL]:.3f} A)")
print(f"Saved top {len(copied_files)} designs + selection.csv to: {BEST_DIR}\n")

aa = parse(str(best_cif))["assemblies"]["1"][0]
motif = aa.mask("sa_mask('SA_CRDSRC') & ~hetero")
src_labels = np.asarray(SA_CRDSRC.annotation(aa))
resids = sorted(set(aa.res_id[motif].astype(int).tolist()))

print(f"Catalytic motif mapping in {best_cif.name}: {int(motif.sum())} atoms across {len(resids)} residues")
for rid in resids:
    m = motif & (aa.res_id == rid)
    print(f"  scaffold res {rid:>4} {aa.res_name[m][0]:<4} <- input {src_labels[m][0]}"
          f"  atoms={sorted(aa.atom_name[m].tolist())}")
print(f"  ligand (DUT) atoms: {int(aa.mask('hetero').sum())}\n")

print("Run on your LOCAL PowerShell to inspect in local PyMOL/ChimeraX:")
print(f"scp -r gomest@euler.ethz.ch:{BEST_DIR} .")
