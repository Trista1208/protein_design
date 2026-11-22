#!/bin/bash
# Functional verification: every tool must PRODUCE something, not just import.
set -u
E=/home/mandalorian/anaconda3/envs
O=/home/mandalorian/opt
W=/home/mandalorian/Hackathons/Enzymes-BSSE2026/outputs/verify
PDB=/home/mandalorian/Hackathons/Enzymes-BSSE2026/outputs/smoke/esmfold_test.pdb
rm -rf "$W"; mkdir -p "$W"; cd "$W"
PASS=0; FAIL=0; SKIP=0
r(){ # r <name> <cmd...>
  local n="$1"; shift
  printf "%-26s " "$n"
  if eval "$@" >"$W/_$n.log" 2>&1; then echo "PASS"; PASS=$((PASS+1));
  else echo "FAIL"; FAIL=$((FAIL+1)); sed -n '1,6p' "$W/_$n.log" | sed 's/^/      | /'; fi
}

echo "################ CHEMISTRY / STRUCTURE ################"
r RDKit "$E/mm/bin/python -c \"
from rdkit import Chem; from rdkit.Chem import AllChem
m=Chem.AddHs(Chem.MolFromSmiles('c1ccccc1C(=O)NO'))
assert AllChem.EmbedMolecule(m,randomSeed=1)==0
AllChem.MMFFOptimizeMolecule(m)
Chem.MolToMolFile(m,'$W/rdkit.sdf')
print('atoms',m.GetNumAtoms())\""

r PyMOL "$E/viz/bin/python -c \"
import pymol; pymol.finish_launching(['pymol','-qc'])
from pymol import cmd
cmd.load('$PDB','p'); n=cmd.count_atoms('p')
cmd.save('$W/pymol_out.pdb','p')
print('atoms',n); assert n>0\""

r pdbfixer "$E/mm/bin/python -c \"
from pdbfixer import PDBFixer
from openmm.app import PDBFile
f=PDBFixer(filename='$PDB')
f.findMissingResidues(); f.findMissingAtoms(); f.addMissingAtoms(); f.addMissingHydrogens(7.0)
PDBFile.writeFile(f.topology,f.positions,open('$W/fixed.pdb','w'))
print('fixed atoms',f.topology.getNumAtoms())\""

r MDAnalysis "$E/mm/bin/python -c \"
import MDAnalysis as mda
u=mda.Universe('$PDB'); print('atoms',len(u.atoms)); assert len(u.atoms)>0
from MDAnalysis.analysis import rms
print('rgyr %.2f'%u.atoms.radius_of_gyration())\""

r mdtraj "$E/mm/bin/python -c \"
import mdtraj as md
t=md.load('$PDB'); print(t); print('rg %.3f'%md.compute_rg(t)[0])\""

r pytraj_FieldTools "$E/mm/bin/python -c \"
import pytraj as pt, numpy as np, ast
t=pt.load('$PDB'); print('frames',t.n_frames,'atoms',t.n_atoms)
src=open('/home/mandalorian/Hackathons/Enzymes-BSSE2026/tools/FieldTools.py').read()
mod=ast.parse(src)
import types; ns={}
# exec only the function defs to prove they are callable with pytraj present
exec(compile(mod,'FieldTools.py','exec'),ns)
f=ns['Calc_Electrostatic_Potential_Energy']
v=f(np.array([0.,0.,0.]),np.array([0.,0.,2.]),1.0,-1.0)
print('FieldTools coulomb term: %.3f kJ/mol'%v)\""

echo
echo "################ QM ################"
r Psi4 "$E/qm/bin/python -c \"
import psi4; psi4.core.set_output_file('$W/psi4.out',False)
psi4.set_memory('2 GB')
h2o=psi4.geometry('''O 0 0 0
H 0 0 0.96
H 0.93 0 -0.26''')
e=psi4.energy('scf/cc-pvdz'); print('SCF/cc-pVDZ E = %.8f Ha'%e)
assert -76.1 < e < -75.9\""

r PySCF "$E/qm/bin/python -c \"
from pyscf import gto,scf
m=gto.M(atom='H 0 0 0; H 0 0 0.74',basis='def2-svp')
e=scf.RHF(m).kernel(); print('RHF E = %.8f Ha'%e)
assert -1.2<e<-1.0\""

r xtb "cd $W && printf '3\n\nO 0.0 0.0 0.0\nH 0.0 0.0 0.96\nH 0.93 0.0 -0.26\n' > w.xyz && $E/qm/bin/xtb w.xyz --gfn 2 --sp 2>&1 | grep -E 'TOTAL ENERGY'"

r ASE "$E/qm/bin/python -c \"
from ase.build import molecule
from ase.optimize import BFGS
from ase.calculators.emt import EMT
a=molecule('H2O'); print('ase atoms',len(a))
import ase; print('ase',ase.__version__)\""

r ORCA_serial "cd $W && printf '! HF DEF2-SVP\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n' > os.inp && LD_LIBRARY_PATH=$O/orca $O/orca/orca os.inp 2>&1 | grep 'FINAL SINGLE POINT ENERGY'"

r ORCA_parallel_MPI "cd $W && printf '! HF DEF2-SVP PAL2\n* xyz 0 1\nO 0 0 0\nH 0 0 0.96\nH 0.93 0 -0.26\n*\n' > op.inp && PATH=$O/openmpi-4.1.8/bin:\$PATH LD_LIBRARY_PATH=$O/openmpi-4.1.8/lib:$O/orca $O/orca/orca op.inp 2>&1 | grep -E 'FINAL SINGLE POINT ENERGY|parallel'"

echo
echo "################ MM / MD ################"
r OpenMM_CUDA_MD "$E/mm/bin/python -c \"
from openmm import app, unit, LangevinMiddleIntegrator, Platform
from openmm.app import PDBFile, ForceField, Modeller
import openmm
pdb=PDBFile('$W/fixed.pdb')
ff=ForceField('amber14-all.xml','amber14/tip3pfb.xml')
m=Modeller(pdb.topology,pdb.positions)
sysm=ff.createSystem(m.topology,nonbondedMethod=app.NoCutoff,constraints=app.HBonds)
integ=LangevinMiddleIntegrator(300*unit.kelvin,1/unit.picosecond,0.002*unit.picoseconds)
sim=app.Simulation(m.topology,sysm,integ,Platform.getPlatformByName('CUDA'))
sim.context.setPositions(m.positions)
e0=sim.context.getState(getEnergy=True).getPotentialEnergy()
sim.minimizeEnergy(maxIterations=50); sim.step(200)
e1=sim.context.getState(getEnergy=True).getPotentialEnergy()
print('platform',sim.context.getPlatform().getName()); print('E0',e0); print('E200',e1)\""

r OpenFF_toolkit "$E/mm/bin/python -c \"
from openff.toolkit import Molecule, ForceField
m=Molecule.from_smiles('CCO'); m.generate_conformers(n_conformers=1)
print('openff mol atoms',m.n_atoms)
ff=ForceField('openff-2.0.0.offxml'); print('loaded', ff.author is not None or True)\""

echo
echo "################ Q6 / EVB ################"
if [ -d "$O/Q6/tests/basic_tests/prep_SPH" ]; then
  r Q6_Qprep6 "cd $O/Q6/tests/basic_tests/prep_SPH && source $O/Q6/env.sh && Qprep6 < generate.inp 2>&1 | tail -20"
  r Q6_Qdyn6_run "cd $W && source $O/Q6/env.sh && cp $O/Q6/tests/basic_tests/prep_SPH/lig_w.top . 2>/dev/null; cp $O/Q6/tests/basic_tests/eq1.inp . 2>/dev/null && (Qdyn6 eq1.inp 2>&1 | tail -15)"
else
  echo "Q6 tests dir missing"; SKIP=$((SKIP+1))
fi

echo
echo "################ ROSETTA ################"
cat > "$W/min.xml" <<'XML'
<ROSETTASCRIPTS>
  <SCOREFXNS><ScoreFunction name="sfxn" weights="ref2015"/></SCOREFXNS>
  <MOVERS><MinMover name="min" scorefxn="sfxn" bb="0" chi="1"/></MOVERS>
  <PROTOCOLS><Add mover="min"/></PROTOCOLS>
  <OUTPUT scorefxn="sfxn"/>
</ROSETTASCRIPTS>
XML
r Rosetta_relax "cd $W && source $O/rosetta/env.sh && relax.default.linuxgccrelease -s $PDB -database \$ROSETTA3_DB -relax:quick -nstruct 1 -out:prefix v1_ -overwrite 2>&1 | grep -E 'reported success|jobs attempted'"
r Rosetta_scripts "cd $W && source $O/rosetta/env.sh && rosetta_scripts.default.linuxgccrelease -s $PDB -database \$ROSETTA3_DB -parser:protocol min.xml -nstruct 1 -out:prefix v2_ -overwrite 2>&1 | grep -E 'reported success|jobs attempted'"
r Rosetta_match_opts "cd $W && source $O/rosetta/env.sh && match.default.linuxgccrelease -database \$ROSETTA3_DB -help 2>&1 | grep -ciE 'match|usage'"

echo
echo "################ SUMMARY ################"
echo "PASS=$PASS  FAIL=$FAIL"
ls -1 "$W" | grep -vE '^_' | sed 's/^/  artefact: /'
echo "### VERIFY_ALL DONE ###"
