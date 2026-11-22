"""Stage 3 AIzymes Directed Evolution Refinement Package.

Bridges top-performing de novo scaffolds from Stage 1/Stage 2 into the native AIzymes
evolutionary optimization loop:
- Freezes catalytic residues to prevent destructive mutations.
- Automatically detects and targets second-shell residues (4-8 Angstroms from active site).
- Generates Rosetta `enzdes.cst` catalytic constraints and FieldTools `field_target.dat`.
- Orchestrates multi-round Boltzmann selection for epistatic optimization and affinity maturation.
"""

from hackathon_extensions.stage3_aizymes.aizymes_bridge import AIzymesBridge
from hackathon_extensions.stage3_aizymes.runner import Stage3AIzymesRunner

__all__ = [
    "AIzymesBridge",
    "Stage3AIzymesRunner",
]
