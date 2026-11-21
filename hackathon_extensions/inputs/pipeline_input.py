"""Pipeline Input Data Structures for Enzyme Design Ingestion.

Defines unified representations for all hackathon input scenarios:
- Case 1: Annotated mmCIF / PDB with transition state + fixed catalytic motif.
- Case 2a: Unannotated holo-structure or active site crop with partial diffusion masks.
- Case 2b: SMILES + known catalytic mechanism (no starting 3D protein structure).
- Case 3: SMILES only (catalytic hypothesis or small-molecule binder-first).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union


class CaseType(str, Enum):
    """Supported hackathon input scenario types."""
    CASE1_ANNOTATED = "CASE1_ANNOTATED"
    CASE2A_UNANNOTATED_HOLO = "CASE2A_UNANNOTATED_HOLO"
    CASE2B_SMILES_MECHANISM = "CASE2B_SMILES_MECHANISM"
    CASE3_SMILES_ONLY = "CASE3_SMILES_ONLY"


@dataclass
class LigandSpec:
    """Specification of a ligand or transition-state analog."""
    smiles: str = ""
    conformer_pdb: Optional[str] = None
    conformer_path: Optional[str] = None
    resname: str = "LIG"
    chain_id: str = "Z"
    resnum: int = 1
    reactive_atoms: List[str] = field(default_factory=list)
    charge: int = 0
    atom_metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> LigandSpec:
        return cls(**data)


@dataclass
class TheozymeResidueSpec:
    """Specification of a single catalytic amino acid residue in a theozyme."""
    resname: str  # e.g., CYS, HIS, ASP, GLU, SER
    chain_id: str = "A"
    resnum: int = 1
    target_atoms: List[str] = field(default_factory=list)  # e.g. ["SG"] or ["NE2"]
    target_ligand_atom: Optional[str] = None  # e.g. "C1"
    target_distance_range: Tuple[float, float] = (2.7, 3.2)  # Angstroms
    target_angle_range: Optional[Tuple[float, float]] = None  # Degrees (e.g. 100.0, 115.0)
    target_dihedral_range: Optional[Tuple[float, float]] = None  # Degrees
    pdb_coordinates: Optional[str] = None  # PDB block for this residue if modeled
    role: Optional[str] = None  # e.g. "nucleophile", "general_base", "transition_state_stabilizer"

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> TheozymeResidueSpec:
        # Tuple deserialization handling if passed as list from JSON
        if "target_distance_range" in data and isinstance(data["target_distance_range"], list):
            data["target_distance_range"] = tuple(data["target_distance_range"])
        if "target_angle_range" in data and isinstance(data["target_angle_range"], list):
            data["target_angle_range"] = tuple(data["target_angle_range"])
        if "target_dihedral_range" in data and isinstance(data["target_dihedral_range"], list):
            data["target_dihedral_range"] = tuple(data["target_dihedral_range"])
        return cls(**data)


@dataclass
class TheozymeSpec:
    """Specification of the entire theozyme arrangement."""
    residues: List[TheozymeResidueSpec] = field(default_factory=list)
    description: str = ""
    mechanism_type: str = "custom"  # e.g., cysteine_hydrolase, serine_hydrolase, kemp_eliminase
    geometric_restraints: List[str] = field(default_factory=list)  # Rosetta .cst or RFD4 constraints
    theozyme_pdb: Optional[str] = None  # Complete PDB string containing ligand + catalytic residues
    theozyme_path: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        res_dicts = [r.to_dict() for r in self.residues]
        d = asdict(self)
        d["residues"] = res_dicts
        return d

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> TheozymeSpec:
        residues = [TheozymeResidueSpec.from_dict(r) for r in data.get("residues", [])]
        data_copy = dict(data)
        data_copy["residues"] = residues
        return cls(**data_copy)


@dataclass
class PipelineInput:
    """Unified container for all input parameters consumed by generative design stages."""
    case_type: CaseType
    ligand: LigandSpec
    theozyme: Optional[TheozymeSpec] = None
    scaffold_source_pdb: Optional[str] = None
    scaffold_source_path: Optional[str] = None
    partial_diffusion_mask: Optional[Dict[str, Any]] = None
    rfd4_condition_spec: Dict[str, Any] = field(default_factory=dict)
    rfd3_config: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def validate(self) -> Tuple[bool, List[str]]:
        """Validate pipeline input completeness and structural consistency."""
        errors: List[str] = []

        if not isinstance(self.case_type, CaseType):
            try:
                self.case_type = CaseType(self.case_type)
            except ValueError:
                errors.append(f"Invalid case_type: {self.case_type}. Must be one of {[c.value for c in CaseType]}")

        # Check ligand
        if not self.ligand:
            errors.append("LigandSpec is missing.")
        else:
            if not self.ligand.smiles and not self.ligand.conformer_pdb and not self.ligand.conformer_path:
                errors.append("Ligand must have at least one of: smiles, conformer_pdb, or conformer_path.")
            if len(self.ligand.resname) > 3:
                errors.append(f"Ligand residue name '{self.ligand.resname}' exceeds 3 characters.")
            if len(self.ligand.chain_id) != 1:
                errors.append(f"Ligand chain ID '{self.ligand.chain_id}' must be exactly 1 character.")

        # Check scenario requirements
        if self.case_type == CaseType.CASE1_ANNOTATED:
            if not self.theozyme or not self.theozyme.residues:
                errors.append("Case 1 requires an annotated theozyme with at least one catalytic residue.")
        elif self.case_type == CaseType.CASE2A_UNANNOTATED_HOLO:
            if not self.scaffold_source_pdb and not self.scaffold_source_path:
                errors.append("Case 2a requires scaffold_source_pdb or scaffold_source_path.")
        elif self.case_type == CaseType.CASE2B_SMILES_MECHANISM:
            if not self.theozyme or not self.theozyme.residues:
                errors.append("Case 2b requires catalytic residue specifications for theozyme building.")
        elif self.case_type == CaseType.CASE3_SMILES_ONLY:
            if not self.ligand.smiles:
                errors.append("Case 3 requires a valid SMILES string.")

        return (len(errors) == 0, errors)

    def to_dict(self) -> Dict[str, Any]:
        """Serialize PipelineInput to a Python dictionary."""
        return {
            "case_type": self.case_type.value if isinstance(self.case_type, CaseType) else str(self.case_type),
            "ligand": self.ligand.to_dict() if self.ligand else None,
            "theozyme": self.theozyme.to_dict() if self.theozyme else None,
            "scaffold_source_pdb": self.scaffold_source_pdb,
            "scaffold_source_path": self.scaffold_source_path,
            "partial_diffusion_mask": self.partial_diffusion_mask,
            "rfd4_condition_spec": self.rfd4_condition_spec,
            "rfd3_config": self.rfd3_config,
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> PipelineInput:
        """Construct PipelineInput from a dictionary."""
        case_type = CaseType(data["case_type"])
        ligand = LigandSpec.from_dict(data["ligand"]) if data.get("ligand") else LigandSpec()
        theozyme = TheozymeSpec.from_dict(data["theozyme"]) if data.get("theozyme") else None
        return cls(
            case_type=case_type,
            ligand=ligand,
            theozyme=theozyme,
            scaffold_source_pdb=data.get("scaffold_source_pdb"),
            scaffold_source_path=data.get("scaffold_source_path"),
            partial_diffusion_mask=data.get("partial_diffusion_mask"),
            rfd4_condition_spec=data.get("rfd4_condition_spec", {}),
            rfd3_config=data.get("rfd3_config", {}),
            metadata=data.get("metadata", {}),
        )

    def save_summary(self, output_path: Union[str, Path]) -> None:
        """Write a formatted JSON summary to disk."""
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        is_valid, errors = self.validate()
        summary = {
            "valid": is_valid,
            "validation_errors": errors,
            "data": self.to_dict(),
        }
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)
