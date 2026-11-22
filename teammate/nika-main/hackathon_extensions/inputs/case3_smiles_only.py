"""Case 3 Parser: SMILES Only (Zero 3D Starting Protein Structure).

Handles raw SMILES input molecules under two primary pathways:
1. Theozyme-First Mode (Step 2A):
   User supplies a hypothesized catalytic mechanism or triad type (e.g. "cysteine_hydrolase").
   The parser generates a 3D conformer, attaches standardized catalytic theozyme geometry,
   and emits catalytic scaffolding specifications for RFD4 and RFD3.

2. Binder-First Shortcut (Step 2B):
   When zero catalytic hypothesis exists, configures RFD4 and RFD3 in tight
   'small_molecule_binder' mode with all ligand atoms buried and fixed.
   Prepares handover metadata and configurations for Stage 3 (AIzymes) to evolve
   catalytic residues into the pocket shell post-generation.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from hackathon_extensions.inputs.case1_annotated import Case1AnnotatedParser
from hackathon_extensions.inputs.conformer_generator import ConformerGenerator
from hackathon_extensions.inputs.pipeline_input import (
    CaseType,
    LigandSpec,
    PipelineInput,
    TheozymeSpec,
)
from hackathon_extensions.inputs.theozyme_builder import TheozymeBuilder


class Case3SmilesOnlyParser:
    """Ingestion parser for SMILES-only inputs."""

    def __init__(
        self,
        conformer_generator: Optional[ConformerGenerator] = None,
        theozyme_builder: Optional[TheozymeBuilder] = None,
        default_scaffold_length_range: Tuple[int, int] = (150, 300),
    ):
        self.conformer_generator = conformer_generator or ConformerGenerator()
        self.theozyme_builder = theozyme_builder or TheozymeBuilder()
        self.default_scaffold_length_range = default_scaffold_length_range

    def parse(
        self,
        smiles: str,
        mode: str = "auto",  # "auto", "theozyme_first", "binder_first"
        hypothesized_mechanism: Optional[str] = None,  # e.g. "cysteine_hydrolase", "serine_hydrolase"
        reactive_atoms: Optional[List[str]] = None,
        scaffold_length_range: Optional[Tuple[int, int]] = None,
        output_dir: Optional[Union[str, Path]] = None,
    ) -> PipelineInput:
        """Parse raw SMILES into a generative design pipeline input.

        Parameters
        ----------
        smiles : str
            Input chemical SMILES string.
        mode : str
            "auto" chooses theozyme_first if hypothesized_mechanism is provided, else binder_first.
            "theozyme_first" requires or infers a catalytic triad.
            "binder_first" generates a tight pocket binder first, deferring catalytic triad introduction to Stage 3 AIzymes.
        hypothesized_mechanism : Optional[str]
            Reaction mechanism name ("cysteine_hydrolase", "serine_hydrolase", "kemp_eliminase").
        reactive_atoms : Optional[List[str]]
            Ligand reactive atoms (default inferred or "C1").
        scaffold_length_range : Optional[Tuple[int, int]]
            Length range of de novo scaffold (default 150-300).
        output_dir : Optional[Union[str, Path]]
            Destination folder for output files.

        Returns
        -------
        PipelineInput
            Unified pipeline input specification.
        """
        smiles = smiles.strip()
        if not smiles:
            raise ValueError("SMILES string cannot be empty.")

        scaffold_range = scaffold_length_range or self.default_scaffold_length_range
        out_p = Path(output_dir) if output_dir else None

        # Step 1: Generate 3D conformer with unique atom names
        lig_pdb_path = (out_p / "ligand_conformer.pdb") if out_p else None
        conformer_pdb, meta = self.conformer_generator.generate_conformer(
            smiles=smiles,
            output_pdb=lig_pdb_path,
            resname="LIG",
            chain_id="Z",
            resnum=1,
        )

        chosen_reactive = reactive_atoms or meta.get("reactive_candidates") or ["C1"]
        ligand_spec = LigandSpec(
            smiles=smiles,
            conformer_pdb=conformer_pdb,
            conformer_path=str(lig_pdb_path) if lig_pdb_path else None,
            resname="LIG",
            chain_id="Z",
            resnum=1,
            reactive_atoms=chosen_reactive,
            atom_metadata=meta,
        )

        # Determine strategy
        if mode == "auto":
            is_theozyme = hypothesized_mechanism is not None
        elif mode == "theozyme_first":
            is_theozyme = True
        elif mode == "binder_first":
            is_theozyme = False
        else:
            raise ValueError(f"Unknown mode: {mode}. Choose from 'auto', 'theozyme_first', 'binder_first'.")

        if is_theozyme:
            # Step 2A: Theozyme-First Mode
            mechanism = hypothesized_mechanism or "cysteine_hydrolase"
            theozyme_pdb_path = (out_p / "theozyme_input.pdb") if out_p else None
            theozyme_cif_path = (out_p / "theozyme_input.cif") if out_p else None
            theozyme_cst_path = (out_p / "enzdes.cst") if out_p else None

            theozyme_spec = self.theozyme_builder.build_theozyme(
                ligand=ligand_spec,
                mechanism_type=mechanism,
                output_pdb=theozyme_pdb_path,
                output_cif=theozyme_cif_path,
                output_cst=theozyme_cst_path,
            )

            annot_parser = Case1AnnotatedParser()
            rfd4_yaml = annot_parser._build_rfd4_condition_spec(
                theozyme_spec.residues,
                ligand_spec,
                scaffold_range,
            )
            rfd3_json = annot_parser._build_rfd3_config(
                theozyme_spec.residues,
                ligand_spec,
                scaffold_range,
            )

            metadata = {
                "strategy": "theozyme_first",
                "hypothesized_mechanism": mechanism,
                "smiles": smiles,
                "notes": "Scaffolding around idealized catalytic triad + ligand",
            }

        else:
            # Step 2B: Binder-First Shortcut
            theozyme_spec = None

            rfd4_yaml = {
                "version": "rfd4-proteina-v1",
                "target": {
                    "design_mode": "small_molecule_binder",
                    "scaffold_length_range": list(scaffold_range),
                    "default_length": int((scaffold_range[0] + scaffold_range[1]) / 2),
                },
                "ligand": {
                    "chain": ligand_spec.chain_id,
                    "resname": ligand_spec.resname,
                    "resnum": ligand_spec.resnum,
                    "fixed_atoms": ["ALL"],
                    "select_buried": ["ALL"],
                    "reactive_atoms": ligand_spec.reactive_atoms,
                },
                "fixed_motifs": [],
                "diffusion_params": {
                    "num_timesteps": 50,
                    "noise_scale": 1.0,
                },
            }

            rfd3_json = {
                "contigmap": {
                    "contigs": [f"{scaffold_range[0]}-{scaffold_range[1]}"],
                    "inpaint_seq": None,
                    "provide_seq": None,
                },
                "inference": {
                    "ligand": f"{ligand_spec.chain_id}{ligand_spec.resnum}",
                    "select_fixed_atoms": {f"{ligand_spec.chain_id}{ligand_spec.resnum}": "ALL"},
                    "select_buried": {f"{ligand_spec.chain_id}{ligand_spec.resnum}": "ALL"},
                },
                "scaffold_length_bounds": list(scaffold_range),
            }

            # Downstream Stage 3 AIzymes Handover Plan
            aizymes_handover = {
                "pipeline_step": "stage3_aizymes_catalytic_activation",
                "strategy": "mutate_contact_shell_to_catalytic_triad",
                "pocket_scan_radius_angstroms": 5.5,
                "recommended_catalytic_triads": [
                    {
                        "mechanism": "cysteine_hydrolase",
                        "residues": ["CYS", "HIS", "ASP"],
                        "target_ligand_atom": chosen_reactive[0] if chosen_reactive else "C1",
                        "ideal_distances": {"CYS:SG": 3.0, "HIS:NE2": 3.4, "ASP:OD1": 2.8},
                    },
                    {
                        "mechanism": "serine_hydrolase",
                        "residues": ["SER", "HIS", "ASP"],
                        "target_ligand_atom": chosen_reactive[0] if chosen_reactive else "C1",
                        "ideal_distances": {"SER:OG": 2.9, "HIS:NE2": 3.0, "ASP:OD1": 2.8},
                    },
                ],
                "instructions": (
                    "Stage 1 will generate a high-affinity binding scaffold for this ligand without pre-placed catalytic side chains. "
                    "Stage 3 (AIzymes) must identify binder pocket positions within 5.5 A of the reactive atom and mutate them to the catalytic triad."
                ),
            }

            metadata = {
                "strategy": "binder_first",
                "design_type": "tight_small_molecule_binder",
                "smiles": smiles,
                "aizymes_evolution_plan": aizymes_handover,
            }

        if out_p:
            with open(out_p / "rfd4_condition_spec.yaml", "w", encoding="utf-8") as f:
                f.write(Case1AnnotatedParser._dump_yaml(rfd4_yaml))
            with open(out_p / "rfd3_config.json", "w", encoding="utf-8") as f:
                json.dump(rfd3_json, f, indent=2)

        return PipelineInput(
            case_type=CaseType.CASE3_SMILES_ONLY,
            ligand=ligand_spec,
            theozyme=theozyme_spec,
            rfd4_condition_spec=rfd4_yaml,
            rfd3_config=rfd3_json,
            metadata=metadata,
        )
