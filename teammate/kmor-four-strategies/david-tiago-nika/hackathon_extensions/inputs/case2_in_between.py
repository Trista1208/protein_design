"""Case 2 Parser: Intermediate Inputs (Holo Active Site Crops & SMILES with Known Mechanism).

Handles the two 'in-between' scenarios encountered in enzyme engineering:
- Subcase 2a: Unannotated holo-structure or active site crop:
    * Strips non-catalytic scaffold residues for pure de novo re-scaffolding, OR
    * Retains scaffold and configures partial diffusion masks (redesigning loops / pocket shells
      while freezing catalytic triad and ligand).
- Subcase 2b: SMILES + known catalytic mechanism (no starting 3D protein structure):
    * Runs ConformerGenerator to produce 3D ligand with unique atom names.
    * Invokes TheozymeBuilder to construct 3D active site geometry (Bürgi-Dunitz trajectory, dyads/triads).
    * Emits RFD4 condition-spec YAML and RFD3 JSON configs.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from hackathon_extensions.inputs.case1_annotated import Case1AnnotatedParser
from hackathon_extensions.inputs.conformer_generator import ConformerGenerator
from hackathon_extensions.inputs.pipeline_input import (
    CaseType,
    LigandSpec,
    PipelineInput,
    TheozymeResidueSpec,
    TheozymeSpec,
)
from hackathon_extensions.inputs.theozyme_builder import TheozymeBuilder


class Case2InBetweenParser:
    """Parser and generator for Case 2a (holo-structures/crops) and Case 2b (SMILES + mechanism)."""

    def __init__(
        self,
        conformer_generator: Optional[ConformerGenerator] = None,
        theozyme_builder: Optional[TheozymeBuilder] = None,
        default_scaffold_length_range: Tuple[int, int] = (150, 300),
    ):
        self.conformer_generator = conformer_generator or ConformerGenerator()
        self.theozyme_builder = theozyme_builder or TheozymeBuilder()
        self.default_scaffold_length_range = default_scaffold_length_range

    # -------------------------------------------------------------------------
    # Subcase 2a: Unannotated Holo-Structure or Active Site Crop
    # -------------------------------------------------------------------------

    def parse_case2a_holo(
        self,
        holo_structure: Union[str, Path],
        catalytic_residues: List[Dict[str, Any]],
        ligand_chain: str = "Z",
        ligand_resname: str = "LIG",
        mode: str = "strip_scaffold",  # "strip_scaffold" or "partial_diffusion"
        loops_to_diffuse: Optional[List[Tuple[str, int, int]]] = None,  # [(chain, start, end)]
        partial_diffusion_timesteps: int = 20,
        scaffold_length_range: Optional[Tuple[int, int]] = None,
        output_dir: Optional[Union[str, Path]] = None,
    ) -> PipelineInput:
        """Parse holo-structure and configure de novo scaffolding or partial diffusion.

        Parameters
        ----------
        holo_structure : Union[str, Path]
            PDB or mmCIF file path or text.
        catalytic_residues : List[Dict[str, Any]]
            User-identified catalytic residues, e.g.:
            [{"chain": "A", "resnum": 25, "resname": "CYS"}, {"chain": "A", "resnum": 162, "resname": "HIS"}]
        ligand_chain : str
            Ligand chain ID (default "Z").
        ligand_resname : str
            Ligand residue name (default "LIG").
        mode : str
            "strip_scaffold" to build a brand new de novo scaffold around catalytic residues,
            "partial_diffusion" to retain backbone and diffuse designated loop regions.
        loops_to_diffuse : Optional[List[Tuple[str, int, int]]]
            Residue index ranges to unfreeze/diffuse, e.g. [("A", 45, 60), ("A", 120, 135)].
        partial_diffusion_timesteps : int
            Number of diffusion steps for loop remodeling (typically 15-25).
        scaffold_length_range : Optional[Tuple[int, int]]
            Scaffold length bounds for de novo generation.
        output_dir : Optional[Union[str, Path]]
            Destination directory to save outputs.

        Returns
        -------
        PipelineInput
            Configured input object for generative design.
        """
        raw_text = Case1AnnotatedParser._load_structure_text(holo_structure)
        is_cif = "_atom_site" in raw_text or str(holo_structure).endswith((".cif", ".mmcif"))

        annot_parser = Case1AnnotatedParser()
        protein_atoms, ligand_atoms = annot_parser._parse_atoms(
            raw_text=raw_text,
            is_cif=is_cif,
            ligand_chain=ligand_chain,
            ligand_resnum=1,
            ligand_resname=ligand_resname,
        )

        if not ligand_atoms:
            protein_atoms, ligand_atoms = annot_parser._auto_locate_ligand(raw_text, is_cif)

        # Parse catalytic residue specs
        cat_specs: List[TheozymeResidueSpec] = []
        cat_keys = set()
        for c in catalytic_residues:
            chain = c.get("chain", "A")
            resnum = int(c.get("resnum", 1))
            resname = c.get("resname", "CYS")
            fixed = c.get("fixed_atoms", ["SG", "CB", "CA"] if resname == "CYS" else ["NE2", "ND1", "CA"])
            cat_specs.append(TheozymeResidueSpec(
                resname=resname,
                chain_id=chain,
                resnum=resnum,
                target_atoms=fixed,
                role=c.get("role", "catalytic_residue"),
            ))
            cat_keys.add((chain, resnum))

        ligand_pdb = annot_parser._atoms_to_pdb(ligand_atoms)
        reactive_atoms = annot_parser._find_reactive_ligand_atoms(ligand_atoms)
        ligand_spec = LigandSpec(
            conformer_pdb=ligand_pdb,
            resname=ligand_resname,
            chain_id=ligand_chain,
            resnum=1,
            reactive_atoms=reactive_atoms,
            atom_metadata={"source": "holo_structure", "atom_count": len(ligand_atoms)},
        )

        scaffold_range = scaffold_length_range or self.default_scaffold_length_range

        if mode == "strip_scaffold":
            # Strip all protein atoms EXCEPT catalytic residues
            catalytic_protein_atoms = [
                a for a in protein_atoms if (a["chain"], a["resnum"]) in cat_keys
            ]
            stripped_theozyme_pdb = annot_parser._atoms_to_pdb(catalytic_protein_atoms + ligand_atoms)

            theozyme_spec = TheozymeSpec(
                residues=cat_specs,
                description="Stripped active site crop from holo-enzyme for de novo scaffolding",
                mechanism_type="case2a_stripped_crop",
                theozyme_pdb=stripped_theozyme_pdb,
            )

            rfd4_yaml = annot_parser._build_rfd4_condition_spec(cat_specs, ligand_spec, scaffold_range)
            rfd3_json = annot_parser._build_rfd3_config(cat_specs, ligand_spec, scaffold_range)
            diffusion_mask = None

        elif mode == "partial_diffusion":
            # Retain entire scaffold and build partial diffusion mask
            all_residues_set = {(a["chain"], a["resnum"]) for a in protein_atoms}
            mask_entries = []

            diffuse_ranges = loops_to_diffuse or []
            for chain, resnum in sorted(all_residues_set):
                is_cat = (chain, resnum) in cat_keys
                is_diffuse_loop = any(
                    chain == d_ch and d_start <= resnum <= d_end
                    for d_ch, d_start, d_end in diffuse_ranges
                )
                mask_entries.append({
                    "chain": chain,
                    "resnum": resnum,
                    "fixed": is_cat or (not is_diffuse_loop and len(diffuse_ranges) > 0),
                    "diffuse": is_diffuse_loop,
                })

            diffusion_mask = {
                "mode": "partial_diffusion",
                "diffusion_timesteps": partial_diffusion_timesteps,
                "loops_diffused": diffuse_ranges,
                "residue_mask": mask_entries,
            }

            theozyme_spec = TheozymeSpec(
                residues=cat_specs,
                description="Holo-scaffold with partial diffusion mask",
                mechanism_type="case2a_partial_diffusion",
                theozyme_pdb=annot_parser._atoms_to_pdb(protein_atoms + ligand_atoms),
            )

            rfd4_yaml = annot_parser._build_rfd4_condition_spec(cat_specs, ligand_spec, scaffold_range)
            rfd4_yaml["partial_diffusion"] = diffusion_mask

            rfd3_json = annot_parser._build_rfd3_config(cat_specs, ligand_spec, scaffold_range)
            rfd3_json["diffuser"] = {
                "partial_T": partial_diffusion_timesteps,
            }
        else:
            raise ValueError(f"Unknown mode: {mode}. Must be 'strip_scaffold' or 'partial_diffusion'.")

        if output_dir:
            out_p = Path(output_dir)
            out_p.mkdir(parents=True, exist_ok=True)
            with open(out_p / "theozyme_core.pdb", "w", encoding="utf-8") as f:
                f.write(theozyme_spec.theozyme_pdb or "")
            with open(out_p / "rfd4_condition_spec.yaml", "w", encoding="utf-8") as f:
                f.write(Case1AnnotatedParser._dump_yaml(rfd4_yaml))
            with open(out_p / "rfd3_config.json", "w", encoding="utf-8") as f:
                json.dump(rfd3_json, f, indent=2)

        return PipelineInput(
            case_type=CaseType.CASE2A_UNANNOTATED_HOLO,
            ligand=ligand_spec,
            theozyme=theozyme_spec,
            scaffold_source_pdb=raw_text,
            partial_diffusion_mask=diffusion_mask,
            rfd4_condition_spec=rfd4_yaml,
            rfd3_config=rfd3_json,
            metadata={"mode": mode, "num_catalytic_residues": len(cat_specs)},
        )

    # -------------------------------------------------------------------------
    # Subcase 2b: SMILES + Known Catalytic Mechanism (No 3D Protein Structure)
    # -------------------------------------------------------------------------

    def parse_case2b_smiles_mechanism(
        self,
        smiles: str,
        mechanism_type: str = "cysteine_hydrolase",
        custom_residues: Optional[List[TheozymeResidueSpec]] = None,
        reactive_atoms: Optional[List[str]] = None,
        scaffold_length_range: Optional[Tuple[int, int]] = None,
        output_dir: Optional[Union[str, Path]] = None,
    ) -> PipelineInput:
        """Construct 3D active site de novo from SMILES and mechanism type.

        Parameters
        ----------
        smiles : str
            Ligand / substrate SMILES string.
        mechanism_type : str
            Catalytic mechanism: "cysteine_hydrolase", "serine_hydrolase", "kemp_eliminase".
        custom_residues : Optional[List[TheozymeResidueSpec]]
            Explicit residue parameters if customizing distances/angles.
        reactive_atoms : Optional[List[str]]
            Ligand reactive atom names (e.g. ["C1"]).
        scaffold_length_range : Optional[Tuple[int, int]]
            Scaffold length range (default 150-300).
        output_dir : Optional[Union[str, Path]]
            Output directory for generated files.

        Returns
        -------
        PipelineInput
            PipelineInput ready for Stage 1 RFD4/RFD3 generation.
        """
        out_p = Path(output_dir) if output_dir else None
        lig_pdb_path = (out_p / "ligand_conformer.pdb") if out_p else None

        # Step 1: Generate 3D conformer with unique atom names
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

        # Step 2: Build 3D theozyme geometry
        theozyme_pdb_path = (out_p / "theozyme_input.pdb") if out_p else None
        theozyme_cif_path = (out_p / "theozyme_input.cif") if out_p else None
        theozyme_cst_path = (out_p / "enzdes.cst") if out_p else None

        theozyme_spec = self.theozyme_builder.build_theozyme(
            ligand=ligand_spec,
            mechanism_type=mechanism_type,
            custom_residues=custom_residues,
            output_pdb=theozyme_pdb_path,
            output_cif=theozyme_cif_path,
            output_cst=theozyme_cst_path,
        )

        # Step 3: Build RFD4 condition-spec & RFD3 config
        scaffold_range = scaffold_length_range or self.default_scaffold_length_range
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

        if out_p:
            with open(out_p / "rfd4_condition_spec.yaml", "w", encoding="utf-8") as f:
                f.write(Case1AnnotatedParser._dump_yaml(rfd4_yaml))
            with open(out_p / "rfd3_config.json", "w", encoding="utf-8") as f:
                json.dump(rfd3_json, f, indent=2)

        return PipelineInput(
            case_type=CaseType.CASE2B_SMILES_MECHANISM,
            ligand=ligand_spec,
            theozyme=theozyme_spec,
            rfd4_condition_spec=rfd4_yaml,
            rfd3_config=rfd3_json,
            metadata={
                "mechanism_type": mechanism_type,
                "smiles": smiles,
                "conformer_atoms": meta.get("num_atoms", 0),
            },
        )
