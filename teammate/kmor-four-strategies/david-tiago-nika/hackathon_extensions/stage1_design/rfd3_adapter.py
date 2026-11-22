"""RFD3 Generative Design Adapter with Automated LigandMPNN Interfacing.

Implements structural scaffolding with RFD3 (RFdiffusion3):
- Generates contig definitions (e.g. '[15-35/A25-25/25-45/A68-68/15-35]').
- Specifies fixed catalytic residue motifs and ligand chain Z constraints.
- Configures 'infer_ori_strategy: "com"' and 'select_fixed_atoms'.
- Emits backbone-only 3D structures.
- Automatically invokes LigandMPNNAdapter to design sequences and pack all-atom
  side chains onto the de novo backbones, outputting complete candidate designs.
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from hackathon_extensions.stage1_design.base_engine import (
    BaseDesignEngine,
    CandidateDesign,
)
from hackathon_extensions.stage1_design.ligand_mpnn_adapter import (
    LigandMPNNAdapter,
    LigandMPNNConfig,
)

logger = logging.getLogger(__name__)


@dataclass
class RFD3Config:
    """Configuration options for RFD3 scaffolding."""
    contig: Optional[str] = None
    infer_ori_strategy: str = "com"  # Center of mass orientation strategy
    select_fixed_atoms: str = "CA,N,C,O"  # Fixed backbone atoms for catalytic residues
    ligand_chain: str = "Z"
    diffusion_steps: int = 50
    noise_level: float = 0.0
    num_designs: int = 4
    run_ligandmpnn_sequence_design: bool = True
    mpnn_config: LigandMPNNConfig = field(default_factory=LigandMPNNConfig)

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["mpnn_config"] = self.mpnn_config.to_dict()
        return d


class RFD3Adapter(BaseDesignEngine):
    """Adapter for RFD3 (RFdiffusion3) structural scaffolding with LigandMPNN sequence design."""

    engine_name: str = "RFD3"

    def __init__(
        self,
        name: str = "RFD3",
        mpnn_adapter: Optional[LigandMPNNAdapter] = None,
    ):
        super().__init__(name=name)
        self.mpnn_adapter = mpnn_adapter or LigandMPNNAdapter()

    def validate_input(self, pipeline_input: Any) -> Tuple[bool, List[str]]:
        """Validate input for RFD3 scaffolding."""
        errors: List[str] = []

        if pipeline_input is None:
            return False, ["Pipeline input cannot be None."]

        ligand = getattr(pipeline_input, "ligand", None)
        if not ligand and isinstance(pipeline_input, dict):
            ligand = pipeline_input.get("ligand")

        if not ligand:
            errors.append("RFD3 requires a ligand specification.")

        theozyme = getattr(pipeline_input, "theozyme", None)
        if not theozyme and isinstance(pipeline_input, dict):
            theozyme = pipeline_input.get("theozyme")

        if not theozyme:
            # If no theozyme, check if contig, length, or unindex was explicitly provided
            rfd3_cfg = getattr(pipeline_input, "rfd3_config", {})
            if isinstance(pipeline_input, dict):
                rfd3_cfg = pipeline_input.get("rfd3_config", {})
            if not rfd3_cfg.get("contig") and not rfd3_cfg.get("length") and not rfd3_cfg.get("unindex"):
                errors.append("RFD3 requires either an annotated theozyme, an explicit contig definition, or length/unindex specs.")

        return len(errors) == 0, errors

    def format_contig(self, theozyme: Any, custom_contig: Optional[str] = None) -> str:
        """Construct standard RFdiffusion contig definition string.
        
        Example: '[15-35/A25-25/25-45/A68-68/25-45/A102-102/15-35]'
        """
        if custom_contig:
            return custom_contig

        residues = getattr(theozyme, "residues", []) if theozyme else []
        if isinstance(theozyme, dict):
            residues = theozyme.get("residues", [])

        if not residues:
            return "[120-150]"

        # Default slot positions if not specified in residues
        default_slots = [25, 68, 102]
        elements = []
        elements.append("15-35")

        for idx, r in enumerate(residues):
            r_num = getattr(r, "resnum", None) if not isinstance(r, dict) else r.get("resnum")
            chain = getattr(r, "chain_id", "A") if not isinstance(r, dict) else r.get("chain_id", "A")
            slot = r_num if (r_num and r_num > 0) else default_slots[idx % len(default_slots)]
            elements.append(f"{chain}{slot}-{slot}")
            elements.append("25-45" if idx < len(residues) - 1 else "15-35")

        return f"[{'/'.join(elements)}]"

    def prepare_job_spec(self, pipeline_input: Any, config: Any) -> Dict[str, Any]:
        """Prepare job spec for RFD3 backbone generation."""
        theozyme = getattr(pipeline_input, "theozyme", None)
        if not theozyme and isinstance(pipeline_input, dict):
            theozyme = pipeline_input.get("theozyme")

        rfd3_custom = getattr(pipeline_input, "rfd3_config", {})
        if not rfd3_custom and isinstance(pipeline_input, dict):
            rfd3_custom = pipeline_input.get("rfd3_config", {})

        contig = self.format_contig(theozyme, rfd3_custom.get("contig"))

        num_designs = 1
        if hasattr(config, "num_designs"):
            num_designs = int(config.num_designs)
        elif isinstance(config, dict) and "num_designs" in config:
            num_designs = int(config["num_designs"])

        ligand = getattr(pipeline_input, "ligand", None)
        if not ligand and isinstance(pipeline_input, dict):
            ligand = pipeline_input.get("ligand")

        raw_res = str(rfd3_custom.get("unindex") or rfd3_custom.get("interfaceResidues") or "")
        import re
        clean_interface = " ".join(re.findall(r"\d+", raw_res))

        job_spec = {
            "engine": self.engine_name,
            "contig": contig,
            "infer_ori_strategy": rfd3_custom.get("infer_ori_strategy", "com"),
            "select_fixed_atoms": rfd3_custom.get("select_fixed_atoms", "CA,N,C,O"),
            "ligand_chain": "Z",
            "diffusion_steps": rfd3_custom.get("diffusion_steps", 50),
            "num_designs": num_designs,
            "joint_sequence_structure_codesign": False,
            "requires_separate_sequence_design": True,
            "ligand": ligand.to_dict() if hasattr(ligand, "to_dict") else ligand,
            "theozyme": theozyme.to_dict() if hasattr(theozyme, "to_dict") else theozyme,
            "tamarind_api_payload": {
                "type": "rfdiffusion3",
                "settings": {
                    "task": "Custom Contigs",
                    "pdbFile": rfd3_custom.get("pdbFile") or (Path(rfd3_custom["input"]).name if "input" in rfd3_custom else None),
                    "contigs": contig,
                    "numDesigns": num_designs,
                }
            },
            "modal_job_spec": {
                "app": "enzyme-design-rfd3",
                "function": "generate_backbones",
                "kwargs": {
                    "contig": contig,
                    "infer_ori_strategy": "com",
                    "num_designs": num_designs,
                }
            }
        }
        return job_spec

    def generate_scaffolds(
        self,
        pipeline_input: Any,
        config: Any,
        backend_router: Any = None,
        dry_run: bool = False,
    ) -> List[CandidateDesign]:
        """Generate backbone scaffolds using RFD3, then run LigandMPNN sequence design.
        
        This bridges RFD3's backbone-only generation with LigandMPNN to produce complete
        all-atom designs ready for Stage 2 validation.
        """
        is_valid, errors = self.validate_input(pipeline_input)
        if not is_valid:
            raise ValueError(f"Invalid input for engine {self.engine_name}: {'; '.join(errors)}")

        config_dry_run = False
        if hasattr(config, "dry_run"):
            config_dry_run = bool(config.dry_run)
        elif isinstance(config, dict):
            config_dry_run = bool(config.get("dry_run", False))

        is_dry_run = dry_run or config_dry_run

        output_dir = getattr(config, "output_dir", None)
        if not output_dir and isinstance(config, dict):
            output_dir = config.get("output_dir")
        if not output_dir:
            output_dir = "outputs/stage1_candidates"
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        job_spec = self.prepare_job_spec(pipeline_input, config)

        if backend_router is None:
            raise RuntimeError(f"[{self.engine_name}] Backend router is required for execution.")

        # 1. Generate backbones.
        # Dry-run follows this identical path with the identical job_spec; the backend
        # stops after Tamarind POST /api/validate-job instead of /submit-job.
        logger.info(
            "[%s] Submitting backbone generation to backend router (dry_run=%s).",
            self.engine_name,
            is_dry_run,
        )
        backend_resp = backend_router.submit_and_poll(
            engine_name=self.engine_name,
            job_spec=job_spec,
        )
        backbone_pdbs: List[Path] = [Path(item.get("pdb_path")) for item in backend_resp.get("designs", []) if item.get("pdb_path")]

        # 2. Automatically run LigandMPNN sequence design onto backbones
        logger.info("[%s] Interfacing with LigandMPNN to design sequences onto %d backbones.", self.engine_name, len(backbone_pdbs))
        theozyme = getattr(pipeline_input, "theozyme", None)
        if not theozyme and isinstance(pipeline_input, dict):
            theozyme = pipeline_input.get("theozyme")

        designed_candidates = self.mpnn_adapter.design_sequences_for_backbones(
            backbone_pdbs=backbone_pdbs,
            theozyme=theozyme,
            output_dir=out_path,
            backend_router=backend_router,
            dry_run=is_dry_run,
        )

        return designed_candidates

    def parse_outputs(self, output_data: Dict[str, Any], output_dir: Union[str, Path]) -> List[CandidateDesign]:
        """Parse raw output dictionary into CandidateDesign instances."""
        out_path = Path(output_dir)
        designs: List[CandidateDesign] = []
        raw_designs = output_data.get("designs", [])

        for item in raw_designs:
            candidate = CandidateDesign.from_dict(item)
            candidate.engine = "RFD3+LigandMPNN"
            designs.append(candidate)

        return designs
