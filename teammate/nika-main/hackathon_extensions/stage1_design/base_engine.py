"""Base Engine Interface and Abstractions for Stage 1 Generative Protein Design.

Provides:
- BaseDesignEngine: Abstract base class defining the lifecycle for all generative
  scaffold engines (RFD4-Proteina, RFD3, etc.).
- CandidateDesign: Standardized representation of generated candidates.
- Real dry-run validation (payload check & export without fabricating fake structures).
"""

from __future__ import annotations

import abc
import json
import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

logger = logging.getLogger(__name__)

# Standard 1-letter to 3-letter amino acid code mapping
AA_1_TO_3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
AA_3_TO_1 = {v: k for k, v in AA_1_TO_3.items()}


@dataclass
class CandidateDesign:
    """Standardized representation of a generated de novo enzyme design."""

    candidate_id: str
    sequence: str
    pdb_path: str
    fasta_path: Optional[str] = None
    plddt: float = 0.0
    motif_rmsd: float = 0.0
    catalytic_contacts: Dict[str, float] = field(default_factory=dict)
    engine: str = "BaseDesignEngine"
    design_mode: str = "enzyme_design"
    thinking_budget: str = "medium"
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert candidate to dictionary representation."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> CandidateDesign:
        """Instantiate CandidateDesign from a dictionary."""
        return cls(**data)


class BaseDesignEngine(abc.ABC):
    """Abstract Base Class for Stage 1 Generative Scaffold Engines.

    Subclasses implement engine-specific job preparation, remote backend dispatch,
    and output parsing.
    """

    engine_name: str = "BaseDesignEngine"

    def __init__(self, name: Optional[str] = None):
        if name:
            self.engine_name = name

    @abc.abstractmethod
    def validate_input(self, pipeline_input: Any) -> Tuple[bool, List[str]]:
        """Validate whether the given pipeline_input is compatible with this engine."""
        pass

    @abc.abstractmethod
    def prepare_job_spec(self, pipeline_input: Any, config: Any) -> Dict[str, Any]:
        """Convert high-level pipeline input and configuration into an engine-specific job spec."""
        pass

    @abc.abstractmethod
    def parse_outputs(self, output_data: Dict[str, Any], output_dir: Union[str, Path]) -> List[CandidateDesign]:
        """Parse raw engine output data into standardized CandidateDesign objects."""
        pass

    def generate_scaffolds(
        self,
        pipeline_input: Any,
        config: Any,
        backend_router: Any = None,
        dry_run: Optional[bool] = None,
    ) -> List[CandidateDesign]:
        """Generate scaffolds using this engine.

        If dry_run is True, validates configuration and writes the exact payload
        to disk without submitting GPU compute or generating fake structures.
        """
        # Validate inputs
        is_valid, error_msgs = self.validate_input(pipeline_input)
        if not is_valid:
            raise ValueError(f"Input validation failed for {self.engine_name}: {'; '.join(error_msgs)}")

        # Determine dry-run status
        is_dry_run = dry_run
        if is_dry_run is None:
            if hasattr(config, "dry_run"):
                is_dry_run = config.dry_run
            elif isinstance(config, dict):
                is_dry_run = config.get("dry_run", False)
            else:
                is_dry_run = False

        # Destination directory
        output_dir = getattr(config, "output_dir", None)
        if not output_dir and isinstance(config, dict):
            output_dir = config.get("output_dir")
        if not output_dir:
            output_dir = "outputs/stage1_candidates"
        out_path = Path(output_dir)
        out_path.mkdir(parents=True, exist_ok=True)

        job_spec = self.prepare_job_spec(pipeline_input, config)

        if backend_router is None:
            raise RuntimeError(
                f"[{self.engine_name}] Backend router is required for execution. "
                "Specify a compute backend (e.g., Tamarind or Modal)."
            )

        # Route through ComputeBackendRouter.
        # Dry-run takes this EXACT same path with the EXACT same job_spec; the only
        # difference is that the backend stops after Tamarind POST /api/validate-job
        # instead of continuing to /submit-job. No separate dry-run code path, no
        # dry-run-only artifacts.
        logger.info(
            "[%s] Submitting job to compute backend router (dry_run=%s).",
            self.engine_name,
            is_dry_run,
        )
        backend_response = backend_router.submit_and_poll(
            engine_name=self.engine_name,
            job_spec=job_spec,
        )
        return self.parse_outputs(backend_response, out_path)

