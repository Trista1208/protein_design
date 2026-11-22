"""Stage 1 Generative De Novo Enzyme Design Runner.

Provides:
- Stage1Runner: High-level orchestrator executing generative backbone and all-atom
  enzyme design using either RFD4-Proteina (primary default) or RFD3 (+ LigandMPNN).
- VRAM Guardrail Enforcement: Blocks execution on low-memory GPU partitions
  (e.g. ETH Euler RTX 2080 with 8GB VRAM) and routes to high-capacity backends
  (Tamarind Bio or Modal AI).
- Standardized manifest generator producing 'candidates_manifest.json' formatted
  for downstream Stage 2 validation.
- Standalone CLI and Python API interface.
"""

from __future__ import annotations

import argparse
import datetime
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from hackathon_extensions.backends.vram_guardrails import (
    VRAMGuardrailValidator,
    VRAMGuardrailViolation,
)
from hackathon_extensions.inputs.pipeline_input import CaseType, PipelineInput
from hackathon_extensions.stage1_design.base_engine import (
    BaseDesignEngine,
    CandidateDesign,
)
from hackathon_extensions.stage1_design.ligand_mpnn_adapter import LigandMPNNAdapter
from hackathon_extensions.stage1_design.rfd3_adapter import RFD3Adapter
from hackathon_extensions.stage1_design.rfd4_proteina_adapter import RFD4ProteinaAdapter

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
logger = logging.getLogger("Stage1Runner")


class Stage1Runner:
    """Orchestrator for Stage 1 de novo enzyme generative design."""

    def __init__(
        self,
        default_engine: str = "rfd4",
        rfd4_adapter: Optional[RFD4ProteinaAdapter] = None,
        rfd3_adapter: Optional[RFD3Adapter] = None,
        mpnn_adapter: Optional[LigandMPNNAdapter] = None,
    ):
        self.rfd4_adapter = rfd4_adapter or RFD4ProteinaAdapter()
        self.mpnn_adapter = mpnn_adapter or LigandMPNNAdapter()
        self.rfd3_adapter = rfd3_adapter or RFD3Adapter(mpnn_adapter=self.mpnn_adapter)
        self.default_engine = default_engine.lower()

    def select_engine(
        self,
        config: Any = None,
        engine_override: Optional[str] = None
    ) -> BaseDesignEngine:
        """Select appropriate generative engine based on priority:
        1. Explicit override argument
        2. Config stage1 engine setting
        3. Default engine ('rfd4')
        """
        choice = engine_override
        if not choice and config:
            if hasattr(config, "stage1_engine"):
                choice = str(config.stage1_engine)
            elif isinstance(config, dict) and "stage1_engine" in config:
                choice = str(config["stage1_engine"])
            elif isinstance(config, dict) and "engine" in config:
                choice = str(config["engine"])

        choice = (choice or self.default_engine).lower().strip()

        if "rfd3" in choice:
            logger.info("Selected engine: RFD3 with LigandMPNN automated sequence design.")
            return self.rfd3_adapter
        else:
            logger.info("Selected engine: RFD4-Proteina (32D force-field-guided joint all-atom codesign).")
            return self.rfd4_adapter

    def enforce_vram_guardrails(self, engine_name: str, config: Any) -> None:
        """Check compute route against hardware memory safety guardrails.
        
        Blocks all-atom diffusion (RFD4) and structural diffusion (RFD3) from being
        dispatched to low-memory nodes such as the ETH Euler GeForce RTX 2080 (8GB VRAM).
        """
        target_backend = "tamarind"  # Safe default

        if config:
            if hasattr(config, "compute_routes"):
                routes = config.compute_routes
                if "rfd4" in engine_name.lower():
                    target_backend = getattr(routes, "stage1_rfd4", "tamarind")
                else:
                    target_backend = getattr(routes, "stage1_rfd3", "tamarind")
            elif isinstance(config, dict) and "compute_routes" in config:
                routes = config["compute_routes"]
                if "rfd4" in engine_name.lower():
                    target_backend = routes.get("stage1_rfd4", "tamarind")
                else:
                    target_backend = routes.get("stage1_rfd3", "tamarind")

        if hasattr(target_backend, "value"):
            target_backend = target_backend.value
        target_backend_str = str(target_backend).lower()

        try:
            VRAMGuardrailValidator.check_route(operation_name=engine_name, backend=target_backend_str)
        except VRAMGuardrailViolation as e:
            logger.error("\n=======================================================")
            logger.error("!!! CRITICAL HARDWARE VRAM GUARDRAIL TRIGGERED !!!")
            logger.error(str(e))
            logger.error("=======================================================\n")
            raise

    def run(
        self,
        pipeline_input: Union[PipelineInput, Dict[str, Any], str, Path],
        config: Optional[Union[Any, Dict[str, Any], str, Path]] = None,
        backend_router: Any = None,
        engine_override: Optional[str] = None,
        output_dir: Optional[Union[str, Path]] = None,
        dry_run: Optional[bool] = None,
    ) -> Dict[str, Any]:
        """Execute Stage 1 de novo enzyme design end-to-end.
        
        Args:
            pipeline_input: PipelineInput instance, dictionary, or path to input JSON.
            config: PipelineConfig, dictionary, or path to config JSON.
            backend_router: Optional compute router (Tamarind/Modal/Euler).
            engine_override: Explicitly choose 'rfd4' or 'rfd3'.
            output_dir: Destination folder for output artifacts.
            dry_run: Override dry-run setting (defaults to True if not specified).
            
        Returns:
            Dictionary containing execution summary, manifest path, and list of candidates.
        """
        # 1. Ingest PipelineInput
        p_in = self._load_pipeline_input(pipeline_input)
        is_valid, errors = p_in.validate()
        if not is_valid:
            raise ValueError(f"Input validation failed: {'; '.join(errors)}")

        # 2. Ingest Config
        cfg = self._load_config(config)

        # Apply dry_run override if specified
        if dry_run is not None:
            if hasattr(cfg, "dry_run"):
                cfg.dry_run = dry_run
            elif isinstance(cfg, dict):
                cfg["dry_run"] = dry_run

        is_dry = getattr(cfg, "dry_run", True) if not isinstance(cfg, dict) else cfg.get("dry_run", True)

        # 3. Resolve destination directory
        if not output_dir:
            output_dir = getattr(cfg, "output_dir", "outputs/stage1_candidates") if not isinstance(cfg, dict) else cfg.get("output_dir", "outputs/stage1_candidates")
        out_path = Path(output_dir)
        if not is_dry:
            out_path.mkdir(parents=True, exist_ok=True)

        if hasattr(cfg, "output_dir"):
            try:
                setattr(cfg, "output_dir", str(out_path))
            except Exception:
                pass
        elif isinstance(cfg, dict):
            cfg["output_dir"] = str(out_path)

        # 4. Select engine & verify VRAM guardrails
        engine = self.select_engine(config=cfg, engine_override=engine_override)
        if not is_dry:
            self.enforce_vram_guardrails(engine.engine_name, cfg)

        if backend_router is None:
            from hackathon_extensions.backends.router import ComputeBackendRouter
            backend_router = ComputeBackendRouter(config=cfg)

        logger.info("Executing %s with dry_run=%s in directory: %s", engine.engine_name, is_dry, out_path.resolve())

        # 5. Generate scaffolds
        candidates: List[CandidateDesign] = engine.generate_scaffolds(
            pipeline_input=p_in,
            config=cfg,
            backend_router=backend_router,
            dry_run=is_dry,
        )

        # 6. Save Manifest (real runs only; a dry run is validation-only and writes nothing)
        manifest_path = out_path / "candidates_manifest.json"
        if not is_dry:
            manifest_data = {
                "stage": "stage1_generative_design",
                "timestamp": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "engine": engine.engine_name,
                "dry_run": is_dry,
                "total_candidates": len(candidates),
                "output_directory": str(out_path.resolve()),
                "input_summary": {
                    "case_type": p_in.case_type.value if hasattr(p_in.case_type, "value") else str(p_in.case_type),
                    "ligand_resname": p_in.ligand.resname if p_in.ligand else "LIG",
                    "ligand_smiles": p_in.ligand.smiles if p_in.ligand else "",
                    "num_catalytic_residues": len(p_in.theozyme.residues) if p_in.theozyme else 0,
                },
                "candidates": [c.to_dict() for c in candidates],
            }

            with open(manifest_path, "w", encoding="utf-8") as f:
                json.dump(manifest_data, f, indent=2)

            logger.info("Successfully generated %d candidate designs.", len(candidates))
            logger.info("Saved manifest to: %s", manifest_path.resolve())

        return {
            "status": "SUCCESS",
            "engine": engine.engine_name,
            "candidates_count": len(candidates),
            "manifest_path": str(manifest_path.resolve()),
            "candidates": candidates,
        }

    def _load_pipeline_input(self, inp: Union[PipelineInput, Dict[str, Any], str, Path]) -> PipelineInput:
        """Helper to coerce diverse input formats into PipelineInput."""
        if isinstance(inp, PipelineInput):
            return inp
        elif isinstance(inp, (str, Path)):
            path = Path(inp)
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
            # Support wrapped { "data": { ... } } or direct dict
            if "data" in data and "case_type" in data["data"]:
                data = data["data"]
            return PipelineInput.from_dict(data)
        elif isinstance(inp, dict):
            if "data" in inp and "case_type" in inp["data"]:
                inp = inp["data"]
            return PipelineInput.from_dict(inp)
        else:
            raise TypeError(f"Unsupported pipeline_input type: {type(inp)}")

    def _load_config(self, cfg: Optional[Union[Any, Dict[str, Any], str, Path]]) -> Any:
        """Helper to coerce diverse config formats."""
        if cfg is None:
            return {"dry_run": True, "num_designs": 4}
        if isinstance(cfg, (str, Path)):
            path = Path(cfg)
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        return cfg


def main():
    """Command-line entry point for Stage 1 Runner."""
    parser = argparse.ArgumentParser(
        description="Stage 1 De Novo Generative Enzyme Design Runner (RFD4-Proteina / RFD3)."
    )
    parser.add_argument("--input", "-i", type=str, required=True, help="Path to pipeline input JSON.")
    parser.add_argument("--config", "-c", type=str, default=None, help="Path to pipeline config JSON.")
    parser.add_argument("--engine", "-e", type=str, default="rfd4", choices=["rfd4", "rfd3"], help="Generative engine choice (default: rfd4).")
    parser.add_argument("--output", "-o", type=str, default="outputs/stage1_candidates", help="Output directory for generated designs.")
    parser.add_argument("--dry-run", action="store_true", default=False, help="Perform pre-flight validation and generate job spec without submitting real jobs.")

    args = parser.parse_args()

    runner = Stage1Runner(default_engine=args.engine)
    result = runner.run(
        pipeline_input=args.input,
        config=args.config,
        engine_override=args.engine,
        output_dir=args.output,
        dry_run=args.dry_run,
    )

    print("\n=======================================================")
    print(f"Stage 1 Generation Complete: {result['candidates_count']} designs produced.")
    print(f"Manifest written to: {result['manifest_path']}")
    print("=======================================================\n")


if __name__ == "__main__":
    main()
