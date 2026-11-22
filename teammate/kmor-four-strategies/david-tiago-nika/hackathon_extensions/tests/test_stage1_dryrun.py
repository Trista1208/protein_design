"""Unit tests for Stage 1 De Novo Generative Design dry-run execution."""

import json
import tempfile
import unittest
from pathlib import Path

from hackathon_extensions.inputs.case1_annotated import Case1AnnotatedParser
from hackathon_extensions.backends.vram_guardrails import VRAMGuardrailViolation
from hackathon_extensions.config.config_schema import BackendChoice, ComputeRouteConfig, PipelineConfig
from hackathon_extensions.inputs.pipeline_input import (
    CaseType,
    LigandSpec,
    PipelineInput,
    TheozymeResidueSpec,
    TheozymeSpec,
)
from hackathon_extensions.stage1_design.ligand_mpnn_adapter import LigandMPNNAdapter
from hackathon_extensions.stage1_design.rfd3_adapter import RFD3Adapter
from hackathon_extensions.stage1_design.rfd4_proteina_adapter import RFD4ProteinaAdapter
from hackathon_extensions.stage1_design.runner import Stage1Runner


from hackathon_extensions.backends.router import ComputeBackendRouter

class TestStage1DryRun(unittest.TestCase):
    """Test suite for Stage 1 de novo generative design adapters and runner."""

    def setUp(self):
        # Real annotated input (1QJI theozyme): the Tamarind validator checks that the
        # referenced structure exists, so a mock input without a file cannot validate.
        repo_root = Path(__file__).resolve().parents[2]
        task_json = repo_root / "test_pipeline" / "1qji_json.json"
        self.pipeline_input = Case1AnnotatedParser().parse_rfd3_json(task_json)
        self.cfg = PipelineConfig(dry_run=True)

    def test_rfd4_proteina_adapter_dryrun(self):
        """Verify RFD4ProteinaAdapter payload preparation and mock generation."""
        adapter = RFD4ProteinaAdapter()
        job_spec = adapter.prepare_job_spec(
            pipeline_input=self.pipeline_input,
            config=self.cfg,
        )
        self.assertIn("tamarind_api_payload", job_spec)
        self.assertIn("thinking_budget", job_spec)

        # Generate candidates in dry-run mode (pre-flight validation, no fake PDBs)
        with tempfile.TemporaryDirectory() as tmp_dir:
            config_dict = {"dry_run": True, "num_designs": 2, "output_dir": tmp_dir}
            candidates = adapter.generate_scaffolds(
                pipeline_input=self.pipeline_input,
                config=config_dict,
                backend_router=ComputeBackendRouter(config=self.cfg),
                dry_run=True,
            )
            self.assertEqual(len(candidates), 0)
            # Dry-run is validation-only (Tamarind /api/validate-job): nothing is written.
            self.assertEqual(list(Path(tmp_dir).iterdir()), [])

    def test_rfd3_and_ligandmpnn_adapter_dryrun(self):
        """Verify RFD3Adapter pre-flight validation in dry-run mode."""
        mpnn = LigandMPNNAdapter()
        rfd3 = RFD3Adapter(mpnn_adapter=mpnn)

        with tempfile.TemporaryDirectory() as tmp_dir:
            config_dict = {"dry_run": True, "num_designs": 2, "output_dir": tmp_dir}
            candidates = rfd3.generate_scaffolds(
                pipeline_input=self.pipeline_input,
                config=config_dict,
                backend_router=ComputeBackendRouter(config=self.cfg),
                dry_run=True,
            )
            self.assertEqual(len(candidates), 0)
            self.assertEqual(list(Path(tmp_dir).iterdir()), [])

    def test_stage1_runner_full_dryrun(self):
        """Verify Stage1Runner end-to-end execution and manifest export."""
        runner = Stage1Runner(default_engine="rfd4")
        cfg = PipelineConfig(dry_run=True)

        with tempfile.TemporaryDirectory() as tmp_dir:
            result = runner.run(
                pipeline_input=self.pipeline_input,
                config=cfg,
                engine_override="rfd4",
                output_dir=tmp_dir,
                dry_run=True,
            )
            self.assertEqual(result["status"], "SUCCESS")
            self.assertEqual(result["candidates_count"], 0)
            # Dry-run is validation-only: no manifest, no artifacts.
            self.assertFalse(Path(result["manifest_path"]).exists())
            self.assertEqual(list(Path(tmp_dir).iterdir()), [])

    def test_stage1_runner_enforces_vram_guardrail(self):
        """Verify that routing RFD4 diffusion to ETH Euler RTX 2080 raises VRAMGuardrailViolation."""
        runner = Stage1Runner(default_engine="rfd4")
        bad_config = PipelineConfig(
            compute_routes=ComputeRouteConfig(stage1_rfd4=BackendChoice.ETH_CLUSTER)
        )
        with self.assertRaises(VRAMGuardrailViolation):
            runner.enforce_vram_guardrails("stage1_rfd4", bad_config)


if __name__ == "__main__":
    unittest.main()
