"""Unit tests for Credentials Manager, Compute Router, and Hardware VRAM Guardrails."""

import unittest

from hackathon_extensions.backends.router import ComputeBackendRouter
from hackathon_extensions.backends.vram_guardrails import (
    VRAMGuardrailValidator,
    VRAMGuardrailViolation,
)
from hackathon_extensions.config.config_schema import (
    BackendChoice,
    ComputeRouteConfig,
    CredentialsConfig,
    PipelineConfig,
)
from hackathon_extensions.config.credentials_manager import CredentialsManager


class TestBackendsAndGuardrails(unittest.TestCase):
    """Test suite for compute routing, authentication, and GPU VRAM constraints."""

    def setUp(self):
        self.creds_mgr = CredentialsManager()

    def test_credentials_manager_dry_run(self):
        """Verify CredentialsManager generates valid mock credentials in dry-run mode."""
        creds = self.creds_mgr.get_credentials(dry_run=True)
        self.assertIsNotNone(creds)
        self.assertIn("mock", creds.tamarind_api_key.lower())
        self.assertIn("mock", creds.modal_token_id.lower())
        self.assertEqual(creds.eth_euler_host, "euler.ethz.ch")

    def test_vram_guardrail_blocks_heavy_models_on_euler(self):
        """Verify strict blocking of high-VRAM models (>=24GB) on ETH Euler 8GB RTX 2080."""
        heavy_ops = [
            "stage1_rfd4",
            "rfd4-proteina",
            "stage1_rfd3",
            "boltz-2",
            "gate1_cofolding",
            "alphafold3",
            "placer_100",
            "gate2_placer",
        ]
        for op in heavy_ops:
            with self.assertRaises(
                VRAMGuardrailViolation,
                msg=f"Expected VRAMGuardrailViolation when scheduling '{op}' on ETH Euler RTX 2080!",
            ):
                VRAMGuardrailValidator.check_route(
                    operation_name=op,
                    backend=BackendChoice.ETH_CLUSTER,
                )

    def test_vram_guardrail_permits_lightweight_models_on_euler(self):
        """Verify that lightweight tasks (<=8GB VRAM) are permitted on ETH Euler RTX 2080."""
        lightweight_ops = [
            "stage1_ligandmpnn",
            "ligandmpnn",
            "gate3_thermo",
            "thermompnn",
            "rosettarelax",
            "gate4_openmm",
            "openmm",
            "gate4_fieldtools",
            "fieldtools",
            "gate4_orca",
            "orca",
            "stage3_aizymes",
            "aizymes",
        ]
        for op in lightweight_ops:
            try:
                VRAMGuardrailValidator.check_route(
                    operation_name=op,
                    backend=BackendChoice.ETH_CLUSTER,
                )
            except VRAMGuardrailViolation as e:
                self.fail(f"Operation '{op}' should be permitted on ETH Euler 8GB node, but was blocked: {e}")

    def test_vram_guardrail_permits_heavy_models_on_tamarind_and_modal(self):
        """Verify that heavy diffusion and co-folding are fully allowed on Tamarind and Modal."""
        for backend in [BackendChoice.TAMARIND, BackendChoice.MODAL]:
            try:
                VRAMGuardrailValidator.check_route("stage1_rfd4", backend=backend)
                VRAMGuardrailValidator.check_route("gate1_cofolding", backend=backend)
                VRAMGuardrailValidator.check_route("gate2_placer", backend=backend)
            except VRAMGuardrailViolation as e:
                self.fail(f"Heavy operations should pass on {backend.value}: {e}")

    def test_default_compute_route_validation(self):
        """Verify that the production default routing table strictly passes VRAM guardrails."""
        default_routes = ComputeRouteConfig()
        validated_mapping = VRAMGuardrailValidator.validate_pipeline_routes(default_routes)
        self.assertEqual(validated_mapping["stage1_rfd4"], "tamarind")
        self.assertEqual(validated_mapping["gate1_cofolding"], "tamarind")
        self.assertEqual(validated_mapping["gate2_placer"], "modal")
        self.assertEqual(validated_mapping["gate3_thermo"], "eth_cluster")
        self.assertEqual(validated_mapping["gate4_openmm"], "eth_cluster")
        self.assertEqual(validated_mapping["stage3_aizymes"], "eth_cluster")

    def test_compute_router_dispatch(self):
        """Verify ComputeBackendRouter instantiates backends and can dispatch dry-run jobs."""
        cfg = PipelineConfig(dry_run=True)
        router = ComputeBackendRouter(config=cfg)

        self.assertIn(BackendChoice.LOCAL, router.backends)
        self.assertIn(BackendChoice.TAMARIND, router.backends)
        self.assertIn(BackendChoice.MODAL, router.backends)
        self.assertIn(BackendChoice.ETH_CLUSTER, router.backends)

        # Dispatch dry-run stage1_rfd4
        job_id, chosen_backend = router.submit_stage_job(
            stage_or_gate_name="stage1_rfd4",
            job_spec={"theozyme": "mock", "ligand": "LIG"},
        )
        self.assertIsNotNone(job_id)
        self.assertEqual(chosen_backend, BackendChoice.TAMARIND)


if __name__ == "__main__":
    unittest.main()
