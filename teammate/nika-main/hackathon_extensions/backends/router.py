"""
Central Compute Backend Router for the Enzyme Design Pipeline.
Dispatches jobs across Local, Tamarind Bio, Modal AI, and ETH Euler backends.
Enforces VRAM guardrails and supports seamless dry-run execution and fallbacks.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional, Tuple, Union

from hackathon_extensions.backends.base import ComputeBackend
from hackathon_extensions.backends.eth_cluster_backend import ETHClusterBackend
from hackathon_extensions.backends.local_backend import LocalBackend
from hackathon_extensions.backends.modal_backend import ModalBackend
from hackathon_extensions.backends.tamarind_backend import TamarindBackend
from hackathon_extensions.backends.vram_guardrails import (
    VRAMGuardrailValidator,
    VRAMGuardrailViolation,
)
from hackathon_extensions.config.config_schema import (
    BackendChoice,
    CredentialsConfig,
    PipelineConfig,
)

logger = logging.getLogger(__name__)


class ComputeBackendRouter:
    """
    Central dispatcher managing compute routing, VRAM safety checks,
    and job delegation across available backends.
    """

    def __init__(self, config: Optional[PipelineConfig] = None):
        if config is None:
            # Fallback to default pipeline config
            from hackathon_extensions.config import load_pipeline_config
            self.config = load_pipeline_config()
        else:
            self.config = config

        self.dry_run = self.config.dry_run
        creds = self.config.credentials or CredentialsConfig()

        # Instantiate backend clients
        self.backends: Dict[BackendChoice, ComputeBackend] = {
            BackendChoice.LOCAL: LocalBackend(dry_run=self.dry_run),
            BackendChoice.TAMARIND: TamarindBackend(
                api_key=creds.tamarind_api_key,
                dry_run=self.dry_run,
            ),
            BackendChoice.MODAL: ModalBackend(
                token_id=creds.modal_token_id,
                token_secret=creds.modal_token_secret,
                dry_run=self.dry_run,
            ),
            BackendChoice.ETH_CLUSTER: ETHClusterBackend(
                host=creds.eth_euler_host,
                user=creds.eth_euler_user,
                port=creds.eth_euler_port,
                key_filename=creds.eth_euler_key_filename,
                password=creds.eth_euler_password,
                slurm_partition=creds.eth_euler_slurm_partition,
                scratch_dir=creds.eth_euler_scratch_dir,
                dry_run=self.dry_run,
            ),
        }

        # Verify all configured routes against memory guardrails on init
        self.validate_all_routes()

    def validate_all_routes(self) -> Dict[str, str]:
        """Validates all configured routes against VRAM guardrails."""
        return VRAMGuardrailValidator.validate_pipeline_routes(self.config.compute_routes)

    def get_backend(self, backend_choice: Union[BackendChoice, str]) -> ComputeBackend:
        """Returns initialized backend instance for chosen target."""
        if isinstance(backend_choice, str):
            choice = BackendChoice(backend_choice.lower())
        else:
            choice = backend_choice

        if choice not in self.backends:
            raise ValueError(f"Unknown backend choice: {choice}")
        return self.backends[choice]

    def get_route_for_stage(self, stage_or_gate_name: str) -> BackendChoice:
        """Looks up the assigned backend for a specific pipeline stage or gate."""
        routes = self.config.compute_routes
        norm_name = stage_or_gate_name.lower().replace("-", "_")

        route_mapping = {
            "stage1_rfd4": routes.stage1_rfd4,
            "rfd4": routes.stage1_rfd4,
            "stage1_rfd3": routes.stage1_rfd3,
            "rfd3": routes.stage1_rfd3,
            "stage1_ligandmpnn": routes.stage1_ligandmpnn,
            "ligandmpnn": routes.stage1_ligandmpnn,
            "gate1_cofolding": routes.gate1_cofolding,
            "cofolding": routes.gate1_cofolding,
            "boltz2": routes.gate1_cofolding,
            "gate2_placer": routes.gate2_placer,
            "placer": routes.gate2_placer,
            "gate3_thermo": routes.gate3_thermo,
            "thermompnn": routes.gate3_thermo,
            "gate3_esm3": routes.gate3_esm3,
            "esm3": routes.gate3_esm3,
            "gate4_openmm": routes.gate4_openmm,
            "openmm": routes.gate4_openmm,
            "gate4_fieldtools": routes.gate4_fieldtools,
            "fieldtools": routes.gate4_fieldtools,
            "gate4_orca": routes.gate4_orca,
            "orca": routes.gate4_orca,
            "stage3_aizymes": routes.stage3_aizymes,
            "aizymes": routes.stage3_aizymes,
        }

        for k, v in route_mapping.items():
            if k in norm_name:
                return v if isinstance(v, BackendChoice) else BackendChoice(str(v).lower())

        # Default fallback to local
        return BackendChoice.LOCAL

    def submit_stage_job(
        self,
        stage_or_gate_name: str,
        job_spec: Dict[str, Any],
        override_backend: Optional[Union[BackendChoice, str]] = None,
    ) -> Tuple[str, BackendChoice]:
        """
        Submits job for stage/gate asynchronously after validating VRAM safety.
        Returns (job_id, chosen_backend).
        """
        if override_backend:
            target_backend = (
                override_backend
                if isinstance(override_backend, BackendChoice)
                else BackendChoice(str(override_backend).lower())
            )
        else:
            target_backend = self.get_route_for_stage(stage_or_gate_name)

        # Enforce VRAM guardrails
        VRAMGuardrailValidator.check_route(operation_name=stage_or_gate_name, backend=target_backend)

        # Submit job to target backend
        backend_instance = self.get_backend(target_backend)
        spec = dict(job_spec)
        if "operation" not in spec:
            spec["operation"] = stage_or_gate_name

        job_id = backend_instance.submit_job(spec)
        return job_id, target_backend

    def dispatch_job(
        self,
        stage_or_gate_name: str,
        job_spec: Dict[str, Any],
        output_dir: Optional[str] = None,
        override_backend: Optional[Union[BackendChoice, str]] = None,
    ) -> Dict[str, Any]:
        """
        Dispatches and synchronously waits for job execution on the configured backend.
        Guarantees VRAM safety check before submission.
        """
        if override_backend:
            target_backend = (
                override_backend
                if isinstance(override_backend, BackendChoice)
                else BackendChoice(str(override_backend).lower())
            )
        else:
            target_backend = self.get_route_for_stage(stage_or_gate_name)

        # Enforce VRAM guardrails
        VRAMGuardrailValidator.check_route(operation_name=stage_or_gate_name, backend=target_backend)

        spec = dict(job_spec)
        if "operation" not in spec:
            spec["operation"] = stage_or_gate_name

        backend_instance = self.get_backend(target_backend)
        results = backend_instance.run_synchronous(
            job_spec=spec,
            output_dir=output_dir,
        )

        # Attach execution metadata
        results["dispatched_backend"] = target_backend.value
        results["stage_or_gate"] = stage_or_gate_name
        return results

    def submit_and_poll(
        self,
        engine_name: str,
        job_spec: Dict[str, Any],
        output_dir: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Synchronously submits and polls job on the designated backend."""
        return self.dispatch_job(stage_or_gate_name=engine_name, job_spec=job_spec, output_dir=output_dir)
