"""
VRAM Guardrails for compute routing.
Prevents catastrophic CUDA Out-Of-Memory (OOM) failures on low-VRAM nodes
(e.g., ETH Euler GeForce RTX 2080 nodes with only 8GB VRAM) by strictly blocking
heavy all-atom diffusion, large-scale co-folding, and deep ensemble sampling.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional, Set, Union

from hackathon_extensions.config.config_schema import BackendChoice, ComputeRouteConfig

logger = logging.getLogger(__name__)


class VRAMGuardrailViolation(Exception):
    """Raised when an operation exceeding hardware VRAM capacity is scheduled on a restricted node."""
    pass


class VRAMGuardrailValidator:
    """
    Enforces hardware memory safety constraints across execution targets.

    Specifically:
    ETH Zurich Euler cluster GPU partitions commonly feature NVIDIA GeForce RTX 2080 GPUs
    with only 8.0 GB of VRAM. Attempting all-atom diffusion (RFD4-Proteina) or heavy
    co-folding (Boltz-2 / AlphaFold3 / Protenix-v1) or 100-sample PLACER ensembles on 8GB
    causes fatal CUDA out-of-memory kernel termination.
    """

    # Estimated minimum VRAM required per operation category (in Gigabytes)
    OPERATION_MIN_VRAM_GB: Dict[str, float] = {
        # High-VRAM operations (>= 24GB required -> Tamarind Bio or Modal A100/H100)
        "rfd4": 24.0,
        "rfd4-proteina": 24.0,
        "stage1_rfd4": 24.0,
        "rfd3": 24.0,
        "stage1_rfd3": 24.0,
        "boltz2": 24.0,
        "boltz-2": 24.0,
        "af3": 24.0,
        "alphafold3": 24.0,
        "protenix-v1": 24.0,
        "gate1_cofolding": 24.0,
        "placer_100": 24.0,
        "placer_ensemble": 24.0,
        "gate2_placer": 24.0,
        "esm3": 24.0,
        "gate3_esm3": 24.0,

        # Low-VRAM operations (<= 8GB compatible -> ETH Euler RTX 2080 / Local)
        "ligandmpnn": 4.0,
        "stage1_ligandmpnn": 4.0,
        "thermompnn": 4.0,
        "gate3_thermo": 4.0,
        "rosettarelax": 2.0,
        "openmm": 6.0,
        "gate4_openmm": 6.0,
        "fieldtools": 2.0,
        "gate4_fieldtools": 2.0,
        "orca": 4.0,
        "gate4_orca": 4.0,
        "stage3_aizymes": 4.0,
        "aizymes": 4.0,
    }

    # Strict blocked operation identifiers for 8GB nodes
    BLOCKED_ON_8GB: Set[str] = {
        "rfd4", "rfd4-proteina", "stage1_rfd4",
        "rfd3", "stage1_rfd3",
        "boltz2", "boltz-2", "af3", "alphafold3", "protenix-v1", "chai-1", "gate1_cofolding",
        "placer_100", "gate2_placer", "placer_ensemble",
    }

    # Explicitly allowed operations on 8GB nodes (ETH Euler RTX 2080)
    ALLOWED_ON_8GB: Set[str] = {
        "ligandmpnn", "stage1_ligandmpnn",
        "thermompnn", "gate3_thermo",
        "rosettarelax",
        "openmm", "gate4_openmm",
        "fieldtools", "gate4_fieldtools",
        "orca", "gate4_orca",
        "stage3_aizymes", "aizymes",
    }

    @classmethod
    def get_required_vram(cls, operation_name: str) -> float:
        """Returns estimated minimum VRAM in GB for an operation."""
        op_norm = operation_name.lower().replace(" ", "_").replace("-", "_")
        for k, v in cls.OPERATION_MIN_VRAM_GB.items():
            if k.replace("-", "_") in op_norm:
                return v
        return 8.0  # Conservative default

    @classmethod
    def check_route(
        cls,
        operation_name: str,
        backend: Union[BackendChoice, str],
        hardware_vram_limit_gb: Optional[float] = None
    ) -> None:
        """
        Validates whether operation is memory-safe for target backend.
        Raises VRAMGuardrailViolation if unsafe.
        """
        backend_str = backend.value if isinstance(backend, BackendChoice) else str(backend).lower()
        op_norm = operation_name.lower().replace(" ", "_").replace("-", "_")

        # Determine target node VRAM limit:
        # ETH Euler default GPU partition has 8GB RTX 2080.
        # Local laptop may have 8GB.
        # Tamarind and Modal provide high-VRAM A100/H100 (40GB - 80GB).
        if hardware_vram_limit_gb is not None:
            effective_vram = hardware_vram_limit_gb
        elif backend_str == BackendChoice.ETH_CLUSTER.value:
            effective_vram = 8.0
        elif backend_str in (BackendChoice.TAMARIND.value, BackendChoice.MODAL.value):
            effective_vram = 80.0
        else:
            effective_vram = 16.0  # Local default

        is_constrained = effective_vram <= 8.0

        if not is_constrained:
            # High-VRAM targets (Tamarind REST, Modal A100/H100) are unconstrained
            return

        # Check if operation is in the blocked list or requires > effective_vram
        is_blocked = any(b.replace("-", "_") in op_norm for b in cls.BLOCKED_ON_8GB)
        req_vram = cls.get_required_vram(op_norm)

        if is_blocked or req_vram > effective_vram:
            raise VRAMGuardrailViolation(
                f"[VRAM Guardrail Violation] Operation '{operation_name}' requires ~{req_vram:.1f} GB VRAM, "
                f"which exceeds the {effective_vram:.1f} GB capacity of backend '{backend_str}' "
                f"(ETH Zurich Euler RTX 2080 GPU partition).\n"
                f"Running deep all-atom generative diffusion (RFD4-Proteina), co-folding (Boltz-2 / AlphaFold3), "
                f"or 100-sample PLACER ensembles on an 8GB GPU will trigger fatal CUDA Out-Of-Memory (OOM) errors.\n"
                f"Remedy: Re-route '{operation_name}' to 'tamarind' (REST API) or 'modal' (Serverless A100/H100) "
                f"in pipeline_config.compute_routes."
            )

    @classmethod
    def validate_pipeline_routes(cls, routes: ComputeRouteConfig) -> Dict[str, str]:
        """
        Validates all routes in ComputeRouteConfig against VRAM guardrails.
        Raises VRAMGuardrailViolation if any invalid mapping is discovered.
        """
        mappings = {
            "stage1_rfd4": routes.stage1_rfd4,
            "stage1_rfd3": routes.stage1_rfd3,
            "stage1_ligandmpnn": routes.stage1_ligandmpnn,
            "gate1_cofolding": routes.gate1_cofolding,
            "gate2_placer": routes.gate2_placer,
            "gate3_thermo": routes.gate3_thermo,
            "gate3_esm3": routes.gate3_esm3,
            "gate4_openmm": routes.gate4_openmm,
            "gate4_fieldtools": routes.gate4_fieldtools,
            "gate4_orca": routes.gate4_orca,
            "stage3_aizymes": routes.stage3_aizymes,
        }

        for stage_name, target_backend in mappings.items():
            cls.check_route(operation_name=stage_name, backend=target_backend)

        return {k: (v.value if isinstance(v, BackendChoice) else str(v)) for k, v in mappings.items()}
