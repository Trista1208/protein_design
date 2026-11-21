"""
Backends package for the Enzyme Design Pipeline.
Exports compute execution engines, VRAM safety guardrails, and central router.
"""

from __future__ import annotations

from hackathon_extensions.backends.base import ComputeBackend, JobStatus
from hackathon_extensions.backends.eth_cluster_backend import ETHClusterBackend
from hackathon_extensions.backends.local_backend import LocalBackend
from hackathon_extensions.backends.modal_backend import ModalBackend
from hackathon_extensions.backends.router import ComputeBackendRouter
from hackathon_extensions.backends.tamarind_backend import TamarindBackend
from hackathon_extensions.backends.vram_guardrails import (
    VRAMGuardrailValidator,
    VRAMGuardrailViolation,
)

__all__ = [
    "JobStatus",
    "ComputeBackend",
    "LocalBackend",
    "TamarindBackend",
    "ModalBackend",
    "ETHClusterBackend",
    "VRAMGuardrailValidator",
    "VRAMGuardrailViolation",
    "ComputeBackendRouter",
]
