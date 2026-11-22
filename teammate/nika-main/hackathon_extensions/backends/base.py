"""
Base abstract compute backend for the Enzyme Design Pipeline.
Defines the standard execution lifecycle: job specification, submission,
status monitoring, and artifact retrieval.
"""

from __future__ import annotations

import time
import uuid
from abc import ABC, abstractmethod
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Optional


class JobStatus(str, Enum):
    """Standardized lifecycle status across all execution backends."""
    SUBMITTED = "SUBMITTED"
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"


class ComputeBackend(ABC):
    """
    Abstract base class for all compute execution backends
    (Local, Tamarind Bio, Modal AI, ETH Euler Cluster).
    """

    def __init__(self, backend_name: str, dry_run: bool = True):
        self.backend_name = backend_name
        self.dry_run = dry_run

    @abstractmethod
    def submit_job(self, job_spec: Dict[str, Any]) -> str:
        """
        Submit a compute job specification asynchronously.
        Returns unique job_id string.
        """
        pass

    @abstractmethod
    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        """
        Check the execution status of a submitted job.
        Returns a dictionary containing at least:
        - 'job_id': str
        - 'status': JobStatus (or status string)
        - 'backend': str
        - 'message': Optional[str]
        - 'progress': Optional[float]
        """
        pass

    @abstractmethod
    def fetch_results(self, job_id: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        """
        Download or retrieve the output artifacts and metrics for a completed job.
        Returns dictionary containing output files, metrics, and logs.
        """
        pass

    def run_synchronous(
        self,
        job_spec: Dict[str, Any],
        output_dir: Optional[str] = None,
        poll_interval_sec: float = 1.0,
        max_timeout_sec: float = 300.0,
    ) -> Dict[str, Any]:
        """
        Submits job and waits until completion or timeout, then fetches results.
        In dry_run mode, returns simulated mock results immediately without blocking.
        """
        job_id = self.submit_job(job_spec)

        if self.dry_run:
            # Immediate return for mock execution without waiting
            return self.fetch_results(job_id=job_id, output_dir=output_dir)

        # Polling loop for real backend execution
        start_time = time.time()
        while time.time() - start_time < max_timeout_sec:
            status_info = self.get_job_status(job_id)
            current_status = status_info.get("status")

            if current_status == JobStatus.COMPLETED or current_status == "COMPLETED":
                return self.fetch_results(job_id=job_id, output_dir=output_dir)
            elif current_status == JobStatus.FAILED or current_status == "FAILED":
                err_msg = status_info.get("message", "Job failed without explicit message")
                raise RuntimeError(f"Backend {self.backend_name} job {job_id} failed: {err_msg}")

            time.sleep(poll_interval_sec)

        raise TimeoutError(f"Backend {self.backend_name} job {job_id} timed out after {max_timeout_sec}s")

    @staticmethod
    def generate_job_id(prefix: str = "job") -> str:
        """Generate clean, unique identifier for jobs."""
        short_uuid = uuid.uuid4().hex[:8]
        return f"{prefix}_{short_uuid}"
