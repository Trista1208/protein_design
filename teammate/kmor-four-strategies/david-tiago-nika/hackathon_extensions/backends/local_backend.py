"""
Local Compute Backend.
Executes lightweight local calculations (orchestration, geometric alignment,
motif extraction, and rapid local dry-run simulation).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Callable, Dict, Optional

from hackathon_extensions.backends.base import ComputeBackend, JobStatus


class LocalBackend(ComputeBackend):
    """
    Executes tasks locally on the host machine.
    Used for orchestration, fast geometric checks, and immediate mock responses.
    """

    def __init__(self, dry_run: bool = True):
        super().__init__(backend_name="local", dry_run=dry_run)
        self._jobs: Dict[str, Dict[str, Any]] = {}

    def submit_job(self, job_spec: Dict[str, Any]) -> str:
        job_id = self.generate_job_id("local")
        operation = job_spec.get("operation", "generic_local")

        # Record job
        self._jobs[job_id] = {
            "job_id": job_id,
            "operation": operation,
            "status": JobStatus.RUNNING,
            "spec": job_spec,
            "results": None,
            "error": None,
        }

        # If a callable task is passed, execute it directly
        task_fn: Optional[Callable[..., Any]] = job_spec.get("callable")
        if callable(task_fn):
            try:
                args = job_spec.get("args", ())
                kwargs = job_spec.get("kwargs", {})
                res = task_fn(*args, **kwargs)
                self._jobs[job_id]["results"] = res
                self._jobs[job_id]["status"] = JobStatus.COMPLETED
            except Exception as e:
                self._jobs[job_id]["error"] = str(e)
                self._jobs[job_id]["status"] = JobStatus.FAILED
        else:
            # Generate local execution result (dry-run or built-in operation)
            res = self._execute_local_operation(operation, job_spec)
            self._jobs[job_id]["results"] = res
            self._jobs[job_id]["status"] = JobStatus.COMPLETED

        return job_id

    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        job_data = self._jobs.get(job_id)
        if not job_data:
            return {"job_id": job_id, "status": JobStatus.FAILED, "message": "Unknown job_id"}
        return {
            "job_id": job_id,
            "backend": self.backend_name,
            "status": job_data["status"],
            "message": job_data.get("error"),
        }

    def fetch_results(self, job_id: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        job_data = self._jobs.get(job_id)
        if not job_data:
            raise KeyError(f"Job {job_id} not found in LocalBackend registry.")

        if job_data["status"] == JobStatus.FAILED:
            raise RuntimeError(f"Local job {job_id} failed: {job_data.get('error')}")

        results = job_data.get("results") or {}

        # Optionally persist outputs if output_dir provided
        if output_dir:
            out_p = Path(output_dir)
            out_p.mkdir(parents=True, exist_ok=True)
            res_summary_file = out_p / f"{job_id}_summary.json"
            try:
                import json
                res_summary_file.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
            except Exception:
                pass

        return results

    def _execute_local_operation(self, operation: str, job_spec: Dict[str, Any]) -> Dict[str, Any]:
        """Built-in handlers for local operations and mocks."""
        op_norm = operation.lower()

        if "geometry" in op_norm or "rmsd" in op_norm:
            return {
                "operation": operation,
                "rmsd_angstrom": 0.82,
                "attack_distance_angstrom": 3.12,
                "catalytic_clash_count": 0,
                "status": "PASS",
            }
        elif "ligand" in op_norm:
            return {
                "operation": operation,
                "heavy_atom_count": 18,
                "rotatable_bonds": 4,
                "formal_charge": -1,
                "status": "PARSED",
            }
        else:
            return {
                "operation": operation,
                "backend": "local",
                "success": True,
                "message": f"Local execution completed for {operation}",
            }
