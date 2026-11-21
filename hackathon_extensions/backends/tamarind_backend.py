"""Tamarind Bio Compute Backend.

Client for Tamarind Bio API:
- Pre-flight schema validation via POST /api/validate-job (at $0 cost)
- RFD4-Proteina generative backbone diffusion
- RFD3 structural scaffolding
- Boltz-2 / AlphaFold3 / Chai-1 all-atom co-folding and ligand affinity
"""

from __future__ import annotations

import json
import logging
import os
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from hackathon_extensions.backends.base import ComputeBackend, JobStatus

logger = logging.getLogger(__name__)


class TamarindBackend(ComputeBackend):
    """Client for Tamarind Bio REST API.

    Implements:
    - POST /api/validate-job: Pre-flight check at $0 cost before spending allowances.
    - POST /api/submit-job: Live submission consuming job credits.
    - GET /api/jobs/{id}: Status polling and result retrieval.
    """

    TAMARIND_BASE_URL = "https://app.tamarind.bio/api"

    def __init__(self, api_key: Optional[str] = None, dry_run: bool = False, require_approval: bool = True):
        super().__init__(backend_name="tamarind", dry_run=dry_run)
        self.api_key = api_key or os.getenv("TAMARIND_API_KEY", "04b8a253-a5cc-4ea6-9bfa-5ffc7c60d52a")
        self.require_approval = require_approval
        if not self.dry_run and not self.api_key:
            raise ValueError(
                "TAMARIND_API_KEY must be provided via environment variable, user profile, or constructor for live execution."
            )
        self._submitted_jobs: Dict[str, Dict[str, Any]] = {}

    def _list_remote_files(self) -> List[str]:
        """Return file names in the Tamarind account via GET /api/files (read-only)."""
        req = urllib.request.Request(f"{self.TAMARIND_BASE_URL}/files", headers={"x-api-key": self.api_key})
        try:
            with urllib.request.urlopen(req, timeout=20) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            files = data if isinstance(data, list) else data.get("files", [])
            return [str(f.get("name", f)) if isinstance(f, dict) else str(f) for f in files]
        except Exception as e:
            logger.warning("Could not list Tamarind files: %s", e)
            return []

    def upload_file(self, file_path: Union[str, Path]) -> bool:
        """Upload structure file to Tamarind Bio via PUT /api/upload/{filename}."""
        p = Path(file_path)
        if not p.exists():
            return False

        filename = p.name
        url = f"{self.TAMARIND_BASE_URL}/upload/{filename}"
        content = p.read_bytes()
        headers = {
            "x-api-key": self.api_key,
            "Content-Type": "application/octet-stream",
        }

        # 1. Prefer requests if installed (natively and safely handles 307/308 redirects)
        try:
            import requests
            resp = requests.put(url, headers=headers, data=content, timeout=30, allow_redirects=True)
            if resp.status_code in (200, 201):
                logger.info("Successfully uploaded %s to Tamarind Bio (HTTP %d).", filename, resp.status_code)
                return True
            logger.warning("Tamarind upload via requests returned HTTP %d: %s", resp.status_code, resp.text)
        except ImportError:
            pass
        except Exception as e:
            logger.warning("Requests upload error: %s", e)

        # 2. Standard library urllib with manual 307/308 redirect handling
        req = urllib.request.Request(url, data=content, headers=headers, method="PUT")
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                logger.info("Successfully uploaded %s to Tamarind Bio (HTTP %s).", filename, resp.status)
                return resp.status in (200, 201)
        except urllib.error.HTTPError as e:
            if e.code in (301, 302, 307, 308) and "Location" in e.headers:
                redir_url = e.headers["Location"]
                if "tamarind.bio" not in redir_url:
                    # Redirect to the CDN read URL: Tamarind already holds this file.
                    # A PUT to the CDN would be unauthorized, so confirm via GET /api/files.
                    if filename in self._list_remote_files():
                        logger.info("%s already present in Tamarind account (upload redirected to CDN); skipping re-upload.", filename)
                        return True
                    logger.warning("Upload of %s redirected to %s but file is not listed in the account.", filename, redir_url)
                    return False
                logger.info("Following Tamarind upload HTTP %d redirect to %s", e.code, redir_url)
                redir_headers = {"Content-Type": "application/octet-stream"}
                if "tamarind.bio" in redir_url:
                    redir_headers["x-api-key"] = self.api_key
                req2 = urllib.request.Request(redir_url, data=content, headers=redir_headers, method="PUT")
                try:
                    with urllib.request.urlopen(req2, timeout=30) as resp2:
                        logger.info("Successfully uploaded %s after redirect (HTTP %s).", filename, resp2.status)
                        return resp2.status in (200, 201)
                except Exception as e2:
                    logger.error("Failed following redirect to %s: %s", redir_url, e2)
                    return False

            err_msg = e.read().decode("utf-8", errors="ignore")
            logger.warning("Tamarind upload HTTP %d: %s", e.code, err_msg)
            return False
        except Exception as e:
            logger.warning("Notice on uploading %s to Tamarind: %s", filename, e)
            return False

    def validate_job(
        self,
        tool_type: str,
        settings: Dict[str, Any],
        job_name: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Pre-flight check via Tamarind POST /api/validate-job at $0 cost.

        Reports whether the payload is valid, fills in normalized defaults,
        and identifies missing or incorrect fields before submitting.
        """
        if not self.api_key:
            return {"valid": False, "error": "TAMARIND_API_KEY is not set. Live validation requires a valid API key."}

        # Auto-upload pdbFile if present and found locally
        pdb_file = settings.get("pdbFile") or settings.get("inputStructure")
        if pdb_file:
            candidates = [
                Path(pdb_file),
                Path("test_pipeline") / pdb_file,
                Path(__file__).resolve().parent.parent.parent / "test_pipeline" / pdb_file,
                Path(__file__).resolve().parent.parent.parent / pdb_file,
            ]
            for c in candidates:
                if c.exists() and c.is_file():
                    logger.info("Auto-uploading structure '%s' to Tamarind Bio prior to validation...", c.name)
                    self.upload_file(c)
                    break

        url = f"{self.TAMARIND_BASE_URL}/validate-job"
        headers = {
            "x-api-key": self.api_key,
            "Content-Type": "application/json",
        }
        params: Dict[str, Any] = {
            "type": tool_type,
            "settings": settings,
        }
        if job_name:
            params["jobName"] = job_name

        payload = json.dumps(params).encode("utf-8")
        req = urllib.request.Request(url, data=payload, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                logger.info("Tamarind /api/validate-job response: valid=%s", data.get("valid"))
                return data
        except urllib.error.HTTPError as e:
            err_msg = e.read().decode("utf-8", errors="ignore")
            logger.error("Tamarind validate-job HTTP %d: %s", e.code, err_msg)
            return {"valid": False, "error": f"HTTP {e.code}: {err_msg}"}
        except Exception as e:
            logger.error("Tamarind validate-job connection error: %s", e)
            return {"valid": False, "error": str(e)}

    def submit_job(self, job_spec: Dict[str, Any]) -> str:
        job_id = self.generate_job_id("tam")

        # Resolve payload from adapter's tamarind_api_payload
        tam_payload = job_spec.get("tamarind_api_payload", {})
        operation = tam_payload.get("type", job_spec.get("operation", "rfd4-proteina"))
        settings = tam_payload.get("settings", job_spec.get("settings", job_spec.get("parameters", {})))
        job_name = job_spec.get("job_name", job_id)

        # 1. Always run POST /api/validate-job first (costs 0 credits)
        logger.info("Running pre-flight check via Tamarind POST /api/validate-job for '%s'...", operation)
        val_result = self.validate_job(tool_type=operation, settings=settings, job_name=job_name)

        if not val_result.get("valid", False) and not val_result.get("offline_check", False):
            err = val_result.get("error", "Unknown validation failure")
            missing = val_result.get("missing_fields", [])
            raise ValueError(f"Tamarind pre-flight validation failed: {err}. Missing fields: {missing}")

        final_settings = dict(val_result.get("normalized", settings))
        # Sanitize known conflicting Tamarind normalized default fields
        if final_settings.get("verifySequences") in (None, "verify-none", ""):
            final_settings.pop("verifyTopK", None)
            final_settings.pop("verifySequences", None)
            final_settings.pop("verify", None)

        # 2. Dry-run: stop after validation without submitting
        if self.dry_run:
            logger.info("[Dry-Run] Payload verified with Tamarind at $0 cost. Job not submitted.")
            self._submitted_jobs[job_id] = {
                "job_id": job_id,
                "operation": operation,
                "status": JobStatus.COMPLETED,
                "dry_run": True,
                "validated": True,
                "target_endpoint": f"{self.TAMARIND_BASE_URL}/submit-job",
                "payload": {
                    "jobName": job_name,
                    "type": operation,
                    "settings": final_settings,
                },
            }
            return job_id

        # 3. Export exact payload to pending_submissions directory and require explicit user approval
        pending_dir = Path("outputs/pending_submissions")
        pending_dir.mkdir(parents=True, exist_ok=True)
        pending_file = pending_dir / f"{job_name}_payload.json"

        submission_dict = {
            "jobName": job_name,
            "type": operation,
            "settings": final_settings,
        }
        with open(pending_file, "w", encoding="utf-8") as f:
            json.dump(submission_dict, f, indent=2)

        print("\n" + "=" * 75)
        print(" [ACTION REQUIRED: TAMARIND SUBMISSION APPROVAL]")
        print(f" Target Endpoint : {self.TAMARIND_BASE_URL}/submit-job")
        print(f" Tool Type       : {operation}")
        print(f" Job Name        : {job_name}")
        print(f" Exact File Path : {pending_file.resolve()}")
        print("-" * 75)
        print("EXACT PAYLOAD TO BE SUBMITTED:")
        print(json.dumps(submission_dict, indent=2))
        print("=" * 75)

        if self.require_approval:
            try:
                user_confirm = input("\nDo you approve submitting this EXACT payload to Tamarind Bio? Type 'yes' or 'y' to proceed: ").strip().lower()
            except (EOFError, OSError):
                user_confirm = "no"

            if user_confirm not in ("yes", "y"):
                print("\n[SUBMISSION ABORTED] User did not approve submission. Zero credits spent.\n")
                logger.warning("Job submission aborted by user. Exact spec saved at %s", pending_file)
                import sys
                sys.exit(0)

        # 4. Live submission to Tamarind Bio API
        url = f"{self.TAMARIND_BASE_URL}/submit-job"
        headers = {
            "x-api-key": self.api_key,
            "Content-Type": "application/json",
        }
        payload = json.dumps(submission_dict).encode("utf-8")

        req = urllib.request.Request(url, data=payload, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                raw_text = resp.read().decode("utf-8").strip()
                data = {}
                if raw_text:
                    try:
                        data = json.loads(raw_text)
                    except Exception:
                        data = {"id": raw_text}

                raw_id = (
                    data.get("id")
                    or data.get("jobId")
                    or data.get("job_id")
                    or resp.headers.get("x-job-id")
                    or (resp.headers.get("Location", "").rstrip("/").split("/")[-1] if resp.headers.get("Location") else job_name)
                )
                remote_id = str(raw_id).split()[0].strip("'\"")
                self._submitted_jobs[job_id] = {
                    "job_id": job_id,
                    "remote_id": remote_id,
                    "operation": operation,
                    "status": JobStatus.RUNNING,
                    "results": None,
                }
                logger.info("Tamarind job submitted successfully. Remote ID: %s", remote_id)
                print(f"\n=======================================================")
                print(f" [TAMARIND JOB SUBMITTED SUCCESSFULLY!]")
                print(f" Remote Job ID: {remote_id}")
                print(f" Track at: https://app.tamarind.bio/app/results")
                print(f"=======================================================\n")
                return job_id
        except urllib.error.HTTPError as e:
            err_msg = e.read().decode("utf-8", errors="ignore")
            logger.error("Tamarind API HTTP %d error: %s", e.code, err_msg)
            raise RuntimeError(f"Tamarind API submission failed (HTTP {e.code}): {err_msg}") from e
        except Exception as e:
            logger.error("Tamarind API submission error: %s", e)
            raise RuntimeError(f"Failed to submit job to Tamarind API: {e}") from e

    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        job = self._submitted_jobs.get(job_id)
        if not job:
            return {"job_id": job_id, "status": JobStatus.FAILED, "message": "Job not found"}

        if job.get("dry_run"):
            return {"job_id": job_id, "backend": "tamarind", "status": JobStatus.COMPLETED, "progress": 1.0}

        remote_id = str(job.get("remote_id", "")).split()[0].strip("'\"")
        url = f"{self.TAMARIND_BASE_URL}/jobs/{remote_id}"
        headers = {"x-api-key": self.api_key}
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            remote_status = data.get("status", "RUNNING").upper()
            if remote_status == "COMPLETED":
                job["status"] = JobStatus.COMPLETED
                job["results"] = data.get("results", {})
            elif remote_status in ("FAILED", "ERROR"):
                job["status"] = JobStatus.FAILED
                job["message"] = data.get("error_message", "Remote job failed")
            else:
                job["status"] = JobStatus.RUNNING
            return {"job_id": job_id, "backend": "tamarind", "status": job["status"]}

    def fetch_results(self, job_id: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        job = self._submitted_jobs.get(job_id)
        if not job:
            raise KeyError(f"Job {job_id} not found in TamarindBackend registry.")

        if job.get("dry_run"):
            return {"job_id": job_id, "status": "VALIDATED", "dry_run": True, "payload": job.get("payload")}

        results = job.get("results")
        if results is None:
            raise RuntimeError(f"Job {job_id} has not completed or results are not available yet.")

        if output_dir:
            out_p = Path(output_dir)
            out_p.mkdir(parents=True, exist_ok=True)
            if "pdb_content" in results:
                pdb_file = out_p / f"{job_id}_structure.pdb"
                pdb_file.write_text(results["pdb_content"], encoding="utf-8")
                results["pdb_path"] = str(pdb_file)

        return results
