"""
ETH Euler Cluster Compute Backend.
Manages job submission, Slurm script generation, and SSH execution
for 8GB RTX 2080 GPU and CPU nodes on the ETH Zurich Euler HPC cluster:
- LigandMPNN sequence redesign
- ThermoMPNN thermodynamic stability screening
- OpenMM all-atom MD & NAC persistence
- FieldTools electrostatic transition-state alignment
- ORCA QM/MM activation barrier calculations
- AIzymes in silico directed evolution
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

from hackathon_extensions.backends.base import ComputeBackend, JobStatus

logger = logging.getLogger(__name__)


class ETHClusterBackend(ComputeBackend):
    """
    Backend interface for ETH Zurich Euler cluster.
    Generates Slurm sbatch scripts, handles SSH/SFTP transfer,
    and emulates execution with realistic biophysical metrics in dry-run mode.
    """

    def __init__(
        self,
        host: str = "euler.ethz.ch",
        user: Optional[str] = None,
        port: int = 22,
        key_filename: Optional[str] = None,
        password: Optional[str] = None,
        slurm_partition: str = "gpu.2080",
        scratch_dir: str = "/cluster/scratch",
        dry_run: bool = True
    ):
        super().__init__(backend_name="eth_cluster", dry_run=dry_run)
        self.host = host
        self.user = user or os.getenv("ETH_EULER_USER")
        self.port = port
        self.key_filename = key_filename or os.getenv("ETH_EULER_KEY")
        self.password = password or os.getenv("ETH_EULER_PASSWORD")
        self.slurm_partition = slurm_partition
        self.scratch_dir = scratch_dir

        if not self.user and not self.dry_run:
            logger.warning("No ETH Euler user credentials provided. Defaulting to dry_run mode.")
            self.dry_run = True

        self._submitted_jobs: Dict[str, Dict[str, Any]] = {}

    def submit_job(self, job_spec: Dict[str, Any]) -> str:
        job_id = self.generate_job_id("euler")
        operation = job_spec.get("operation", "thermompnn")

        # Generate the Slurm sbatch script
        sbatch_script = self.generate_sbatch_script(operation, job_spec, job_id)

        self._submitted_jobs[job_id] = {
            "job_id": job_id,
            "operation": operation,
            "spec": job_spec,
            "sbatch_script": sbatch_script,
            "status": JobStatus.PENDING,
            "results": None,
        }

        if self.dry_run:
            logger.info("[Dry-Run] Pre-flight generated Slurm sbatch script for Euler task '%s'", operation)
            self._submitted_jobs[job_id]["status"] = JobStatus.COMPLETED
            self._submitted_jobs[job_id]["results"] = {
                "job_id": job_id,
                "operation": operation,
                "status": "VALIDATED",
                "dry_run": True,
                "sbatch_script": sbatch_script,
            }
            return job_id

        # Real SSH / Slurm sbatch execution on Euler via Paramiko
        try:
            import paramiko  # type: ignore

            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(
                hostname=self.host,
                port=self.port,
                username=self.user,
                key_filename=self.key_filename,
                password=self.password,
                timeout=15,
            )

            remote_script_path = f"{self.scratch_dir}/{self.user}/jobs/{job_id}.sbatch"
            sftp = ssh.open_sftp()
            try:
                sftp.mkdir(f"{self.scratch_dir}/{self.user}")
                sftp.mkdir(f"{self.scratch_dir}/{self.user}/jobs")
            except IOError:
                pass

            with sftp.file(remote_script_path, "w") as f:
                f.write(sbatch_script)
            sftp.close()

            stdin, stdout, stderr = ssh.exec_command(f"sbatch {remote_script_path}")
            out = stdout.read().decode().strip()
            slurm_job_id = out.split()[-1] if out else job_id

            self._submitted_jobs[job_id]["slurm_job_id"] = slurm_job_id
            self._submitted_jobs[job_id]["status"] = JobStatus.RUNNING
            ssh.close()
            return job_id

        except ImportError as e:
            raise RuntimeError("Paramiko is required for live ETH Euler execution. Install via 'pip install paramiko'.") from e
        except Exception as e:
            logger.error(f"SSH Euler submission error: {e}")
            raise RuntimeError(f"Failed to submit Slurm job to ETH Euler cluster: {e}") from e

    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        job = self._submitted_jobs.get(job_id)
        if not job:
            return {"job_id": job_id, "status": JobStatus.FAILED, "message": "Job not found"}

        if self.dry_run or job.get("results") is not None:
            return {"job_id": job_id, "backend": "eth_cluster", "status": JobStatus.COMPLETED, "progress": 1.0}

        slurm_id = job.get("slurm_job_id")
        if not slurm_id:
            return {"job_id": job_id, "backend": "eth_cluster", "status": job["status"]}

        # Query squeue on Euler
        try:
            import paramiko  # type: ignore
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(
                hostname=self.host,
                port=self.port,
                username=self.user,
                key_filename=self.key_filename,
                password=self.password,
                timeout=10,
            )
            stdin, stdout, stderr = ssh.exec_command(f"squeue -j {slurm_id} -h -o %T")
            out = stdout.read().decode().strip()
            ssh.close()

            if not out:  # Job completed and left queue
                job["status"] = JobStatus.COMPLETED
            elif "RUNNING" in out:
                job["status"] = JobStatus.RUNNING
            elif "PENDING" in out:
                job["status"] = JobStatus.PENDING
            else:
                job["status"] = JobStatus.COMPLETED

            return {"job_id": job_id, "backend": "eth_cluster", "status": job["status"]}
        except Exception as e:
            return {"job_id": job_id, "backend": "eth_cluster", "status": job["status"], "error": str(e)}

    def fetch_results(self, job_id: str, output_dir: Optional[str] = None) -> Dict[str, Any]:
        job = self._submitted_jobs.get(job_id)
        if not job:
            raise KeyError(f"Job {job_id} not found in ETHClusterBackend registry.")

        if self.dry_run:
            return {"job_id": job_id, "status": "VALIDATED", "dry_run": True, "sbatch_script": job.get("sbatch_script")}

        results = job.get("results")
        if results is None:
            raise RuntimeError(f"Euler Slurm job {job_id} has not completed yet or output files have not been retrieved.")

        return results

        if output_dir:
            out_p = Path(output_dir)
            out_p.mkdir(parents=True, exist_ok=True)
            # Save sbatch script and results summary
            if "sbatch_script" in job:
                (out_p / f"{job_id}.sbatch").write_text(job["sbatch_script"], encoding="utf-8")
            summary_file = out_p / f"{job_id}_euler_results.json"
            try:
                import json
                summary_file.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
            except Exception:
                pass

        return results

    def generate_sbatch_script(self, operation: str, job_spec: Dict[str, Any], job_id: str) -> str:
        """Generates realistic Slurm sbatch script for Euler GPU / CPU partition."""
        user = self.user or "user"
        partition = self.slurm_partition
        time_limit = job_spec.get("time_limit", "04:00:00")
        mem_gb = job_spec.get("mem_gb", 16)
        num_gpus = 1 if "gpu" in partition else 0

        script_lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={operation}_{job_id}",
            f"#SBATCH --output={self.scratch_dir}/{user}/logs/{job_id}.out",
            f"#SBATCH --error={self.scratch_dir}/{user}/logs/{job_id}.err",
            f"#SBATCH --time={time_limit}",
            f"#SBATCH --mem={mem_gb}G",
            f"#SBATCH --partition={partition}",
        ]

        if num_gpus > 0:
            script_lines.append(f"#SBATCH --gres=gpu:{num_gpus}")

        script_lines.extend([
            "",
            "# Load Euler Environment Modules",
            "module purge",
            "module load gcc/11.2.0 python/3.11.2 cuda/12.1",
            "",
            f"echo 'Starting {operation} job {job_id} on $(hostname) at $(date)'",
            f"WORKDIR={self.scratch_dir}/{user}/work/{job_id}",
            "mkdir -p $WORKDIR && cd $WORKDIR",
            "",
        ])

        op_norm = operation.lower()
        if "ligandmpnn" in op_norm:
            script_lines.extend([
                "# LigandMPNN sequence redesign within 8GB VRAM",
                "python /cluster/apps/ligandmpnn/run.py --model_type ligand_mpnn --temperature 0.15",
            ])
        elif "thermo" in op_norm:
            script_lines.extend([
                "# ThermoMPNN stability assessment",
                "python /cluster/apps/thermompnn/score.py --input structure.pdb --output ddg.json",
            ])
        elif "openmm" in op_norm:
            script_lines.extend([
                "# OpenMM MD relaxation & NAC trajectory analysis",
                "python -m openmm.simulate --pdb structure.pdb --steps 50000 --platform CUDA",
            ])
        elif "field" in op_norm:
            script_lines.extend([
                "# FieldTools electrostatic projection along catalytic dipole",
                "python -m fieldtools.project --pqr structure.pqr --dipole-axis C1-O1",
            ])
        elif "orca" in op_norm:
            script_lines.extend([
                "# ORCA QM/MM Transition State barrier calculation",
                "orca qmmm_input.inp > orca_output.out",
            ])
        elif "aizymes" in op_norm:
            script_lines.extend([
                "# AIzymes in silico directed evolution screening",
                "python -m aizymes.evolve --target enzyme.pdb --rounds 5 --library-size 250",
            ])
        else:
            script_lines.append(f"echo 'Running generic task {operation}'")

        script_lines.extend([
            "",
            "echo 'Job finished successfully at $(date)'",
        ])

        return "\n".join(script_lines)
