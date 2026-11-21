"""Gate 3: Stability and Sequence Fitness Verification.

Uses ThermoMPNN and ESM3 / SaProt to verify thermodynamic stability of folding
and evolutionary sequence plausibility.

Physical & Evolutionary Rationale:
----------------------------------
1. ThermoMPNN:
   Calculates predicted Delta Delta G of folding (kcal/mol) relative to a reference
   or absolute thermodynamic stability score. De novo designs with heavily destabilizing
   mutations (Delta Delta G > 1.5 kcal/mol) are prone to aggregation, misfolding, or low
   expression titers in wet-lab expression.
   - Resource footprint: Low VRAM (< 4 GB), highly suitable for ETH Euler RTX 2080 or local CPU/GPU.

2. ESM3 / SaProt:
   Calculates masked sequence pseudo-log-likelihood / evolutionary fitness score.
   Ensures that generated sequences respect natural co-evolutionary patterns and secondary
   structure propensities, filtering out hallucinated non-biophysical motifs.
   - Resource footprint: Medium-to-high VRAM, routed to Modal AI serverless GPU.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

try:
    from hackathon_extensions.config.config_schema import BackendChoice, GateCutoffsConfig
except ImportError:
    from dataclasses import dataclass
    @dataclass
    class GateCutoffsConfig:
        max_ddg_destabilization: float = 1.5
        min_esm3_log_likelihood: float = -2.5

logger = logging.getLogger(__name__)


class Gate3StabilityEvaluator:
    """Evaluates candidate scaffolds against Gate 3 (Stability & Sequence Fitness)."""

    def __init__(
        self,
        cutoffs: Optional[GateCutoffsConfig] = None,
        backend: Optional[Any] = None,
        dry_run: bool = True,
    ):
        self.cutoffs = cutoffs or GateCutoffsConfig()
        self.backend = backend
        self.dry_run = dry_run

    def evaluate(
        self,
        candidate_pdb: Union[str, Path],
        sequence: Optional[str] = None,
        wt_pdb: Optional[Union[str, Path]] = None,
        dry_run: Optional[bool] = None,
    ) -> Dict[str, Any]:
        """Run Gate 3 stability evaluation.

        Returns:
            Dict containing:
                - 'passed': bool
                - 'thermo_ddg': float (<= 1.5 kcal/mol)
                - 'esm3_fitness_score': float (>= -2.5 or -1.8)
                - 'checks': dict of per-metric booleans
                - 'details': dict
        """
        is_dry_run = self.dry_run if dry_run is None else dry_run
        pdb_path = Path(candidate_pdb)

        # Extract sequence from PDB if not provided
        if not sequence and pdb_path.exists():
            sequence = self._extract_sequence_from_pdb(pdb_path)

        if is_dry_run:
            logger.info("[Gate 3 Dry-Run] Pre-flight validating stability calculation specs for '%s'", pdb_path.name)
            return {
                "passed": True,
                "status": "VALIDATED",
                "dry_run": True,
                "message": f"Pre-flight validated ThermoMPNN / ESM3 stability specs for {pdb_path.name}.",
                "sequence_length": len(sequence or ""),
                "thermo_ddg": None,
                "esm3_fitness_score": None,
                "routing": {
                    "thermo_engine": "ThermoMPNN (ETH Euler RTX 2080 / Local CPU)",
                    "fitness_engine": "ESM3 / SaProt (Modal AI Serverless GPU)",
                },
            }

        return self._run_live_stability(pdb_path, sequence, wt_pdb)

    def _extract_sequence_from_pdb(self, pdb_path: Path) -> str:
        """Parse standard SEQRES or ATOM records to single-letter amino acid sequence."""
        three_to_one = {
            "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
            "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
            "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
            "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
        }
        seq: List[str] = []
        seen_res = set()
        try:
            with open(pdb_path, "r", encoding="utf-8") as f:
                for line in f:
                    if line.startswith("ATOM") and line[12:16].strip() == "CA":
                        res_name = line[17:20].strip()
                        chain = line[21:22].strip()
                        res_num = line[22:26].strip()
                        key = (chain, res_num)
                        if key not in seen_res and res_name in three_to_one:
                            seen_res.add(key)
                            seq.append(three_to_one[res_name])
        except Exception:
            return ""
        return "".join(seq)

    def _run_live_stability(
        self,
        pdb_path: Path,
        sequence: Optional[str],
        wt_pdb: Optional[Union[str, Path]],
    ) -> Dict[str, Any]:
        """Execute live ThermoMPNN on ETH Euler and ESM3 on Modal AI."""
        logger.info("Executing live stability assessment for %s", pdb_path.name)
        if self.backend and hasattr(self.backend, "submit_stability"):
            return self.backend.submit_stability(
                pdb_path=str(pdb_path),
                sequence=sequence,
                wt_pdb=str(wt_pdb) if wt_pdb else None,
            )
        raise RuntimeError(f"No stability backend available to execute live stability assessment for {pdb_path.name}")
