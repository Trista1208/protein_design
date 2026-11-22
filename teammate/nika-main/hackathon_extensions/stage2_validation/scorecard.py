"""Design Scorecard Module for Stage 2 Multi-Gate Validation.

Provides data structures and serialization utilities to track quantitative metrics,
per-gate boolean evaluations, failure diagnoses, and multi-objective ranking.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union


@dataclass
class DesignScorecard:
    """Consolidated scorecard tracking candidate metrics across all 4 validation gates."""

    candidate_id: str
    pdb_path: str
    sequence: str = ""
    scaffold_type: str = "de_novo"
    generation_method: str = "RFD4-Proteina"

    # Gate 1: Co-folding & Binding Affinity (Boltz-2 / AlphaFold3 / Protenix-v1)
    gate1_passed: bool = False
    catalytic_rmsd: Optional[float] = None          # Cutoff: <= 1.2 Angstroms
    attack_distance: Optional[float] = None         # Cutoff: <= 3.3 Angstroms
    affinity_probability: Optional[float] = None    # Cutoff: > 0.72 (Boltz-2 binary affinity)
    plddt: Optional[float] = None                   # Cutoff: > 70.0
    ipae: Optional[float] = None                    # Cutoff: < 10.0
    gate1_details: Dict[str, Any] = field(default_factory=dict)

    # Gate 2: Active-Site Preorganization (PLACER 100 ensemble samples)
    gate2_passed: bool = False
    ensemble_size: int = 100
    conformations_satisfying_geom: int = 0
    preorganization_percentage: Optional[float] = None  # Cutoff: >= 30.0% (30-50%)
    entropy_penalty_kcal_mol: Optional[float] = None     # -R * T * ln(P_active)
    gate2_details: Dict[str, Any] = field(default_factory=dict)

    # Gate 3: Stability & Sequence Fitness (ThermoMPNN & ESM3 / SaProt)
    gate3_passed: bool = False
    thermo_ddg: Optional[float] = None             # Cutoff: <= 1.5 kcal/mol
    esm3_fitness_score: Optional[float] = None     # Cutoff: >= -1.8 log-likelihood
    gate3_details: Dict[str, Any] = field(default_factory=dict)

    # Gate 4: Dynamics & Transition-State Energetics (OpenMM MD, FieldTools, QM/MM)
    gate4_passed: bool = False
    nac_persistence: Optional[float] = None        # Cutoff: >= 0.25 (Near-Attack Conformation fraction)
    electric_field_mv_per_cm: Optional[float] = None  # Cutoff: > 10.0 MV/cm (or > 20 MV/cm)
    dG_barrier_enzyme: Optional[float] = None      # kcal/mol
    dG_barrier_water: Optional[float] = None       # kcal/mol
    dG_barrier_delta: Optional[float] = None       # dG_enzyme - dG_water (Cutoff: < 0)
    gate4_details: Dict[str, Any] = field(default_factory=dict)

    # Global status & ranking
    all_passed: bool = False
    failed_gate: Optional[str] = None
    status: str = "PENDING"                        # PENDING, PASSED, FAILED
    composite_score: float = 0.0
    rank: Optional[int] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def evaluate_overall(self, active_gates: Optional[List[int]] = None) -> bool:
        """Update overall pass/fail status based on active gates."""
        if active_gates is None:
            active_gates = [1, 2, 3, 4]

        gate_checks = {
            1: self.gate1_passed,
            2: self.gate2_passed,
            3: self.gate3_passed,
            4: self.gate4_passed,
        }

        self.failed_gate = None
        for g in active_gates:
            if not gate_checks.get(g, False):
                self.failed_gate = f"gate{g}"
                self.all_passed = False
                self.status = "FAILED"
                return False

        self.all_passed = True
        self.status = "PASSED"
        return True

    def calculate_composite_score(self, weights: Optional[Dict[str, float]] = None) -> float:
        """Compute normalized composite score for multi-objective ranking.
        
        Higher is better.
        """
        if weights is None:
            weights = {
                "affinity": 0.25,
                "plddt": 0.15,
                "preorganization": 0.20,
                "stability": 0.15,
                "dynamics": 0.10,
                "efield": 0.15,
            }

        score = 0.0

        # Affinity (0 to 1)
        if self.affinity_probability is not None:
            score += weights.get("affinity", 0.0) * float(self.affinity_probability)

        # pLDDT (normalized 0 to 1 from 0 to 100)
        if self.plddt is not None:
            score += weights.get("plddt", 0.0) * (float(self.plddt) / 100.0)

        # Preorganization (0 to 1 from 0 to 100%)
        if self.preorganization_percentage is not None:
            score += weights.get("preorganization", 0.0) * (float(self.preorganization_percentage) / 100.0)

        # Stability: lower thermo_ddg is better (e.g. -2 to +2 kcal/mol mapped)
        if self.thermo_ddg is not None:
            # Map -2 kcal/mol -> 1.0, +2 kcal/mol -> 0.0
            norm_stab = max(0.0, min(1.0, (2.0 - float(self.thermo_ddg)) / 4.0))
            score += weights.get("stability", 0.0) * norm_stab

        # Dynamics (NAC persistence: 0 to 1)
        if self.nac_persistence is not None:
            score += weights.get("dynamics", 0.0) * float(self.nac_persistence)

        # Electric field (normalized up to 30 MV/cm)
        if self.electric_field_mv_per_cm is not None:
            norm_efield = max(0.0, min(1.0, float(self.electric_field_mv_per_cm) / 30.0))
            score += weights.get("efield", 0.0) * norm_efield

        self.composite_score = round(score, 4)
        return self.composite_score

    def to_dict(self) -> Dict[str, Any]:
        """Convert scorecard to dictionary."""
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        """Serialize scorecard to JSON string."""
        return json.dumps(self.to_dict(), indent=indent)

    def save_json(self, output_path: Union[str, Path]) -> None:
        """Write scorecard to JSON file."""
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> DesignScorecard:
        """Construct DesignScorecard from a dictionary."""
        return cls(**data)


def filter_passed(candidates: List[DesignScorecard]) -> List[DesignScorecard]:
    """Return candidates that have passed all active gates."""
    return [c for c in candidates if c.all_passed]


def rank_candidates(
    candidates: List[DesignScorecard],
    weights: Optional[Dict[str, float]] = None,
) -> List[DesignScorecard]:
    """Score and rank candidates in descending order of composite score."""
    for c in candidates:
        c.calculate_composite_score(weights=weights)

    # Sort descending by composite_score, prioritizing candidates that passed all gates
    ranked = sorted(
        candidates,
        key=lambda c: (c.all_passed, c.composite_score),
        reverse=True,
    )

    for idx, c in enumerate(ranked, start=1):
        c.rank = idx

    return ranked


def export_summary_csv(
    candidates: List[DesignScorecard],
    output_path: Union[str, Path],
) -> None:
    """Export candidates' metrics and pass/fail statuses to summary CSV."""
    out_p = Path(output_path)
    out_p.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "rank",
        "candidate_id",
        "status",
        "all_passed",
        "failed_gate",
        "composite_score",
        "gate1_passed",
        "catalytic_rmsd",
        "attack_distance",
        "affinity_probability",
        "plddt",
        "ipae",
        "gate2_passed",
        "preorganization_percentage",
        "entropy_penalty_kcal_mol",
        "gate3_passed",
        "thermo_ddg",
        "esm3_fitness_score",
        "gate4_passed",
        "nac_persistence",
        "electric_field_mv_per_cm",
        "dG_barrier_delta",
        "scaffold_type",
        "generation_method",
        "pdb_path",
    ]

    with open(out_p, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for c in candidates:
            row = {
                "rank": c.rank,
                "candidate_id": c.candidate_id,
                "status": c.status,
                "all_passed": c.all_passed,
                "failed_gate": c.failed_gate or "",
                "composite_score": c.composite_score,
                "gate1_passed": c.gate1_passed,
                "catalytic_rmsd": c.catalytic_rmsd,
                "attack_distance": c.attack_distance,
                "affinity_probability": c.affinity_probability,
                "plddt": c.plddt,
                "ipae": c.ipae,
                "gate2_passed": c.gate2_passed,
                "preorganization_percentage": c.preorganization_percentage,
                "entropy_penalty_kcal_mol": c.entropy_penalty_kcal_mol,
                "gate3_passed": c.gate3_passed,
                "thermo_ddg": c.thermo_ddg,
                "esm3_fitness_score": c.esm3_fitness_score,
                "gate4_passed": c.gate4_passed,
                "nac_persistence": c.nac_persistence,
                "electric_field_mv_per_cm": c.electric_field_mv_per_cm,
                "dG_barrier_delta": c.dG_barrier_delta,
                "scaffold_type": c.scaffold_type,
                "generation_method": c.generation_method,
                "pdb_path": c.pdb_path,
            }
            writer.writerow(row)
