"""Stage 2 Multi-Gate Validation Pipeline Orchestrator.

Orchestrates sequential biophysical screening across:
- Gate 1: Co-folding & Binding Affinity (Boltz-2 / AlphaFold3)
- Gate 2: Active-Site Preorganization (PLACER 100 ensemble)
- Gate 3: Thermodynamic Stability & Sequence Fitness (ThermoMPNN & ESM3)
- Gate 4: Catalytic Dynamics, FieldTools & QM/MM Barriers (OpenMM, ORCA)

Features:
- Configurable active gates (toggle 1, 2, 3, 4 individually).
- Early-stopping on failure (skips expensive downstream gates to conserve compute).
- Multi-objective composite scoring and automated ranking.
- Generates consolidated `scorecards.json` and `validation_summary.csv`.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from hackathon_extensions.stage2_validation.gate1_cofolding import Gate1CofoldingEvaluator
from hackathon_extensions.stage2_validation.gate2_preorganization import Gate2PreorganizationEvaluator
from hackathon_extensions.stage2_validation.gate3_stability import Gate3StabilityEvaluator
from hackathon_extensions.stage2_validation.gate4_dynamics_qm import Gate4DynamicsQMEvaluator
from hackathon_extensions.stage2_validation.scorecard import (
    DesignScorecard,
    export_summary_csv,
    filter_passed,
    rank_candidates,
)

try:
    from hackathon_extensions.config.config_schema import (
        ComputeRouteConfig,
        GateCutoffsConfig,
        PipelineConfig,
    )
except ImportError:
    @dataclass
    class GateCutoffsConfig:
        pass

    @dataclass
    class PipelineConfig:
        gate_cutoffs: Any = field(default_factory=GateCutoffsConfig)
        dry_run: bool = True

try:
    from hackathon_extensions.inputs.pipeline_input import LigandSpec, TheozymeSpec
except ImportError:
    LigandSpec = Any
    TheozymeSpec = Any

logger = logging.getLogger(__name__)


class Stage2ValidationPipeline:
    """Master controller orchestrating Stage 2 verification gates."""

    def __init__(
        self,
        config: Optional[PipelineConfig] = None,
        active_gates: Optional[List[int]] = None,
        early_stopping: bool = True,
        dry_run: Optional[bool] = None,
        backend_router: Optional[Any] = None,
    ):
        self.config = config or PipelineConfig()
        self.active_gates = active_gates if active_gates is not None else [1, 2, 3, 4]
        self.early_stopping = early_stopping
        self.dry_run = getattr(self.config, "dry_run", True) if dry_run is None else dry_run
        self.backend_router = backend_router

        cutoffs = getattr(self.config, "gate_cutoffs", GateCutoffsConfig())

        # Initialize gate evaluators
        self.gate1 = Gate1CofoldingEvaluator(cutoffs=cutoffs, dry_run=self.dry_run)
        self.gate2 = Gate2PreorganizationEvaluator(cutoffs=cutoffs, dry_run=self.dry_run)
        self.gate3 = Gate3StabilityEvaluator(cutoffs=cutoffs, dry_run=self.dry_run)
        self.gate4 = Gate4DynamicsQMEvaluator(cutoffs=cutoffs, dry_run=self.dry_run)

    def validate_candidate(
        self,
        candidate_pdb: Union[str, Path],
        candidate_id: Optional[str] = None,
        ligand_spec: Optional[LigandSpec] = None,
        theozyme_spec: Optional[TheozymeSpec] = None,
        sequence: Optional[str] = None,
        scaffold_type: str = "de_novo",
        generation_method: str = "RFD4-Proteina",
    ) -> DesignScorecard:
        """Run configured validation gates on a single design candidate."""
        pdb_p = Path(candidate_pdb)
        cid = candidate_id or pdb_p.stem

        card = DesignScorecard(
            candidate_id=cid,
            pdb_path=str(pdb_p),
            sequence=sequence or "",
            scaffold_type=scaffold_type,
            generation_method=generation_method,
        )

        logger.info("Evaluating candidate %s (active gates: %s)", cid, self.active_gates)

        # Gate 1: Co-folding & Binding Affinity
        if 1 in self.active_gates:
            res1 = self.gate1.evaluate(
                candidate_pdb=pdb_p,
                ligand_spec=ligand_spec,
                theozyme_spec=theozyme_spec,
                dry_run=self.dry_run,
            )
            card.gate1_passed = res1["passed"]
            card.catalytic_rmsd = res1.get("catalytic_rmsd")
            card.attack_distance = res1.get("attack_distance")
            card.affinity_probability = res1.get("affinity_probability")
            card.plddt = res1.get("plddt")
            card.ipae = res1.get("ipae")
            card.gate1_details = res1

            if not card.gate1_passed and self.early_stopping:
                logger.info("Candidate %s failed Gate 1. Early-stopping active gates.", cid)
                card.evaluate_overall(self.active_gates)
                card.calculate_composite_score()
                return card

        # Gate 2: Active-Site Preorganization
        if 2 in self.active_gates:
            res2 = self.gate2.evaluate(
                candidate_pdb=pdb_p,
                ligand_spec=ligand_spec,
                theozyme_spec=theozyme_spec,
                dry_run=self.dry_run,
            )
            card.gate2_passed = res2["passed"]
            card.ensemble_size = res2.get("ensemble_size", 100)
            card.conformations_satisfying_geom = res2.get("conformations_satisfying_geom", 0)
            card.preorganization_percentage = res2.get("preorganization_percentage")
            card.entropy_penalty_kcal_mol = res2.get("entropy_penalty_kcal_mol")
            card.gate2_details = res2

            if not card.gate2_passed and self.early_stopping:
                logger.info("Candidate %s failed Gate 2. Early-stopping active gates.", cid)
                card.evaluate_overall(self.active_gates)
                card.calculate_composite_score()
                return card

        # Gate 3: Stability & Sequence Fitness
        if 3 in self.active_gates:
            res3 = self.gate3.evaluate(
                candidate_pdb=pdb_p,
                sequence=sequence,
                dry_run=self.dry_run,
            )
            card.gate3_passed = res3["passed"]
            card.thermo_ddg = res3.get("thermo_ddg")
            card.esm3_fitness_score = res3.get("esm3_fitness_score")
            card.gate3_details = res3

            if not card.gate3_passed and self.early_stopping:
                logger.info("Candidate %s failed Gate 3. Early-stopping active gates.", cid)
                card.evaluate_overall(self.active_gates)
                card.calculate_composite_score()
                return card

        # Gate 4: Dynamics & QM/MM Energetics
        if 4 in self.active_gates:
            res4 = self.gate4.evaluate(
                candidate_pdb=pdb_p,
                ligand_spec=ligand_spec,
                theozyme_spec=theozyme_spec,
                dry_run=self.dry_run,
            )
            card.gate4_passed = res4["passed"]
            card.nac_persistence = res4.get("nac_persistence")
            card.electric_field_mv_per_cm = res4.get("electric_field_mv_per_cm")
            card.dG_barrier_enzyme = res4.get("dG_barrier_enzyme")
            card.dG_barrier_water = res4.get("dG_barrier_water")
            card.dG_barrier_delta = res4.get("dG_barrier_delta")
            card.gate4_details = res4

        # Final score compilation
        card.evaluate_overall(self.active_gates)
        card.calculate_composite_score()
        return card

    def validate_batch(
        self,
        candidate_pdbs: List[Union[str, Path]],
        ligand_spec: Optional[LigandSpec] = None,
        theozyme_spec: Optional[TheozymeSpec] = None,
        output_dir: Optional[Union[str, Path]] = None,
        ranking_weights: Optional[Dict[str, float]] = None,
    ) -> List[DesignScorecard]:
        """Validate an entire cohort of designs, rank them, and write summary reports."""
        cards: List[DesignScorecard] = []
        for pdb_p in candidate_pdbs:
            card = self.validate_candidate(
                candidate_pdb=pdb_p,
                ligand_spec=ligand_spec,
                theozyme_spec=theozyme_spec,
            )
            cards.append(card)

        # Rank candidates
        ranked_cards = rank_candidates(cards, weights=ranking_weights)

        # Output persistence
        if output_dir:
            out_p = Path(output_dir)
            out_p.mkdir(parents=True, exist_ok=True)

            # JSON export
            json_file = out_p / "scorecards.json"
            with open(json_file, "w", encoding="utf-8") as f:
                json.dump([c.to_dict() for c in ranked_cards], f, indent=2)

            # Summary CSV export
            csv_file = out_p / "validation_summary.csv"
            export_summary_csv(ranked_cards, csv_file)

            logger.info("Stage 2 validation reports exported to %s", out_p)

        return ranked_cards

    def get_top_candidates(
        self,
        ranked_cards: List[DesignScorecard],
        n_top: int = 5,
        only_passed: bool = True,
    ) -> List[DesignScorecard]:
        """Filter and select top candidates for Stage 3 AIzymes refinement handoff."""
        candidates = filter_passed(ranked_cards) if only_passed else ranked_cards
        if not candidates and only_passed:
            logger.warning("No candidates passed all gates! Falling back to highest composite scores.")
            candidates = ranked_cards
        return candidates[:n_top]
