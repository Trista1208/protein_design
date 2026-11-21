"""Unit tests for Stage 2 Multi-Gate Validation Pipeline and Scorecards in dry-run mode."""

import json
import tempfile
import unittest
from pathlib import Path

from hackathon_extensions.stage2_validation.gate1_cofolding import Gate1CofoldingEvaluator
from hackathon_extensions.stage2_validation.gate2_preorganization import Gate2PreorganizationEvaluator
from hackathon_extensions.stage2_validation.gate3_stability import Gate3StabilityEvaluator
from hackathon_extensions.stage2_validation.gate4_dynamics_qm import Gate4DynamicsQMEvaluator
from hackathon_extensions.stage2_validation.pipeline import Stage2ValidationPipeline
from hackathon_extensions.stage2_validation.scorecard import (
    DesignScorecard,
    export_summary_csv,
    filter_passed,
    rank_candidates,
)


class TestStage2DryRun(unittest.TestCase):
    """Test suite for 4-gate biophysical validation and scorecard assessment."""

    def setUp(self):
        self.mock_data_dir = Path(__file__).parent / "mock_data"
        self.mock_pdb = self.mock_data_dir / "mock_rfd4_output.pdb"

    def test_design_scorecard_metrics_and_ranking(self):
        """Verify DesignScorecard evaluation, composite scoring, and candidate ranking."""
        c1 = DesignScorecard(
            candidate_id="design_01",
            pdb_path=str(self.mock_pdb),
            gate1_passed=True,
            affinity_probability=0.88,
            plddt=85.0,
            gate2_passed=True,
            preorganization_percentage=45.0,
            gate3_passed=True,
            thermo_ddg=-0.5,
            gate4_passed=True,
            nac_persistence=0.35,
            electric_field_mv_per_cm=18.0,
        )
        c1.evaluate_overall()
        self.assertTrue(c1.all_passed)
        self.assertEqual(c1.status, "PASSED")

        c2 = DesignScorecard(
            candidate_id="design_02",
            pdb_path=str(self.mock_pdb),
            gate1_passed=False,
            affinity_probability=0.45,
            plddt=62.0,
        )
        c2.evaluate_overall()
        self.assertFalse(c2.all_passed)
        self.assertEqual(c2.status, "FAILED")
        self.assertEqual(c2.failed_gate, "gate1")

        ranked = rank_candidates([c2, c1])
        self.assertEqual(ranked[0].candidate_id, "design_01")
        self.assertEqual(ranked[0].rank, 1)

        passed_only = filter_passed([c1, c2])
        self.assertEqual(len(passed_only), 1)
        self.assertEqual(passed_only[0].candidate_id, "design_01")

    def test_gate1_cofolding_evaluator(self):
        """Verify Gate 1 evaluator metrics, cutoffs, and ESMFold deficiency explanation."""
        evaluator = Gate1CofoldingEvaluator(dry_run=True)
        res = evaluator.evaluate(candidate_pdb=self.mock_pdb)

        self.assertTrue(res["passed"])
        self.assertEqual(res["status"], "VALIDATED")
        self.assertTrue(res["dry_run"])

        # Check ESMFold comparison analysis
        esm_note = Gate1CofoldingEvaluator.explain_esmfold_deficiency()
        self.assertIn("ESMFold", esm_note)
        self.assertIn("Boltz-2", esm_note)

    def test_gate2_preorganization_evaluator(self):
        """Verify Gate 2 PLACER 100 ensemble pre-flight validation."""
        evaluator = Gate2PreorganizationEvaluator(ensemble_size=100, dry_run=True)
        res = evaluator.evaluate(candidate_pdb=self.mock_pdb)

        self.assertTrue(res["passed"])
        self.assertEqual(res["status"], "VALIDATED")
        self.assertEqual(res["ensemble_size"], 100)

    def test_gate3_stability_evaluator(self):
        """Verify Gate 3 ThermoMPNN / ESM3 pre-flight validation."""
        evaluator = Gate3StabilityEvaluator(dry_run=True)
        res = evaluator.evaluate(candidate_pdb=self.mock_pdb)

        self.assertTrue(res["passed"])
        self.assertEqual(res["status"], "VALIDATED")
        self.assertIn("routing", res)

    def test_gate4_dynamics_qm_evaluator(self):
        """Verify Gate 4 sequential OpenMM MD, FieldTools, and ORCA barrier pre-flight validation."""
        evaluator = Gate4DynamicsQMEvaluator(dry_run=True)
        res = evaluator.evaluate(candidate_pdb=self.mock_pdb)

        self.assertTrue(res["passed"])
        self.assertEqual(res["status"], "VALIDATED")
        self.assertIn("routing", res)

    def test_stage2_pipeline_batch_and_reports(self):
        """Verify Stage2ValidationPipeline batch processing, ranking, and report output."""
        pipeline = Stage2ValidationPipeline(dry_run=True, early_stopping=True)

        with tempfile.TemporaryDirectory() as tmp_dir:
            out_dir = Path(tmp_dir)
            candidates = [self.mock_pdb, self.mock_pdb]
            ranked = pipeline.validate_batch(candidates, output_dir=out_dir)

            self.assertEqual(len(ranked), 2)
            self.assertTrue((out_dir / "scorecards.json").exists())
            self.assertTrue((out_dir / "validation_summary.csv").exists())

            with open(out_dir / "scorecards.json", "r", encoding="utf-8") as f:
                scorecards = json.load(f)
            self.assertEqual(len(scorecards), 2)

            top = pipeline.get_top_candidates(ranked, n_top=1)
            self.assertEqual(len(top), 1)


if __name__ == "__main__":
    unittest.main()
