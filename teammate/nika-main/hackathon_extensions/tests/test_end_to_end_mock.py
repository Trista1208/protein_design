"""Complete end-to-end integration test for the enzyme design pipeline.

Verifies the full pipeline flow in mock/dry-run mode:
Input Ingestion -> Stage 1 De Novo Design -> Stage 2 4-Gate Validation -> Stage 3 AIzymes Refinement
"""

import json
import tempfile
import unittest
from pathlib import Path

from hackathon_extensions.inputs import ingest_input
from hackathon_extensions.stage1_design.runner import Stage1Runner
from hackathon_extensions.stage2_validation.pipeline import Stage2ValidationPipeline
from hackathon_extensions.stage3_aizymes.aizymes_bridge import AIzymesBridge
from hackathon_extensions.stage3_aizymes.runner import Stage3AIzymesRunner


class TestEndToEndMockPipeline(unittest.TestCase):
    """Integration test suite executing the complete hybrid pipeline."""

    def setUp(self):
        self.mock_data_dir = Path(__file__).parent / "mock_data"
        self.mock_theozyme_pdb = self.mock_data_dir / "mock_theozyme.pdb"

    def test_full_pipeline_flow(self):
        """Execute complete flow: Ingest -> Stage 1 -> Stage 2 -> Stage 3."""
        with tempfile.TemporaryDirectory() as base_tmp:
            base_dir = Path(base_tmp)
            stage1_out = base_dir / "stage1_candidates"
            stage2_out = base_dir / "stage2_validation"
            stage3_out = base_dir / "stage3_aizymes"

            # -------------------------------------------------------------
            # Step 1: Input Ingestion (Case 1 Annotated)
            # -------------------------------------------------------------
            pipeline_input = ingest_input(
                source=str(self.mock_theozyme_pdb),
                case_type="CASE1_ANNOTATED",
                ligand_resname="LIG",
            )
            is_valid, errors = pipeline_input.validate()
            self.assertTrue(is_valid, f"Input validation failed: {errors}")

            # -------------------------------------------------------------
            # Step 2: Stage 1 Generative De Novo Design Pre-flight Dry-Run
            # -------------------------------------------------------------
            stage1_runner = Stage1Runner(default_engine="rfd4")
            s1_result = stage1_runner.run(
                pipeline_input=pipeline_input,
                output_dir=stage1_out,
                dry_run=True,
            )
            self.assertEqual(s1_result["status"], "SUCCESS")
            self.assertEqual(s1_result["candidates_count"], 0)

            manifest_path = Path(s1_result["manifest_path"])
            self.assertTrue(manifest_path.exists())

            # -------------------------------------------------------------
            # Step 3: Stage 2 4-Gate Biophysical Validation Pre-flight Dry-Run
            # -------------------------------------------------------------
            stage2_pipeline = Stage2ValidationPipeline(dry_run=True, early_stopping=True)
            scorecards = stage2_pipeline.validate_batch(
                candidate_pdbs=[str(self.mock_theozyme_pdb)],
                ligand_spec=pipeline_input.ligand,
                theozyme_spec=pipeline_input.theozyme,
                output_dir=stage2_out,
            )
            self.assertGreaterEqual(len(scorecards), 1)
            self.assertTrue((stage2_out / "scorecards.json").exists())
            self.assertTrue((stage2_out / "validation_summary.csv").exists())

            top_candidates = stage2_pipeline.get_top_candidates(scorecards, n_top=1)
            self.assertEqual(len(top_candidates), 1)
            best_scaffold_pdb = top_candidates[0].pdb_path

            # -------------------------------------------------------------
            # Step 4: Stage 3 AIzymes Directed Evolution Pre-flight Dry-Run
            # -------------------------------------------------------------
            bridge = AIzymesBridge()
            aizymes_runner = Stage3AIzymesRunner(bridge=bridge, dry_run=True)
            s3_result = aizymes_runner.run_refinement(
                scaffold_pdb=best_scaffold_pdb,
                output_dir=stage3_out,
                ligand_spec=pipeline_input.ligand,
                theozyme_spec=pipeline_input.theozyme,
                generations=2,
                designs_per_gen=5,
                dry_run=True,
            )
            self.assertEqual(s3_result["status"], "VALIDATED")
            self.assertTrue(Path(s3_result["spec_file"]).exists())


if __name__ == "__main__":
    unittest.main()
