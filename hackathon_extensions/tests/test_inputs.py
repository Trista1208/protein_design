"""Unit tests for Flexible Input Ingestion module across Cases 1, 2a, 2b, and 3."""

import unittest
from pathlib import Path

from hackathon_extensions.inputs.pipeline_input import (
    CaseType,
    LigandSpec,
    PipelineInput,
    TheozymeResidueSpec,
    TheozymeSpec,
)
from hackathon_extensions.inputs import ingest_input
from hackathon_extensions.inputs.case1_annotated import Case1AnnotatedParser
from hackathon_extensions.inputs.case2_in_between import Case2InBetweenParser
from hackathon_extensions.inputs.case3_smiles_only import Case3SmilesOnlyParser


class TestInputs(unittest.TestCase):
    """Test suite for pipeline input ingestion and data models."""

    def setUp(self):
        self.mock_data_dir = Path(__file__).parent / "mock_data"
        self.mock_theozyme_pdb = self.mock_data_dir / "mock_theozyme.pdb"
        self.mock_smi = self.mock_data_dir / "test_ester_substrate.smi"

    def test_pipeline_input_serialization_and_validation(self):
        """Verify PipelineInput dataclass serialization, deserialization, and schema validation."""
        ligand = LigandSpec(smiles="CC(=O)OC1=CC=CC=C1", resname="LIG", chain_id="Z", resnum=1)
        theozyme = TheozymeSpec(
            residues=[
                TheozymeResidueSpec(resname="CYS", chain_id="A", resnum=25, target_atoms=["SG"]),
                TheozymeResidueSpec(resname="HIS", chain_id="A", resnum=106, target_atoms=["NE2"]),
                TheozymeResidueSpec(resname="ASP", chain_id="A", resnum=125, target_atoms=["OD1"]),
            ],
            mechanism_type="cysteine_hydrolase",
        )
        inp = PipelineInput(
            case_type=CaseType.CASE1_ANNOTATED,
            ligand=ligand,
            theozyme=theozyme,
        )

        is_valid, errors = inp.validate()
        self.assertTrue(is_valid, f"Validation failed with: {errors}")

        d = inp.to_dict()
        self.assertEqual(d["case_type"], "CASE1_ANNOTATED")
        self.assertEqual(d["ligand"]["resname"], "LIG")
        self.assertEqual(len(d["theozyme"]["residues"]), 3)

        reconstructed = PipelineInput.from_dict(d)
        self.assertEqual(reconstructed.case_type, CaseType.CASE1_ANNOTATED)
        self.assertEqual(len(reconstructed.theozyme.residues), 3)

    def test_case1_annotated_ingestion(self):
        """Test Case 1 parser on mock annotated PDB with triad."""
        parser = Case1AnnotatedParser()
        res = parser.parse(
            structure_input=str(self.mock_theozyme_pdb),
            ligand_resname="LIG",
            catalytic_residues=[
                {"chain": "A", "resnum": 25, "resname": "CYS"},
                {"chain": "A", "resnum": 106, "resname": "HIS"},
                {"chain": "A", "resnum": 125, "resname": "ASP"},
            ],
        )
        self.assertEqual(res.case_type, CaseType.CASE1_ANNOTATED)
        self.assertIsNotNone(res.theozyme)
        self.assertGreaterEqual(len(res.theozyme.residues), 3)

    def test_case2a_unannotated_holo(self):
        """Test Case 2a unannotated holo active site cropping / partial diffusion."""
        parser = Case2InBetweenParser()
        res = parser.parse_case2a_holo(
            holo_structure=str(self.mock_theozyme_pdb),
            catalytic_residues=[
                {"chain": "A", "resnum": 25, "resname": "CYS"},
                {"chain": "A", "resnum": 106, "resname": "HIS"},
                {"chain": "A", "resnum": 125, "resname": "ASP"},
            ],
            ligand_resname="LIG",
            mode="strip_scaffold",
        )
        self.assertEqual(res.case_type, CaseType.CASE2A_UNANNOTATED_HOLO)
        self.assertIsNotNone(res.theozyme)

    def test_case2b_smiles_mechanism(self):
        """Test Case 2b SMILES + known catalytic mechanism triad assembly."""
        smiles = "CC(=O)Oc1ccccc1"
        parser = Case2InBetweenParser()
        res = parser.parse_case2b_smiles_mechanism(
            smiles=smiles,
            mechanism_type="cysteine_hydrolase",
        )
        self.assertEqual(res.case_type, CaseType.CASE2B_SMILES_MECHANISM)
        self.assertEqual(len(res.theozyme.residues), 3)

    def test_case3_smiles_only_theozyme_first(self):
        """Test Case 3 SMILES only input with theozyme assembly."""
        smiles = "CC1=CC(=O)OC2=C1C=CC(=C2)OC(=O)CC3=CC=CC=C3"
        parser = Case3SmilesOnlyParser()
        res = parser.parse(smiles=smiles, hypothesized_mechanism="cysteine_hydrolase")
        self.assertEqual(res.case_type, CaseType.CASE3_SMILES_ONLY)
        self.assertIsNotNone(res.theozyme)
        self.assertEqual(len(res.theozyme.residues), 3)

    def test_case3_smiles_only_binder_first(self):
        """Test Case 3 SMILES only input with binder-first shortcut."""
        smiles = "CC1=CC(=O)OC2=C1C=CC(=C2)OC(=O)CC3=CC=CC=C3"
        parser = Case3SmilesOnlyParser()
        res = parser.parse(smiles=smiles, mode="binder_first")
        self.assertEqual(res.case_type, CaseType.CASE3_SMILES_ONLY)
        self.assertIsNone(res.theozyme)
        self.assertIn("aizymes_evolution_plan", res.metadata)
        self.assertEqual(res.metadata["strategy"], "binder_first")

    def test_unified_ingest_input_factory(self):
        """Test the top-level ingest_input dispatch factory."""
        # Test Case 3 via raw SMILES
        inp = ingest_input(source="CC(=O)Oc1ccccc1")
        self.assertEqual(inp.case_type, CaseType.CASE3_SMILES_ONLY)

        # Test Case 1 via input_case override
        inp_c1 = ingest_input(
            source=str(self.mock_theozyme_pdb),
            case_type="CASE1_ANNOTATED",
            ligand_resname="LIG",
        )
        self.assertEqual(inp_c1.case_type, CaseType.CASE1_ANNOTATED)


if __name__ == "__main__":
    unittest.main()
