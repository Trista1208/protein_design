#!/usr/bin/env python3
"""Unified Command-Line Interface for the Enzyme Design Pipeline.

Supports standalone stage execution and end-to-end multi-stage workflows:
- `stage1`: De novo generative design (RFD4-Proteina / RFD3 + LigandMPNN)
- `stage2`: Standalone 4-gate biophysical validation on candidate PDBs
- `stage1-2`: De novo generative design + 4-gate validation triage
- `full`: Complete hybrid pipeline: Stage 1 -> Stage 2 -> Stage 3 AIzymes refinement
- `aizymes-standalone`: Direct AIzymes evolutionary optimization on existing holo-complex

Examples:
    # Full hybrid pipeline in dry-run mode from SMILES:
    python hackathon_extensions/cli.py full --input "CC(=O)Oc1ccccc1" --input-case case3 --dry-run

    # Standalone Stage 1 de novo generation with RFD4:
    python hackathon_extensions/cli.py stage1 --input "hackathon_extensions/tests/mock_data/mock_theozyme.pdb" --input-case case1 --dry-run

    # Standalone Stage 2 validation on existing PDB candidates:
    python hackathon_extensions/cli.py stage2 --candidates "outputs/stage1_candidates/*.pdb" --output-dir "outputs/stage2_validation" --dry-run
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

# Ensure repository root is on sys.path
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from hackathon_extensions.backends.router import ComputeBackendRouter
from hackathon_extensions.config import load_pipeline_config
from hackathon_extensions.config.config_schema import BackendChoice, DesignMode, PipelineConfig, ThinkingBudget
from hackathon_extensions.config.credentials_manager import CredentialsManager
from hackathon_extensions.inputs import ingest_input
from hackathon_extensions.inputs.pipeline_input import CaseType, PipelineInput
from hackathon_extensions.stage1_design.runner import Stage1Runner
from hackathon_extensions.stage2_validation.pipeline import Stage2ValidationPipeline
from hackathon_extensions.stage3_aizymes.aizymes_bridge import AIzymesBridge
from hackathon_extensions.stage3_aizymes.runner import Stage3AIzymesRunner

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("EnzymeDesignCLI")


def print_banner():
    banner = r"""
========================================================================================
   ____   ____     ___  _   _ _____   __  __ _____     ____  _____ ____ ___ ____ _   _ 
  / __ \ / ___|   / _ \| \ | |_   _| |  \/  |  _ \   |  _ \| ____/ ___|_ _/ ___| \ | |
 / / _` | |  _   | | | |  \| | | |   | |\/| | |_) |  | | | |  _| \___ \| | |  _|  \| |
| | (_| | |_| |  | |_| | |\  | | |   | |  | |  __/   | |_| | |___ ___) | | |_| | |\  |
 \ \__,_|\____|   \___/|_| \_| |_|   |_|  |_|_|      |____/|_____|____/___|\____|_| \_|
  \____/                                                                               
   Hybrid RFD4-Proteina / RFD3  ==>  4-Gate Biophysical Validation  ==>  AIzymes Refinement
========================================================================================
"""
    print(banner)


def resolve_case_type(case_str: Optional[str]) -> Optional[str]:
    """Map user-friendly case strings to CaseType string."""
    if not case_str:
        return None
    c = case_str.lower().strip()
    if c in ("case1", "1", "annotated", "case1_annotated"):
        return "CASE1_ANNOTATED"
    elif c in ("case2a", "2a", "holo", "case2a_unannotated_holo"):
        return "CASE2A_UNANNOTATED_HOLO"
    elif c in ("case2b", "2b", "mechanism", "case2b_smiles_mechanism"):
        return "CASE2B_SMILES_MECHANISM"
    elif c in ("case3", "3", "smiles", "case3_smiles_only"):
        return "CASE3_SMILES_ONLY"
    return case_str


def parse_args():
    parser = argparse.ArgumentParser(
        prog="enzyme-design",
        description="Unified CLI for Hybrid De Novo Enzyme Design, Validation, and Evolutionary Refinement.",
    )

    # Subcommands
    subparsers = parser.add_subparsers(dest="command", required=True, help="Subcommand to execute")

    # Common options parent parser
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--input", "-i", type=str, default=None, help="Input specification (PDB/CIF path, JSON, or SMILES)")
    common.add_argument("--input-case", type=str, default=None, help="Input case: case1, case2a, case2b, case3")
    common.add_argument("--user", "-u", type=str, default=None, help="User credential profile name (e.g. 'tiago')")
    common.add_argument("--dry-run", action="store_true", default=False, help="Pre-flight validation check at $0 cost without submitting live GPU jobs")
    common.add_argument("--yes", "-y", action="store_true", default=False, help="Skip interactive confirmation prompts for remote job submissions")
    common.add_argument("--config", "-c", type=str, default=None, help="Custom pipeline configuration YAML/JSON")
    common.add_argument("--output-dir", "-o", type=str, default="outputs", help="Root directory for generated artifacts")

    # --- stage1 ---
    p_s1 = subparsers.add_parser("stage1", parents=[common], help="Run Stage 1 de novo generative design (RFD4-Proteina / RFD3)")
    p_s1.add_argument("--engine", "-e", type=str, default="rfd4", choices=["rfd4", "rfd3"], help="Generative backbone engine")
    p_s1.add_argument("--num-designs", "-n", type=int, default=1, help="Number of de novo designs to generate")
    p_s1.add_argument("--thinking-budget", type=str, default="medium", choices=["low", "medium", "high"], help="RFD4 latent thinking budget")
    p_s1.add_argument("--design-mode", type=str, default="enzyme_design", choices=["enzyme_design", "small_molecule_binder", "partial_diffusion"], help="Design mode")

    # --- stage2 ---
    p_s2 = subparsers.add_parser("stage2", parents=[common], help="Run Stage 2 standalone 4-gate validation on candidate PDBs")
    p_s2.add_argument("--candidates", type=str, required=True, help="Path or glob pattern for candidate PDB files to validate")
    p_s2.add_argument("--gates", type=str, default="1,2,3,4", help="Comma-separated list of active validation gates (e.g. '1,2,3,4')")
    p_s2.add_argument("--no-early-stopping", action="store_true", help="Continue evaluating all gates even after candidate fails a gate")

    # --- stage1-2 ---
    p_s12 = subparsers.add_parser("stage1-2", parents=[common], help="Run Stage 1 de novo design followed immediately by Stage 2 validation triage")
    p_s12.add_argument("--engine", "-e", type=str, default="rfd4", choices=["rfd4", "rfd3"], help="Generative backbone engine")
    p_s12.add_argument("--num-designs", "-n", type=int, default=1, help="Number of designs to generate")
    p_s12.add_argument("--gates", type=str, default="1,2,3,4", help="Comma-separated list of active validation gates")
    p_s12.add_argument("--no-early-stopping", action="store_true", help="Disable early stopping")

    # --- full ---
    p_full = subparsers.add_parser("full", parents=[common], help="Run complete hybrid pipeline: Stage 1 -> Stage 2 -> Stage 3 AIzymes refinement")
    p_full.add_argument("--engine", "-e", type=str, default="rfd4", choices=["rfd4", "rfd3"], help="Generative backbone engine")
    p_full.add_argument("--num-designs", "-n", type=int, default=1, help="Number of initial designs to generate")
    p_full.add_argument("--generations", "-g", type=int, default=3, help="Number of AIzymes evolutionary generations")
    p_full.add_argument("--designs-per-gen", type=int, default=15, help="Designs sampled per evolutionary generation")
    p_full.add_argument("--gates", type=str, default="1,2,3,4", help="Active validation gates")
    p_full.add_argument("--no-early-stopping", action="store_true", help="Disable early stopping")

    # --- aizymes-standalone ---
    p_aiz = subparsers.add_parser("aizymes-standalone", parents=[common], help="Run native AIzymes evolutionary loop on an existing complex")
    p_aiz.add_argument("--structure", "-s", type=str, required=True, help="Path to starting WT PDB structure")
    p_aiz.add_argument("--ligand-name", type=str, default="LIG", help="Ligand 3-letter residue name")
    p_aiz.add_argument("--generations", "-g", type=int, default=3, help="Number of evolutionary generations")
    p_aiz.add_argument("--designs-per-gen", type=int, default=15, help="Designs per generation")

    return parser.parse_args()


def run_stage1_cmd(args: argparse.Namespace, cfg: PipelineConfig) -> Dict[str, Any]:
    """Execute Stage 1."""
    if not args.input:
        logger.error("--input is required for stage1!")
        sys.exit(1)

    case_type = resolve_case_type(args.input_case)
    pipeline_input = ingest_input(source=args.input, case_type=case_type)

    out_dir = Path(args.output_dir) / "stage1_candidates"
    router = ComputeBackendRouter(config=cfg)
    runner = Stage1Runner(default_engine=args.engine)
    res = runner.run(
        pipeline_input=pipeline_input,
        config=cfg,
        backend_router=router,
        engine_override=args.engine,
        output_dir=out_dir,
        dry_run=args.dry_run,
    )
    logger.info("Stage 1 complete. Generated %d candidates in %s", res["candidates_count"], out_dir)
    return res


def run_stage2_cmd(args: argparse.Namespace, cfg: PipelineConfig, candidate_pdbs: Optional[List[str]] = None) -> List[Any]:
    """Execute Stage 2."""
    out_dir = Path(args.output_dir) / "stage2_validation"
    active_gates = [int(g.strip()) for g in args.gates.split(",") if g.strip()]
    early_stopping = not getattr(args, "no_early_stopping", False)

    if not candidate_pdbs:
        if not getattr(args, "candidates", None):
            logger.error("--candidates pattern or path is required for stage2!")
            sys.exit(1)
        raw_pattern = args.candidates
        if "*" in raw_pattern or "?" in raw_pattern:
            candidate_pdbs = sorted(glob.glob(raw_pattern))
        else:
            p = Path(raw_pattern)
            if p.is_dir():
                candidate_pdbs = sorted([str(f) for f in p.glob("*.pdb")])
            else:
                candidate_pdbs = [str(p)]

    if not candidate_pdbs:
        logger.error("No candidate PDB files found matching: %s", getattr(args, "candidates", ""))
        sys.exit(1)

    logger.info("Found %d candidate PDBs for Stage 2 validation.", len(candidate_pdbs))
    pipeline = Stage2ValidationPipeline(
        config=cfg,
        active_gates=active_gates,
        early_stopping=early_stopping,
        dry_run=args.dry_run,
    )

    ranked_cards = pipeline.validate_batch(
        candidate_pdbs=candidate_pdbs,
        output_dir=out_dir,
    )

    logger.info(
        "Stage 2 complete. Evaluated and ranked %d candidates. Summary: %s",
        len(ranked_cards),
        out_dir / "validation_summary.csv",
    )
    return ranked_cards


def run_stage3_cmd(
    args: argparse.Namespace,
    scaffold_pdb: str,
    ligand_name: str = "LIG",
    generations: int = 3,
    designs_per_gen: int = 15,
) -> Dict[str, Any]:
    """Execute Stage 3 AIzymes refinement."""
    out_dir = Path(args.output_dir) / "stage3_aizymes"
    bridge = AIzymesBridge()
    runner = Stage3AIzymesRunner(bridge=bridge, dry_run=args.dry_run)

    res = runner.run_refinement(
        scaffold_pdb=scaffold_pdb,
        output_dir=out_dir,
        generations=generations,
        designs_per_gen=designs_per_gen,
        dry_run=args.dry_run,
    )
    logger.info("Stage 3 complete. Refined variants persisted in %s", out_dir)
    return res


def main():
    print_banner()
    args = parse_args()

    # Load configuration and user credentials profile
    cfg = load_pipeline_config(
        config_path=args.config,
        user_profile=args.user,
        dry_run=args.dry_run,
    )

    logger.info("Operational Mode: %s", "PRE-FLIGHT DRY-RUN (ZERO COST)" if args.dry_run else "LIVE REMOTE COMPUTE")
    if args.user:
        logger.info("Authenticated User Profile: %s", args.user)

    cmd = args.command

    if cmd == "stage1":
        run_stage1_cmd(args, cfg)

    elif cmd == "stage2":
        run_stage2_cmd(args, cfg)

    elif cmd == "stage1-2":
        s1_res = run_stage1_cmd(args, cfg)
        with open(s1_res["manifest_path"], "r", encoding="utf-8") as f:
            manifest = json.load(f)
        pdbs = [c["pdb_path"] for c in manifest["candidates"]]
        run_stage2_cmd(args, cfg, candidate_pdbs=pdbs)

    elif cmd == "full":
        logger.info("=== STEP 1: De Novo Generative Design ===")
        s1_res = run_stage1_cmd(args, cfg)
        with open(s1_res["manifest_path"], "r", encoding="utf-8") as f:
            manifest = json.load(f)
        pdbs = [c["pdb_path"] for c in manifest["candidates"]]

        logger.info("=== STEP 2: 4-Gate Biophysical Validation & Triage ===")
        scorecards = run_stage2_cmd(args, cfg, candidate_pdbs=pdbs)

        pipeline = Stage2ValidationPipeline(config=cfg, dry_run=args.dry_run)
        top = pipeline.get_top_candidates(scorecards, n_top=1)
        best_pdb = top[0].pdb_path
        logger.info("Selected top validated candidate for AIzymes refinement: %s (rank #%d)", best_pdb, top[0].rank)

        logger.info("=== STEP 3: Stage 3 AIzymes Directed Evolution Refinement ===")
        run_stage3_cmd(
            args=args,
            scaffold_pdb=best_pdb,
            generations=args.generations,
            designs_per_gen=args.designs_per_gen,
        )
        logger.info("\n>>> Full Hybrid Pipeline Execution Successfully Completed! <<<")

    elif cmd == "aizymes-standalone":
        run_stage3_cmd(
            args=args,
            scaffold_pdb=args.structure,
            ligand_name=args.ligand_name,
            generations=args.generations,
            designs_per_gen=args.designs_per_gen,
        )

    sys.exit(0)


if __name__ == "__main__":
    main()
