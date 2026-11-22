"""AIzymes Bridge Module for Evolutionary Optimization.

Translates top-ranked de novo scaffolds from Stage 1/Stage 2 into the exact configuration,
directory layout, and constraint files expected by native AIzymes:
1. Formats the design scaffold as the starting WT PDB structure.
2. Identifies ligand residue code (e.g., LIG, 5TS, MUA).
3. Freezes catalytic triad residues to prevent destructive mutations.
4. Identifies second-shell residues (4.0 - 8.0 Angstroms from catalytic residues and ligand)
   as mutating positions for sequence optimization and pocket packing.
5. Generates Rosetta catalytic constraints (`enzdes.cst`) and FieldTools targets (`field_target.dat`).
6. Produces ready-to-run parameter dictionaries for `AIzymes_MAIN.setup()`.
"""

from __future__ import annotations

import logging
import math
import shutil
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

try:
    from hackathon_extensions.inputs.pipeline_input import LigandSpec, TheozymeResidueSpec, TheozymeSpec
except ImportError:
    LigandSpec = Any
    TheozymeResidueSpec = Any
    TheozymeSpec = Any

logger = logging.getLogger(__name__)


class AIzymesBridge:
    """Prepares and bridges design candidates for AIzymes directed evolution."""

    def __init__(
        self,
        second_shell_min_dist: float = 4.0,
        second_shell_max_dist: float = 8.0,
    ):
        self.second_shell_min_dist = second_shell_min_dist
        self.second_shell_max_dist = second_shell_max_dist

    def prepare_aizymes_workspace(
        self,
        scaffold_pdb: Union[str, Path],
        output_workspace: Union[str, Path],
        ligand_spec: Optional[LigandSpec] = None,
        theozyme_spec: Optional[TheozymeSpec] = None,
        catalytic_residues: Optional[List[int]] = None,
        design_methods: Optional[List[Any]] = None,
    ) -> Dict[str, Any]:
        """Build AIzymes project folder and return setup() parameters."""
        src_pdb = Path(scaffold_pdb)
        ws_path = Path(output_workspace)
        ws_path.mkdir(parents=True, exist_ok=True)

        # 1. Standardize WT PDB inside workspace
        wt_pdb_name = f"{src_pdb.stem}_WT.pdb"
        dest_wt_pdb = ws_path / wt_pdb_name
        shutil.copyfile(src_pdb, dest_wt_pdb)

        # 2. Extract ligand residue name
        ligand_name = ligand_spec.resname if ligand_spec and ligand_spec.resname else "LIG"

        # 3. Parse coordinates from PDB
        atoms_by_res, ligand_coords = self._parse_pdb_atoms(dest_wt_pdb, ligand_name)

        # 4. Resolve catalytic residue numbers
        cat_resnums: Set[int] = set()
        if catalytic_residues:
            cat_resnums.update(catalytic_residues)
        elif theozyme_spec and theozyme_spec.residues:
            for r in theozyme_spec.residues:
                cat_resnums.add(r.resnum)
        else:
            # Auto-detect catalytic candidates (CYS, SER, HIS, ASP) closest to ligand
            cat_resnums = self._detect_catalytic_residues(atoms_by_res, ligand_coords)

        # 5. Detect second-shell design residues (4 - 8 A from ligand & catalytic core)
        second_shell_res = self._find_second_shell_residues(
            atoms_by_res=atoms_by_res,
            ligand_coords=ligand_coords,
            cat_resnums=cat_resnums,
            min_dist=self.second_shell_min_dist,
            max_dist=self.second_shell_max_dist,
        )

        design_str = ",".join(str(r) for r in sorted(second_shell_res))
        logger.info(
            "Identified %d second-shell design residues for AIzymes: %s",
            len(second_shell_res),
            design_str,
        )

        # 6. Freeze catalytic residues in RESTRICT_RESIDUES
        restrict_residues: List[List[Any]] = []
        for r_num in sorted(cat_resnums):
            if r_num in atoms_by_res:
                res_3 = atoms_by_res[r_num]["resname"]
                res_1 = self._three_to_one(res_3)
                restrict_residues.append([r_num, res_1])

        # 7. Generate Rosetta catalytic constraint file enzdes.cst
        cst_path = ws_path / "enzdes.cst"
        self._write_enzdes_cst(
            output_path=cst_path,
            ligand_name=ligand_name,
            theozyme_spec=theozyme_spec,
            cat_resnums=cat_resnums,
            atoms_by_res=atoms_by_res,
        )

        # 8. Generate FieldTools target file field_target.dat
        field_target_path = ws_path / "field_target.dat"
        self._write_field_target_dat(
            output_path=field_target_path,
            ligand_name=ligand_name,
            cat_resnums=cat_resnums,
        )

        # Default design methods combining Rosetta and ProteinMPNN + FieldTools
        if design_methods is None:
            design_methods = [
                [0.7, "RosettaDesign", "ElectricFields"],
                [0.3, "ProteinMPNN", "ESMfold", "RosettaRelax", "ElectricFields"],
            ]

        # Consolidated setup parameters for AIzymes_MAIN.setup()
        setup_params = {
            "WT": str(dest_wt_pdb),
            "LIGAND": ligand_name,
            "DESIGN": design_str,
            "RESTRICT_RESIDUES": restrict_residues,
            "CST_NAME": str(cst_path),
            "FIELD_TARGET": str(field_target_path),
            "FOLDER_HOME": str(ws_path),
            "PARENT_DES_MED": ["RosettaDesign", "ElectricFields"],
            "DESIGN_METHODS": design_methods,
            "SELECTED_SCORES": ["total", "catalytic", "interface", "efield"],
            "WEIGHT_TOTAL": 1.0,
            "WEIGHT_CATALYTIC": 1.5,
            "WEIGHT_INTERFACE": 1.0,
            "WEIGHT_EFIELD": 1.5,
            "MAX_DESIGNS": 50,
            "KBT_BOLTZMANN": [1.0, 0.5, 0.02],
        }

        # Write manifest file for workspace tracking
        manifest_path = ws_path / "aizymes_bridge_manifest.json"
        import json
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump({
                "source_pdb": str(src_pdb),
                "wt_pdb": str(dest_wt_pdb),
                "ligand_name": ligand_name,
                "catalytic_residues_frozen": sorted(list(cat_resnums)),
                "second_shell_design_residues": sorted(list(second_shell_res)),
                "cst_file": str(cst_path),
                "field_target_file": str(field_target_path),
                "setup_params": setup_params,
            }, f, indent=2)

        return setup_params

    def _parse_pdb_atoms(
        self,
        pdb_path: Path,
        ligand_name: str,
    ) -> Tuple[Dict[int, Dict[str, Any]], List[List[float]]]:
        """Extract all atom coordinates grouped by residue number, and ligand atom coordinates."""
        atoms_by_res: Dict[int, Dict[str, Any]] = {}
        ligand_coords: List[List[float]] = []

        with open(pdb_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    res_num_str = line[22:26].strip()
                    try:
                        res_num = int(res_num_str)
                    except ValueError:
                        continue

                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                    except ValueError:
                        continue

                    coord = [x, y, z]

                    if res_name == ligand_name:
                        ligand_coords.append(coord)
                    else:
                        if res_num not in atoms_by_res:
                            atoms_by_res[res_num] = {
                                "resname": res_name,
                                "atoms": {},
                                "ca_coord": None,
                            }
                        atoms_by_res[res_num]["atoms"][atom_name] = coord
                        if atom_name == "CA":
                            atoms_by_res[res_num]["ca_coord"] = coord

        return atoms_by_res, ligand_coords

    def _detect_catalytic_residues(
        self,
        atoms_by_res: Dict[int, Dict[str, Any]],
        ligand_coords: List[List[float]],
    ) -> Set[int]:
        """Find catalytic candidate residues (CYS, SER, HIS, ASP) closest to ligand."""
        cat_candidates = {"CYS", "SER", "HIS", "ASP", "GLU", "LYS"}
        res_distances: List[Tuple[float, int]] = []

        for res_num, data in atoms_by_res.items():
            if data["resname"] in cat_candidates and data["atoms"]:
                # Compute min distance to any ligand atom
                min_d = 999.0
                if ligand_coords:
                    for atom_coord in data["atoms"].values():
                        for l_coord in ligand_coords:
                            d = self._dist(atom_coord, l_coord)
                            if d < min_d:
                                min_d = d
                else:
                    min_d = 3.0  # default proximity if no ligand atoms found
                res_distances.append((min_d, res_num))

        res_distances.sort(key=lambda x: x[0])
        # Return the top 3 closest catalytic residues (triad)
        return {r[1] for r in res_distances[:3]}

    def _find_second_shell_residues(
        self,
        atoms_by_res: Dict[int, Dict[str, Any]],
        ligand_coords: List[List[float]],
        cat_resnums: Set[int],
        min_dist: float,
        max_dist: float,
    ) -> Set[int]:
        """Find residues within min_dist - max_dist of catalytic residues or ligand."""
        second_shell: Set[int] = set()

        # Collect catalytic core reference coordinates
        core_coords: List[List[float]] = list(ligand_coords)
        for c_res in cat_resnums:
            if c_res in atoms_by_res:
                core_coords.extend(atoms_by_res[c_res]["atoms"].values())

        if not core_coords:
            # Fallback: take 8 residues around catalytic numbers if no coordinates
            for c_res in cat_resnums:
                for offset in [-3, -2, -1, 1, 2, 3]:
                    candidate = c_res + offset
                    if candidate in atoms_by_res and candidate not in cat_resnums:
                        second_shell.add(candidate)
            return second_shell

        for res_num, data in atoms_by_res.items():
            if res_num in cat_resnums:
                continue  # Never mutate frozen catalytic triad

            # Check min distance of any atom in this residue to any core coordinate
            min_r_dist = 999.0
            for atom_coord in data["atoms"].values():
                for c_coord in core_coords:
                    d = self._dist(atom_coord, c_coord)
                    if d < min_r_dist:
                        min_r_dist = d
                        if min_r_dist < min_dist:
                            break
                if min_r_dist < min_dist:
                    break

            # If within second shell zone
            if min_dist <= min_r_dist <= max_dist:
                second_shell.add(res_num)

        # Ensure at least 4 second-shell residues are selected
        if len(second_shell) < 4:
            for res_num in sorted(atoms_by_res.keys()):
                if res_num not in cat_resnums and res_num not in second_shell:
                    second_shell.add(res_num)
                    if len(second_shell) >= 6:
                        break

        return second_shell

    @staticmethod
    def _dist(p1: List[float], p2: List[float]) -> float:
        return math.sqrt(
            (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2
        )

    @staticmethod
    def _three_to_one(res3: str) -> str:
        mapping = {
            "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
            "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
            "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
            "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
        }
        return mapping.get(res3.upper(), "A")

    def _write_enzdes_cst(
        self,
        output_path: Path,
        ligand_name: str,
        theozyme_spec: Optional[TheozymeSpec],
        cat_resnums: Set[int],
        atoms_by_res: Dict[int, Dict[str, Any]],
    ) -> None:
        """Write Rosetta enzdes catalytic constraint file."""
        lines = [
            f"# Rosetta enzdes catalytic constraints for {ligand_name}",
            "# Generated by AIzymesBridge",
            "",
        ]

        if theozyme_spec and theozyme_spec.residues:
            for r in theozyme_spec.residues:
                cat_atom = r.target_atoms[0] if r.target_atoms else "SG"
                lig_atom = r.target_ligand_atom or "C1"
                min_d, max_d = r.target_distance_range
                mean_d = round((min_d + max_d) / 2.0, 2)
                lines.extend([
                    "CST::BEGIN",
                    f"  TEMPLATE::   ATOM_MAP: 1 atom_name: {lig_atom}",
                    f"  TEMPLATE::   ATOM_MAP: 1 residue3: {ligand_name}",
                    f"  TEMPLATE::   ATOM_MAP: 2 atom_name: {cat_atom}",
                    f"  TEMPLATE::   ATOM_MAP: 2 residue3: {r.resname}",
                    f"  CONSTRAINT:: distanceAB: {mean_d:.2f} 0.20 100.0 0 0",
                    "  CONSTRAINT:: angle_A: 109.5 10.0 50.0 360.0 0",
                    "  CONSTRAINT:: angle_B: 109.5 10.0 50.0 360.0 0",
                    "CST::END",
                    "",
                ])
        else:
            # Generic constraint block
            lines.extend([
                "CST::BEGIN",
                f"  TEMPLATE::   ATOM_MAP: 1 atom_name: C1",
                f"  TEMPLATE::   ATOM_MAP: 1 residue3: {ligand_name}",
                f"  TEMPLATE::   ATOM_MAP: 2 atom_name: SG",
                f"  TEMPLATE::   ATOM_MAP: 2 residue3: CYS",
                "  CONSTRAINT:: distanceAB: 3.10 0.20 100.0 0 0",
                "  CONSTRAINT:: angle_A: 109.5 10.0 50.0 360.0 0",
                "  CONSTRAINT:: angle_B: 109.5 10.0 50.0 360.0 0",
                "CST::END",
                "",
            ])

        with open(output_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines))

    def _write_field_target_dat(
        self,
        output_path: Path,
        ligand_name: str,
        cat_resnums: Set[int],
    ) -> None:
        """Write FieldTools target definition file defining reaction dipole axis."""
        lines = [
            f"# FieldTools Electric Field Target Axis for {ligand_name}",
            "# Format: atom1_res atom1_name atom2_res atom2_name target_field_MV_per_cm",
            f"1 C1 1 O1 20.0",
        ]
        with open(output_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
