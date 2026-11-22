"""Theozyme 3D Geometry Builder and Constraint Generator.

Constructs idealized catalytic amino acid side chains and active site geometries
around substrates or transition-state analogs. Positions:
- Nucleophiles (e.g. Cys Sγ or Ser Oγ) along the Bürgi–Dunitz trajectory (~107°).
- General bases (e.g. His Nε2/Nδ1, Asp/Glu carboxylates) positioned for proton transfer.
- Catalytic triads (Cys-His-Asp, Ser-His-Asp) and oxyanion hole donors.
- Exports coordinates in PDB and mmCIF formats.
- Generates Rosetta Match/enzdes (.cst) and RFD4 condition-spec geometric constraints.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

from hackathon_extensions.inputs.pipeline_input import (
    LigandSpec,
    TheozymeResidueSpec,
    TheozymeSpec,
)


class TheozymeBuilder:
    """Builds 3D catalytic theozyme geometries and restraint specifications."""

    # Canonical idealized residue coordinate templates centered around catalytic atoms
    _CANONICAL_RESIDUES: Dict[str, List[Tuple[str, str, Tuple[float, float, float]]]] = {
        "CYS": [
            ("N", "N", (-1.8, 1.2, 0.0)),
            ("CA", "C", (-1.2, 0.0, 0.0)),
            ("C", "C", (-1.9, -1.2, 0.2)),
            ("O", "O", (-2.6, -1.2, 1.2)),
            ("CB", "C", (-0.7, -0.1, -1.4)),
            ("SG", "S", (0.0, 0.0, 0.0)),  # Primary catalytic reference
        ],
        "SER": [
            ("N", "N", (-1.8, 1.2, 0.0)),
            ("CA", "C", (-1.2, 0.0, 0.0)),
            ("C", "C", (-1.9, -1.2, 0.2)),
            ("O", "O", (-2.6, -1.2, 1.2)),
            ("CB", "C", (-0.7, -0.1, -1.3)),
            ("OG", "O", (0.0, 0.0, 0.0)),  # Primary catalytic reference
        ],
        "HIS": [
            ("N", "N", (-2.5, 2.0, 0.5)),
            ("CA", "C", (-1.9, 0.7, 0.0)),
            ("C", "C", (-2.6, -0.4, 0.5)),
            ("O", "O", (-3.3, -0.3, 1.5)),
            ("CB", "C", (-1.5, 0.6, -1.5)),
            ("CG", "C", (-0.6, -0.5, -1.8)),
            ("ND1", "N", (-0.8, -1.8, -1.5)),
            ("CD2", "C", (0.6, -0.4, -2.4)),
            ("CE1", "C", (0.2, -2.4, -1.9)),
            ("NE2", "N", (0.0, 0.0, 0.0)),  # Primary catalytic reference
        ],
        "ASP": [
            ("N", "N", (-2.4, 1.8, 0.2)),
            ("CA", "C", (-1.8, 0.5, 0.0)),
            ("C", "C", (-2.4, -0.5, 0.8)),
            ("O", "O", (-3.0, -0.3, 1.8)),
            ("CB", "C", (-1.3, 0.3, -1.4)),
            ("CG", "C", (-0.5, -0.9, -1.7)),
            ("OD2", "O", (-0.6, -1.9, -1.0)),
            ("OD1", "O", (0.0, 0.0, 0.0)),  # Primary catalytic reference
        ],
        "GLU": [
            ("N", "N", (-3.0, 2.0, 0.5)),
            ("CA", "C", (-2.2, 0.8, 0.0)),
            ("C", "C", (-2.8, -0.3, 0.8)),
            ("O", "O", (-3.4, -0.2, 1.8)),
            ("CB", "C", (-1.8, 0.6, -1.4)),
            ("CG", "C", (-1.2, -0.7, -1.8)),
            ("CD", "C", (-0.6, -0.7, -3.2)),
            ("OE2", "O", (-0.8, -1.7, -3.9)),
            ("OE1", "O", (0.0, 0.0, 0.0)),  # Primary catalytic reference
        ],
    }

    def build_theozyme(
        self,
        ligand: LigandSpec,
        mechanism_type: str = "cysteine_hydrolase",
        custom_residues: Optional[List[TheozymeResidueSpec]] = None,
        output_pdb: Optional[Union[str, Path]] = None,
        output_cif: Optional[Union[str, Path]] = None,
        output_cst: Optional[Union[str, Path]] = None,
    ) -> TheozymeSpec:
        """Construct full 3D coordinates and constraints for the catalytic theozyme.

        Parameters
        ----------
        ligand : LigandSpec
            Ligand specification with 3D coordinates and reactive atom names.
        mechanism_type : str
            Catalytic mechanism: "cysteine_hydrolase", "serine_hydrolase", or "kemp_eliminase".
        custom_residues : Optional[List[TheozymeResidueSpec]]
            Explicit residue definitions if overriding standard triad templates.
        output_pdb : Optional[Union[str, Path]]
            Destination path for theozyme PDB.
        output_cif : Optional[Union[str, Path]]
            Destination path for theozyme mmCIF.
        output_cst : Optional[Union[str, Path]]
            Destination path for Rosetta .cst restraint file.

        Returns
        -------
        TheozymeSpec
            Complete theozyme specification with 3D coordinates and restraints.
        """
        if custom_residues:
            residue_specs = custom_residues
        else:
            residue_specs = self._get_default_mechanism_residues(mechanism_type, ligand)

        # Parse ligand coordinates to find reactive anchor
        lig_atoms = self._parse_pdb_atoms(ligand.conformer_pdb or "")
        reactive_atom_name = ligand.reactive_atoms[0] if ligand.reactive_atoms else "C1"
        anchor_coords = lig_atoms.get(reactive_atom_name, (0.0, 0.0, 0.0))

        # Position catalytic residues around anchor
        assembled_residue_blocks: List[str] = []
        cst_lines: List[str] = []
        atom_serial = 1

        # Track placed residues for geometric chaining (e.g. Cys -> His -> Asp)
        placed_positions: Dict[str, Tuple[float, float, float]] = {reactive_atom_name: anchor_coords}

        for r_spec in residue_specs:
            resname = r_spec.resname.upper()
            chain = r_spec.chain_id
            resnum = r_spec.resnum

            # Determine placement offset
            placement_origin = self._compute_residue_placement(
                r_spec=r_spec,
                placed_positions=placed_positions,
                mechanism_type=mechanism_type,
            )
            cat_atom = r_spec.target_atoms[0] if r_spec.target_atoms else "SG"
            placed_positions[f"{resname}:{cat_atom}"] = placement_origin

            # Translate canonical template to placement origin
            template = self._CANONICAL_RESIDUES.get(resname, self._CANONICAL_RESIDUES["CYS"])
            pdb_block_lines = []

            for atom_name, elem, (tx, ty, tz) in template:
                wx = placement_origin[0] + tx
                wy = placement_origin[1] + ty
                wz = placement_origin[2] + tz
                pdb_line = self._format_atom_line(
                    atom_serial=atom_serial,
                    atom_name=atom_name,
                    resname=resname,
                    chain=chain,
                    resnum=resnum,
                    x=wx, y=wy, z=wz,
                    element=elem,
                )
                pdb_block_lines.append(pdb_line)
                atom_serial += 1

            r_spec.pdb_coordinates = "\n".join(pdb_block_lines)
            assembled_residue_blocks.append(r_spec.pdb_coordinates)

        # Combine catalytic residues + ligand into complete theozyme PDB
        theozyme_pdb_lines: List[str] = [
            f"REMARK   Theozyme Model: {mechanism_type}",
            f"REMARK   Designed by Antigravity TheozymeBuilder",
        ]
        theozyme_pdb_lines.extend(assembled_residue_blocks)

        # Append ligand lines
        lig_lines = [
            line for line in (ligand.conformer_pdb or "").splitlines()
            if line.startswith(("HETATM", "ATOM  "))
        ]
        theozyme_pdb_lines.extend(lig_lines)
        theozyme_pdb_lines.append("END")
        full_theozyme_pdb = "\n".join(theozyme_pdb_lines) + "\n"

        # Generate geometric restraints (Rosetta .cst + RFD4 condition constraints)
        cst_blocks = self._generate_rosetta_cst(residue_specs, ligand, mechanism_type)
        rfd4_constraints = self._generate_rfd4_constraints(residue_specs, ligand)

        theozyme_spec = TheozymeSpec(
            residues=residue_specs,
            description=f"Catalytic arrangement for {mechanism_type}",
            mechanism_type=mechanism_type,
            geometric_restraints=cst_blocks + rfd4_constraints,
            theozyme_pdb=full_theozyme_pdb,
            theozyme_path=str(output_pdb) if output_pdb else None,
        )

        # Write files if paths requested
        if output_pdb:
            p = Path(output_pdb)
            p.parent.mkdir(parents=True, exist_ok=True)
            with open(p, "w", encoding="utf-8") as f:
                f.write(full_theozyme_pdb)

        if output_cif:
            cif_content = self.pdb_to_mmcif(full_theozyme_pdb)
            cp = Path(output_cif)
            cp.parent.mkdir(parents=True, exist_ok=True)
            with open(cp, "w", encoding="utf-8") as f:
                f.write(cif_content)

        if output_cst:
            cst_p = Path(output_cst)
            cst_p.parent.mkdir(parents=True, exist_ok=True)
            with open(cst_p, "w", encoding="utf-8") as f:
                f.write("\n".join(cst_blocks) + "\n")

        return theozyme_spec

    def _get_default_mechanism_residues(
        self,
        mechanism_type: str,
        ligand: LigandSpec,
    ) -> List[TheozymeResidueSpec]:
        """Provide standard catalytic residue specifications based on mechanism."""
        reactive_lig_atom = ligand.reactive_atoms[0] if ligand.reactive_atoms else "C1"

        if mechanism_type == "cysteine_hydrolase":
            return [
                TheozymeResidueSpec(
                    resname="CYS",
                    chain_id="A",
                    resnum=10,
                    target_atoms=["SG"],
                    target_ligand_atom=reactive_lig_atom,
                    target_distance_range=(2.7, 3.2),
                    target_angle_range=(105.0, 110.0),  # Bürgi-Dunitz angle
                    role="nucleophile",
                ),
                TheozymeResidueSpec(
                    resname="HIS",
                    chain_id="A",
                    resnum=40,
                    target_atoms=["NE2", "ND1"],
                    target_ligand_atom=None,
                    target_distance_range=(3.2, 3.8),
                    target_angle_range=(140.0, 180.0),
                    role="general_base",
                ),
                TheozymeResidueSpec(
                    resname="ASP",
                    chain_id="A",
                    resnum=85,
                    target_atoms=["OD1", "OD2"],
                    target_ligand_atom=None,
                    target_distance_range=(2.6, 3.0),
                    target_angle_range=(120.0, 160.0),
                    role="transition_state_polarizer",
                ),
            ]
        elif mechanism_type == "serine_hydrolase":
            return [
                TheozymeResidueSpec(
                    resname="SER",
                    chain_id="A",
                    resnum=15,
                    target_atoms=["OG"],
                    target_ligand_atom=reactive_lig_atom,
                    target_distance_range=(2.7, 3.1),
                    target_angle_range=(105.0, 110.0),
                    role="nucleophile",
                ),
                TheozymeResidueSpec(
                    resname="HIS",
                    chain_id="A",
                    resnum=45,
                    target_atoms=["NE2", "ND1"],
                    target_ligand_atom=None,
                    target_distance_range=(2.7, 3.1),
                    role="general_base",
                ),
                TheozymeResidueSpec(
                    resname="ASP",
                    chain_id="A",
                    resnum=90,
                    target_atoms=["OD1", "OD2"],
                    target_ligand_atom=None,
                    target_distance_range=(2.6, 2.9),
                    role="transition_state_polarizer",
                ),
            ]
        elif mechanism_type == "kemp_eliminase":
            return [
                TheozymeResidueSpec(
                    resname="ASP",
                    chain_id="A",
                    resnum=20,
                    target_atoms=["OD1"],
                    target_ligand_atom=reactive_lig_atom,
                    target_distance_range=(2.6, 3.0),
                    role="general_base",
                ),
            ]
        else:
            # Generic single nucleophile fallback
            return [
                TheozymeResidueSpec(
                    resname="CYS",
                    chain_id="A",
                    resnum=10,
                    target_atoms=["SG"],
                    target_ligand_atom=reactive_lig_atom,
                    target_distance_range=(2.8, 3.2),
                    role="catalytic_residue",
                )
            ]

    def _compute_residue_placement(
        self,
        r_spec: TheozymeResidueSpec,
        placed_positions: Dict[str, Tuple[float, float, float]],
        mechanism_type: str,
    ) -> Tuple[float, float, float]:
        """Compute the 3D position where a catalytic residue reference atom is located."""
        resname = r_spec.resname.upper()

        if r_spec.role == "nucleophile" or resname in ("CYS", "SER"):
            # Position along Bürgi-Dunitz trajectory (~107 deg) relative to ligand C1
            lig_key = r_spec.target_ligand_atom or "C1"
            lx, ly, lz = placed_positions.get(lig_key, (0.0, 0.0, 0.0))
            dist = (r_spec.target_distance_range[0] + r_spec.target_distance_range[1]) / 2.0
            # Vector attacking from top/front face
            bd_angle_rad = math.radians(107.0)
            dx = dist * math.sin(bd_angle_rad)
            dy = dist * math.cos(bd_angle_rad)
            dz = 1.0  # slight z-elevation for stereochemical clearance
            return (lx + dx, ly + dy, lz + dz)

        elif r_spec.role == "general_base" or resname == "HIS":
            # Position relative to nucleophile (CYS:SG or SER:OG)
            nuc_key = "CYS:SG" if "CYS:SG" in placed_positions else "SER:OG"
            nx, ny, nz = placed_positions.get(nuc_key, (2.8, 0.5, 1.0))
            # His NE2 placed ~3.4 A away from Cys SG to form catalytic dyad
            return (nx + 2.5, ny + 1.8, nz - 0.5)

        elif r_spec.role == "transition_state_polarizer" or resname in ("ASP", "GLU"):
            # Position relative to HIS:NE2
            hx, hy, hz = placed_positions.get("HIS:NE2", (5.3, 2.3, 0.5))
            # Asp OD1 placed ~2.8 A away from His ND1/CG backbone
            return (hx + 2.3, hy - 1.2, hz - 0.3)

        else:
            # Default offset along x-axis
            return (3.0, 0.0, 0.0)

    @staticmethod
    def _parse_pdb_atoms(pdb_text: str) -> Dict[str, Tuple[float, float, float]]:
        """Parse coordinates of atoms from PDB string."""
        atoms: Dict[str, Tuple[float, float, float]] = {}
        for line in pdb_text.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                name = line[12:16].strip()
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atoms[name] = (x, y, z)
                except ValueError:
                    continue
        return atoms

    @staticmethod
    def _format_atom_line(
        atom_serial: int,
        atom_name: str,
        resname: str,
        chain: str,
        resnum: int,
        x: float,
        y: float,
        z: float,
        element: str,
    ) -> str:
        """Format standard ATOM record for protein residues."""
        if len(atom_name) <= 3:
            formatted_name = f" {atom_name:<3s}"
        else:
            formatted_name = f"{atom_name:<4s}"

        return (
            f"ATOM  {atom_serial:5d} {formatted_name} {resname:3s} {chain:1s}{resnum:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {element:>2s}  "
        )

    def pdb_to_mmcif(self, pdb_text: str) -> str:
        """Convert PDB text into valid mmCIF string for rfd4-proteina-prep compatibility."""
        cif_lines: List[str] = [
            "data_theozyme",
            "#",
            "_entry.id theozyme",
            "#",
            "loop_",
            "_atom_site.group_PDB",
            "_atom_site.id",
            "_atom_site.type_symbol",
            "_atom_site.label_atom_id",
            "_atom_site.label_alt_id",
            "_atom_site.label_comp_id",
            "_atom_site.label_asym_id",
            "_atom_site.label_seq_id",
            "_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z",
            "_atom_site.occupancy",
            "_atom_site.B_iso_or_equiv",
            "_atom_site.auth_seq_id",
            "_atom_site.auth_comp_id",
            "_atom_site.auth_asym_id",
            "_atom_site.auth_atom_id",
        ]

        atom_id = 1
        for line in pdb_text.splitlines():
            if line.startswith(("ATOM  ", "HETATM")):
                group = line[0:6].strip()
                name = line[12:16].strip()
                comp = line[17:20].strip()
                chain = line[21:22].strip() or "A"
                seq_str = line[22:26].strip()
                seq = int(seq_str) if seq_str.isdigit() else 1
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    occ = float(line[54:60])
                    bfac = float(line[60:66])
                except (ValueError, IndexError):
                    x, y, z, occ, bfac = 0.0, 0.0, 0.0, 1.0, 20.0

                elem = line[76:78].strip()
                if not elem:
                    elem = name[0]

                cif_lines.append(
                    f"{group:<6s} {atom_id:5d} {elem:<2s} {name:<4s} . {comp:<3s} {chain:1s} {seq:4d} "
                    f"{x:8.3f} {y:8.3f} {z:8.3f} {occ:5.2f} {bfac:5.2f} {seq:4d} {comp:<3s} {chain:1s} {name:<4s}"
                )
                atom_id += 1

        cif_lines.append("#\n")
        return "\n".join(cif_lines)

    def _generate_rosetta_cst(
        self,
        residues: List[TheozymeResidueSpec],
        ligand: LigandSpec,
        mechanism_type: str,
    ) -> List[str]:
        """Generate Rosetta Match / Enzdes .cst constraint blocks."""
        blocks: List[str] = []
        reactive_lig_atom = ligand.reactive_atoms[0] if ligand.reactive_atoms else "C1"

        # Constraint 1: Ligand to Primary Catalytic Residue (e.g. Cys SG or Ser OG)
        if residues:
            nuc = residues[0]
            cat_atom = nuc.target_atoms[0] if nuc.target_atoms else "SG"
            min_d, max_d = nuc.target_distance_range
            mean_d = (min_d + max_d) / 2.0
            tol_d = (max_d - min_d) / 2.0

            block1 = (
                f"# Rosetta Enzdes Constraint: {ligand.resname} to {nuc.resname}\n"
                f"CST::BEGIN\n"
                f"  TEMPLATE::   ATOM_MAP: 1 atom_name: {reactive_lig_atom}\n"
                f"  TEMPLATE::   ATOM_MAP: 1 residue3: {ligand.resname}\n\n"
                f"  TEMPLATE::   ATOM_MAP: 2 atom_name: {cat_atom} CB CA\n"
                f"  TEMPLATE::   ATOM_MAP: 2 residue3: {nuc.resname}\n\n"
                f"  CONSTRAINT:: distanceAB:    {mean_d:6.2f}   {tol_d:5.2f}  100.00  1\n"
                f"  CONSTRAINT:: angle_A:     107.00  10.00   50.00  360.00  1\n"
                f"  CONSTRAINT:: angle_B:     109.50  10.00   50.00  360.00  1\n"
                f"  CONSTRAINT:: torsion_A:   180.00  30.00   20.00  360.00  1\n"
                f"  CONSTRAINT:: torsion_AB:  180.00  30.00   20.00  360.00  1\n"
                f"  CONSTRAINT:: torsion_B:   180.00  30.00   20.00  360.00  1\n"
                f"CST::END"
            )
            blocks.append(block1)

        # Constraint 2: Catalytic Dyad (e.g. Cys SG to His NE2)
        if len(residues) >= 2:
            r1, r2 = residues[0], residues[1]
            a1 = r1.target_atoms[0] if r1.target_atoms else "SG"
            a2 = r2.target_atoms[0] if r2.target_atoms else "NE2"
            min_d, max_d = r2.target_distance_range
            mean_d = (min_d + max_d) / 2.0
            tol_d = (max_d - min_d) / 2.0

            block2 = (
                f"# Rosetta Enzdes Constraint: {r1.resname} to {r2.resname} (Catalytic Dyad)\n"
                f"CST::BEGIN\n"
                f"  TEMPLATE::   ATOM_MAP: 1 atom_name: {a1} CB CA\n"
                f"  TEMPLATE::   ATOM_MAP: 1 residue3: {r1.resname}\n\n"
                f"  TEMPLATE::   ATOM_MAP: 2 atom_name: {a2} CD2 CG\n"
                f"  TEMPLATE::   ATOM_MAP: 2 residue3: {r2.resname}\n\n"
                f"  CONSTRAINT:: distanceAB:    {mean_d:6.2f}   {tol_d:5.2f}  100.00  1\n"
                f"  CONSTRAINT:: angle_A:     105.00  15.00   50.00  360.00  1\n"
                f"  CONSTRAINT:: angle_B:     125.00  15.00   50.00  360.00  1\n"
                f"CST::END"
            )
            blocks.append(block2)

        return blocks

    def _generate_rfd4_constraints(
        self,
        residues: List[TheozymeResidueSpec],
        ligand: LigandSpec,
    ) -> List[str]:
        """Generate human-readable RFD4 condition-spec restraint descriptions."""
        lines: List[str] = ["# RFD4 Condition-Spec Constraints:"]
        reactive_lig_atom = ligand.reactive_atoms[0] if ligand.reactive_atoms else "C1"

        for r in residues:
            target_atom = r.target_atoms[0] if r.target_atoms else "SG"
            min_d, max_d = r.target_distance_range
            lines.append(
                f"- restraint: {r.chain_id}{r.resnum}:{target_atom} <-> {ligand.chain_id}{ligand.resnum}:{reactive_lig_atom} "
                f"dist=[{min_d:.2f}, {max_d:.2f}] A"
            )
        return lines
