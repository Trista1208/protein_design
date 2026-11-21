"""Case 1 Parser: Annotated mmCIF / PDB with Transition State & Catalytic Motif.

Handles inputs where a solved or modeled transition-state analog structure
contains pre-annotated catalytic residues (e.g., Cys-His-Asp or Ser-His-Asp).
- Parses mmCIF and PDB structures.
- Validates atom naming, chain IDs, and steric clashes.
- Generates:
  1. RFD4-Proteina condition-spec YAML (fixed motif atoms, buried flags, target scaffold length).
  2. Clean annotated mmCIF (.cif) ready for rfd4-proteina-prep.
  3. RFD3 JSON configuration with segment contigs (e.g. [15-35/A12-14/25-45/A50-52/15-35]).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from hackathon_extensions.inputs.pipeline_input import (
    CaseType,
    LigandSpec,
    PipelineInput,
    TheozymeResidueSpec,
    TheozymeSpec,
)


class Case1AnnotatedParser:
    """Parser and configuration generator for pre-annotated holo-complexes."""

    def __init__(
        self,
        default_scaffold_length_range: Tuple[int, int] = (150, 300),
    ):
        self.default_scaffold_length_range = default_scaffold_length_range

    def parse(
        self,
        structure_input: Union[str, Path],
        catalytic_residues: Optional[List[Dict[str, Any]]] = None,
        ligand_chain: str = "Z",
        ligand_resnum: int = 1,
        ligand_resname: str = "LIG",
        scaffold_length_range: Optional[Tuple[int, int]] = None,
        output_dir: Optional[Union[str, Path]] = None,
    ) -> PipelineInput:
        """Parse annotated PDB or mmCIF structure into unified PipelineInput.

        Parameters
        ----------
        structure_input : Union[str, Path]
            PDB / mmCIF file path or raw coordinate string.
        catalytic_residues : Optional[List[Dict[str, Any]]]
            List of dicts specifying catalytic residues, e.g.:
            [{"chain": "A", "resnum": 12, "resname": "CYS", "fixed_atoms": ["SG", "CB", "CA"]}, ...]
            If None, attempts auto-detection of close protein contacts to ligand.
        ligand_chain : str
            Chain ID of the ligand (default "Z").
        ligand_resnum : int
            Residue number of the ligand (default 1).
        ligand_resname : str
            Residue name of the ligand (default "LIG").
        scaffold_length_range : Optional[Tuple[int, int]]
            Min and max amino acid length for de novo scaffold generation.
        output_dir : Optional[Union[str, Path]]
            If provided, saves generated CIF, RFD4 YAML, and RFD3 JSON to disk.

        Returns
        -------
        PipelineInput
            Complete input configuration object.
        """
        raw_text = self._load_structure_text(structure_input)
        is_cif = "_atom_site" in raw_text or str(structure_input).endswith((".cif", ".mmcif"))

        # Parse atoms from structure
        protein_atoms, ligand_atoms = self._parse_atoms(
            raw_text=raw_text,
            is_cif=is_cif,
            ligand_chain=ligand_chain,
            ligand_resnum=ligand_resnum,
            ligand_resname=ligand_resname,
        )

        if not ligand_atoms:
            # Fallback search for any non-standard residue if exact match missing
            protein_atoms, ligand_atoms = self._auto_locate_ligand(raw_text, is_cif)

        # Validate steric clashing
        clashes = self._detect_clashes(protein_atoms, ligand_atoms)

        # Resolve catalytic residues
        cat_residue_specs = self._resolve_catalytic_residues(
            catalytic_residues=catalytic_residues,
            protein_atoms=protein_atoms,
            ligand_atoms=ligand_atoms,
        )

        # Build LigandSpec
        lig_pdb = self._atoms_to_pdb(ligand_atoms)
        reactive_atoms = self._find_reactive_ligand_atoms(ligand_atoms)
        ligand_spec = LigandSpec(
            conformer_pdb=lig_pdb,
            resname=ligand_resname,
            chain_id=ligand_chain,
            resnum=ligand_resnum,
            reactive_atoms=reactive_atoms,
            atom_metadata={"atom_count": len(ligand_atoms), "clashes": len(clashes)},
        )

        # Build TheozymeSpec
        theozyme_pdb = self._atoms_to_pdb(
            [a for a in protein_atoms if any(
                a["chain"] == r.chain_id and a["resnum"] == r.resnum for r in cat_residue_specs
            )] + ligand_atoms
        )
        theozyme_spec = TheozymeSpec(
            residues=cat_residue_specs,
            description="Annotated catalytic theozyme extracted from structure",
            mechanism_type="annotated_case1",
            theozyme_pdb=theozyme_pdb,
        )

        # Generate RFD4 condition-spec
        scaffold_range = scaffold_length_range or self.default_scaffold_length_range
        rfd4_yaml_dict = self._build_rfd4_condition_spec(
            cat_residue_specs=cat_residue_specs,
            ligand_spec=ligand_spec,
            scaffold_range=scaffold_range,
        )

        # Generate RFD3 JSON configuration
        rfd3_json_dict = self._build_rfd3_config(
            cat_residue_specs=cat_residue_specs,
            ligand_spec=ligand_spec,
            scaffold_range=scaffold_range,
        )

        # Format clean mmCIF
        clean_cif = self._build_clean_mmcif(cat_residue_specs, protein_atoms, ligand_atoms)

        # Save to output directory if specified
        if output_dir:
            out_p = Path(output_dir)
            out_p.mkdir(parents=True, exist_ok=True)

            with open(out_p / "annotated_input.cif", "w", encoding="utf-8") as f:
                f.write(clean_cif)

            with open(out_p / "rfd4_condition_spec.yaml", "w", encoding="utf-8") as f:
                f.write(self._dump_yaml(rfd4_yaml_dict))

            with open(out_p / "rfd3_config.json", "w", encoding="utf-8") as f:
                json.dump(rfd3_json_dict, f, indent=2)

        return PipelineInput(
            case_type=CaseType.CASE1_ANNOTATED,
            ligand=ligand_spec,
            theozyme=theozyme_spec,
            scaffold_source_pdb=raw_text if not is_cif else None,
            rfd4_condition_spec=rfd4_yaml_dict,
            rfd3_config=rfd3_json_dict,
            metadata={
                "parsed_protein_atoms": len(protein_atoms),
                "parsed_ligand_atoms": len(ligand_atoms),
                "clash_warnings": clashes,
            },
        )

    def parse_rfd3_json(
        self,
        json_input: Union[str, Path],
        output_dir: Optional[Union[str, Path]] = None,
    ) -> PipelineInput:
        """Parse an RFdiffusion3 task JSON specification (e.g. 1qji_json.json.txt).

        Reads the referenced PDB/CIF, resolves WSL vs Windows paths, extracts
        catalytic residues from 'unindex' and 'select_fixed_atoms', parses ligands,
        and constructs both the RFD3 task config and the RFD4-Proteina condition-spec.
        """
        raw_text = self._load_structure_text(json_input)

        import re
        m = re.search(r"\{.*\}", raw_text, re.DOTALL)
        if not m:
            raise ValueError(f"No JSON object found in {json_input}")
        data = json.loads(m.group(0))

        task_key = list(data.keys())[0]
        task = data[task_key]

        raw_input_path = task.get("input", "")
        clean_path = raw_input_path.replace("/mnt/c/", "C:/").replace("/mnt/d/", "D:/")
        struct_path = Path(clean_path)
        if not struct_path.is_file():
            candidates = [
                Path(json_input).parent / Path(raw_input_path).name,
                Path("test_pipeline") / Path(raw_input_path).name,
                Path(raw_input_path).name,
            ]
            for cand in candidates:
                if cand.is_file():
                    struct_path = cand
                    break

        if not struct_path.is_file():
            raise FileNotFoundError(f"Could not locate structure file referenced in JSON: {raw_input_path}")

        cif_text = self._load_structure_text(struct_path)
        is_cif = "_atom_site" in cif_text or str(struct_path).endswith((".cif", ".mmcif"))

        json_ligands = [s for s in str(task.get("ligand", "")).split(",") if s.strip()]
        all_protein_atoms, all_ligand_atoms = self._auto_locate_ligand(
            cif_text, is_cif, ligand_resnames=json_ligands
        )

        length_str = str(task.get("length", "120-150"))
        if "-" in length_str:
            parts = length_str.split("-")
            scaffold_range = (int(parts[0]), int(parts[1]))
        else:
            scaffold_range = (int(length_str), int(length_str))

        unindex_strs = [s.strip() for s in task.get("unindex", "").split(",") if s.strip()]
        fixed_atoms_map = task.get("select_fixed_atoms", {})

        cat_residue_specs: List[TheozymeResidueSpec] = []
        for u in unindex_strs:
            ch = u[0]
            num = int(u[1:])
            match_atoms = [a for a in all_protein_atoms if a["chain"] == ch and a["resnum"] == num]
            resname = match_atoms[0]["resname"] if match_atoms else "HIS"
            fixed_str = fixed_atoms_map.get(u, "")
            fixed_list = [a.strip() for a in fixed_str.split(",") if a.strip()]
            cat_residue_specs.append(
                TheozymeResidueSpec(
                    resname=resname,
                    chain_id=ch,
                    resnum=num,
                    target_atoms=fixed_list,
                )
            )

        lig_pdb = self._atoms_to_pdb(all_ligand_atoms)
        lig_resname = all_ligand_atoms[0]["resname"] if all_ligand_atoms else "LIG"
        lig_chain = all_ligand_atoms[0]["chain"] if all_ligand_atoms else "Z"
        lig_resnum = all_ligand_atoms[0]["resnum"] if all_ligand_atoms else 1

        ligand_spec = LigandSpec(
            conformer_pdb=lig_pdb,
            resname=lig_resname,
            chain_id=lig_chain,
            resnum=lig_resnum,
            reactive_atoms=self._find_reactive_ligand_atoms(all_ligand_atoms),
            atom_metadata={"atom_count": len(all_ligand_atoms), "source_task": task_key},
        )

        theozyme_spec = TheozymeSpec(
            residues=cat_residue_specs,
            description=f"Catalytic theozyme from RFD3 config ({task_key})",
            mechanism_type="rfd3_json_annotated",
            theozyme_pdb=self._atoms_to_pdb(all_protein_atoms + all_ligand_atoms),
        )

        rfd4_yaml_dict = self._build_rfd4_condition_spec(
            cat_residue_specs=cat_residue_specs,
            ligand_spec=ligand_spec,
            scaffold_range=scaffold_range,
        )
        if "diffusion_params" not in rfd4_yaml_dict:
            rfd4_yaml_dict["diffusion_params"] = {}
        rfd4_yaml_dict["diffusion_params"]["ori_token"] = task.get("ori_token", [])
        rfd4_yaml_dict["diffusion_params"]["allow_ligand_on_existing_chain"] = task.get(
            "allow_ligand_on_existing_chain", True
        )

        if output_dir:
            out_p = Path(output_dir)
            out_p.mkdir(parents=True, exist_ok=True)
            with open(out_p / "rfd4_condition_spec.yaml", "w", encoding="utf-8") as f:
                f.write(self._dump_yaml(rfd4_yaml_dict))
            with open(out_p / "rfd3_config.json", "w", encoding="utf-8") as f:
                json.dump(task, f, indent=2)

        return PipelineInput(
            case_type=CaseType.CASE1_ANNOTATED,
            ligand=ligand_spec,
            theozyme=theozyme_spec,
            scaffold_source_pdb=None,
            rfd4_condition_spec=rfd4_yaml_dict,
            rfd3_config=task,
            metadata={
                "source_json_task": task_key,
                "parsed_protein_atoms": len(all_protein_atoms),
                "parsed_ligand_atoms": len(all_ligand_atoms),
                "unindex": unindex_strs,
                "length_bounds": scaffold_range,
            },
        )

    @staticmethod
    def _load_structure_text(structure_input: Union[str, Path]) -> str:
        """Load structure content from file path or raw string."""
        if isinstance(structure_input, Path):
            with open(structure_input, "r", encoding="utf-8", errors="ignore") as f:
                return f.read()
        if isinstance(structure_input, str):
            try:
                p = Path(structure_input)
                if p.is_file():
                    with open(p, "r", encoding="utf-8", errors="ignore") as f:
                        return f.read()
            except Exception:
                pass
        return str(structure_input)

    def _parse_atoms(
        self,
        raw_text: str,
        is_cif: bool,
        ligand_chain: str,
        ligand_resnum: int,
        ligand_resname: str,
    ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Extract protein and ligand atoms from PDB or mmCIF."""
        protein_atoms: List[Dict[str, Any]] = []
        ligand_atoms: List[Dict[str, Any]] = []

        if is_cif:
            # Header-aware mmCIF parsing: column order differs between writers
            # (e.g. PyMOL emits label_entity_id), so never rely on fixed indices.
            cols: List[str] = []
            for line in raw_text.splitlines():
                s = line.strip()
                if s.startswith("_atom_site."):
                    cols.append(s.split()[0][len("_atom_site."):])
                    continue
                if not (s.startswith("ATOM") or s.startswith("HETATM")):
                    continue
                parts = s.split()
                if cols:
                    if len(parts) < len(cols):
                        continue
                    row = dict(zip(cols, parts))
                    atom_name = row.get("label_atom_id", "").strip('"').strip("'")
                    comp = row.get("label_comp_id", "")
                    chain = row.get("auth_asym_id") or row.get("label_asym_id", "A")
                    seq = row.get("auth_seq_id")
                    if seq in (None, "?", "."):
                        seq = row.get("label_seq_id")
                    elem = row.get("type_symbol", atom_name[:1])
                    try:
                        resnum = int(seq)
                        x = float(row["Cartn_x"])
                        y = float(row["Cartn_y"])
                        z = float(row["Cartn_z"])
                    except (ValueError, KeyError, TypeError):
                        continue
                    group = row.get("group_PDB", parts[0])
                else:
                    # Legacy fixed-column fallback (no _atom_site header present)
                    if len(parts) < 11:
                        continue
                    atom_name = parts[3].strip('"').strip("'")
                    comp = parts[5].strip('"').strip("'")
                    chain = parts[6].strip('"').strip("'")
                    try:
                        resnum = int(parts[7])
                        x, y, z = float(parts[8]), float(parts[9]), float(parts[10])
                    except (ValueError, IndexError):
                        continue
                    elem = parts[2]
                    group = parts[0]

                entry = {
                    "name": atom_name,
                    "resname": comp,
                    "chain": chain,
                    "resnum": resnum,
                    "x": x, "y": y, "z": z,
                    "element": elem.capitalize(),
                    "is_hetatm": group == "HETATM",
                }
                if chain == ligand_chain or comp == ligand_resname:
                    ligand_atoms.append(entry)
                else:
                    protein_atoms.append(entry)
        else:
            for line in raw_text.splitlines():
                if line.startswith(("ATOM  ", "HETATM")):
                    record = line[0:6].strip()
                    atom_name = line[12:16].strip()
                    comp = line[17:20].strip()
                    chain = line[21:22].strip() or "A"
                    seq_str = line[22:26].strip()
                    resnum = int(seq_str) if seq_str.isdigit() else 1
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except ValueError:
                        continue
                    elem = line[76:78].strip() or atom_name[0]

                    entry = {
                        "name": atom_name,
                        "resname": comp,
                        "chain": chain,
                        "resnum": resnum,
                        "x": x, "y": y, "z": z,
                        "element": elem.capitalize(),
                        "is_hetatm": record == "HETATM",
                    }
                    if chain == ligand_chain or comp == ligand_resname or (record == "HETATM" and comp not in ("HOH", "WAT")):
                        ligand_atoms.append(entry)
                    else:
                        protein_atoms.append(entry)

        return protein_atoms, ligand_atoms

    def _auto_locate_ligand(
        self,
        raw_text: str,
        is_cif: bool,
        ligand_resnames: Optional[List[str]] = None,
    ) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """Identify ligand atoms. Explicit residue names (e.g. 'ZN,PKF' from the task JSON)
        take priority; otherwise fall back to the first non-standard residue."""
        p_atoms, l_atoms = self._parse_atoms(raw_text, is_cif, ligand_chain="Z", ligand_resnum=1, ligand_resname="LIG")
        if ligand_resnames:
            wanted = {n.strip().upper() for n in ligand_resnames if n.strip()}
            moved = [a for a in p_atoms if a["resname"].upper() in wanted]
            if moved:
                l_atoms = l_atoms + moved
                p_atoms = [a for a in p_atoms if a["resname"].upper() not in wanted]
                return p_atoms, l_atoms
        if not l_atoms and p_atoms:
            # Check for non-standard residue
            standard_aa = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
            non_standard = [a for a in p_atoms if a["resname"] not in standard_aa and a["resname"] not in ("HOH", "WAT")]
            if non_standard:
                chosen_res = non_standard[0]["resname"]
                l_atoms = [a for a in non_standard if a["resname"] == chosen_res]
                p_atoms = [a for a in p_atoms if a["resname"] != chosen_res]
        return p_atoms, l_atoms

    @staticmethod
    def _detect_clashes(
        protein_atoms: List[Dict[str, Any]],
        ligand_atoms: List[Dict[str, Any]],
        threshold: float = 0.85,
    ) -> List[str]:
        """Detect unphysical steric clashes between atoms."""
        clashes: List[str] = []
        for p in protein_atoms:
            px, py, pz = p["x"], p["y"], p["z"]
            for l in ligand_atoms:
                dist = math.sqrt((px - l["x"])**2 + (py - l["y"])**2 + (pz - l["z"])**2)
                if dist < threshold:
                    clashes.append(
                        f"Clash: {p['chain']}{p['resnum']}:{p['name']} and {l['chain']}{l['resnum']}:{l['name']} (dist={dist:.2f} A)"
                    )
        return clashes

    def _resolve_catalytic_residues(
        self,
        catalytic_residues: Optional[List[Dict[str, Any]]],
        protein_atoms: List[Dict[str, Any]],
        ligand_atoms: List[Dict[str, Any]],
    ) -> List[TheozymeResidueSpec]:
        """Map user-specified or contact-based residues into TheozymeResidueSpec."""
        specs: List[TheozymeResidueSpec] = []

        if catalytic_residues:
            for cat in catalytic_residues:
                specs.append(
                    TheozymeResidueSpec(
                        resname=cat.get("resname", "CYS"),
                        chain_id=cat.get("chain", "A"),
                        resnum=int(cat.get("resnum", 1)),
                        target_atoms=cat.get("fixed_atoms", ["SG", "CB", "CA"]),
                        role=cat.get("role", "catalytic_motif"),
                    )
                )
        else:
            # Auto-detect protein residues closest to ligand (within 4.0 A)
            contact_res_keys: Set[Tuple[str, int, str]] = set()
            for l in ligand_atoms:
                lx, ly, lz = l["x"], l["y"], l["z"]
                for p in protein_atoms:
                    dist = math.sqrt((p["x"] - lx)**2 + (p["y"] - ly)**2 + (p["z"] - lz)**2)
                    if dist <= 3.8:
                        contact_res_keys.add((p["chain"], p["resnum"], p["resname"]))

            for chain, resnum, comp in sorted(contact_res_keys):
                specs.append(
                    TheozymeResidueSpec(
                        resname=comp,
                        chain_id=chain,
                        resnum=resnum,
                        target_atoms=["CA", "CB"],
                        role="contact_motif",
                    )
                )

        return specs

    @staticmethod
    def _find_reactive_ligand_atoms(ligand_atoms: List[Dict[str, Any]]) -> List[str]:
        """Identify candidate reactive atoms (e.g. carbonyl carbons C1, heteroatoms)."""
        reactive = []
        for a in ligand_atoms:
            if a["element"] == "C" and ("1" in a["name"] or a["name"] in ("C1", "C", "CA")):
                reactive.append(a["name"])
        if not reactive and ligand_atoms:
            reactive.append(ligand_atoms[0]["name"])
        return reactive

    @staticmethod
    def _atoms_to_pdb(atoms: List[Dict[str, Any]]) -> str:
        """Convert list of atom dicts to PDB string."""
        lines = []
        for idx, a in enumerate(atoms, start=1):
            record = "HETATM" if a.get("is_hetatm") else "ATOM  "
            name = a["name"]
            fmt_name = f" {name:<3s}" if len(name) <= 3 else f"{name:<4s}"
            lines.append(
                f"{record}{idx:5d} {fmt_name} {a['resname']:3s} {a['chain']:1s}{a['resnum']:4d}    "
                f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}  1.00 20.00          {a['element']:>2s}  "
            )
        lines.append("END")
        return "\n".join(lines) + "\n"

    def _build_rfd4_condition_spec(
        self,
        cat_residue_specs: List[TheozymeResidueSpec],
        ligand_spec: LigandSpec,
        scaffold_range: Tuple[int, int],
    ) -> Dict[str, Any]:
        """Construct RFD4 condition-spec dictionary."""
        fixed_motifs = []
        for r in cat_residue_specs:
            fixed_motifs.append({
                "chain": r.chain_id,
                "resnum": r.resnum,
                "resname": r.resname,
                "fixed_atoms": r.target_atoms or ["CA", "CB"],
                "role": r.role,
            })

        return {
            "version": "rfd4-proteina-v1",
            "target": {
                "scaffold_length_range": list(scaffold_range),
                "default_length": int((scaffold_range[0] + scaffold_range[1]) / 2),
            },
            "ligand": {
                "chain": ligand_spec.chain_id,
                "resname": ligand_spec.resname,
                "resnum": ligand_spec.resnum,
                "fixed_atoms": ["ALL"],
                "select_buried": ["ALL"],
                "reactive_atoms": ligand_spec.reactive_atoms,
            },
            "fixed_motifs": fixed_motifs,
            "diffusion_params": {
                "num_timesteps": 50,
                "noise_scale": 1.0,
            },
        }

    def _build_rfd3_config(
        self,
        cat_residue_specs: List[TheozymeResidueSpec],
        ligand_spec: LigandSpec,
        scaffold_range: Tuple[int, int],
    ) -> Dict[str, Any]:
        """Construct RFD3 / RFdiffusion configuration with segment contigs."""
        # Build contig segments e.g. [15-35/A12-14/25-45/A50-52/15-35]
        contig_parts: List[str] = []
        fixed_atom_map: Dict[str, List[str]] = {}

        contig_parts.append("15-35")
        for r in cat_residue_specs:
            motif_range = f"{r.chain_id}{r.resnum}-{r.resnum}"
            contig_parts.append(motif_range)
            contig_parts.append("25-45")
            fixed_atom_map[f"{r.chain_id}{r.resnum}"] = r.target_atoms or ["CA", "CB"]

        contigs_expr = f"[{'/'.join(contig_parts)}]"

        return {
            "contigmap": {
                "contigs": [contigs_expr],
                "inpaint_seq": None,
                "provide_seq": None,
            },
            "inference": {
                "ligand": f"{ligand_spec.chain_id}{ligand_spec.resnum}",
                "select_fixed_atoms": fixed_atom_map,
                "select_buried": {f"{ligand_spec.chain_id}{ligand_spec.resnum}": "ALL"},
            },
            "scaffold_length_bounds": list(scaffold_range),
        }

    def _build_clean_mmcif(
        self,
        cat_residue_specs: List[TheozymeResidueSpec],
        protein_atoms: List[Dict[str, Any]],
        ligand_atoms: List[Dict[str, Any]],
    ) -> str:
        """Generate annotated mmCIF format compliant with rfd4-proteina-prep."""
        lines = [
            "data_case1_annotated",
            "#",
            "_entry.id case1_annotated",
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
        all_atoms = protein_atoms + ligand_atoms
        for a in all_atoms:
            group = "HETATM" if a.get("is_hetatm") else "ATOM"
            name = a["name"]
            comp = a["resname"]
            chain = a["chain"]
            seq = a["resnum"]
            elem = a["element"]
            lines.append(
                f"{group:<6s} {atom_id:5d} {elem:<2s} {name:<4s} . {comp:<3s} {chain:1s} {seq:4d} "
                f"{a['x']:8.3f} {a['y']:8.3f} {a['z']:8.3f}  1.00 20.00 {seq:4d} {comp:<3s} {chain:1s} {name:<4s}"
            )
            atom_id += 1

        lines.append("#\n")
        return "\n".join(lines)

    @staticmethod
    def _dump_yaml(data: Dict[str, Any], indent_level: int = 0) -> str:
        """Lightweight YAML dumper avoiding external pyyaml dependency."""
        indent = "  " * indent_level
        lines = []
        for k, v in data.items():
            if isinstance(v, dict):
                lines.append(f"{indent}{k}:")
                lines.append(Case1AnnotatedParser._dump_yaml(v, indent_level + 1))
            elif isinstance(v, list):
                if all(isinstance(x, (int, float, str)) for x in v) and len(v) <= 4:
                    # Inline format e.g. [150, 300]
                    lines.append(f"{indent}{k}: {json.dumps(v)}")
                else:
                    lines.append(f"{indent}{k}:")
                    for item in v:
                        if isinstance(item, dict):
                            sub_str = Case1AnnotatedParser._dump_yaml(item, indent_level + 2).lstrip()
                            lines.append(f"{indent}  - {sub_str}")
                        else:
                            lines.append(f"{indent}  - {item}")
            elif isinstance(v, bool):
                lines.append(f"{indent}{k}: {'true' if v else 'false'}")
            elif v is None:
                lines.append(f"{indent}{k}: null")
            elif isinstance(v, (int, float)):
                lines.append(f"{indent}{k}: {v}")
            else:
                lines.append(f"{indent}{k}: \"{v}\"")
        return "\n".join(lines)
