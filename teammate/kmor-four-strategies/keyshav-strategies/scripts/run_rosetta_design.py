"""Enzdes-constrained FastDesign on a folded Kemp design.

    python scripts/run_rosetta_design.py <input.pdb> [nstruct] [label]

This is the tool for the gap the screen exposed: single mutants fix the
oxyanion DISTANCE or the TORSION, never both (LOOP_DESIGN 6.2 and the valid-set
analysis), because no generative tool we have optimises our objective. A
Rosetta enzdes block scores distance + angle + torsion in one energy.

WHAT THE REMARK 666 LINES ARE FOR
---------------------------------
Enzdes does not guess which residues a constraint block refers to. It reads
`REMARK 666 MATCH TEMPLATE <chain> <lig> <num> MATCH MOTIF <chain> <res> <num>
<block> <exgeom>` from the input PDB. Block 1 of our .cst is the Asp base,
block 2 the Gln donor, so those two lines are written here from the known
catalytic indices rather than left to automatic matching - which would happily
attach block 2 to a different glutamine and report a fine energy for the wrong
site. That is the same failure class as the ligand atom-name collision caught
in make_6nt_params.py.

CAVEAT THAT DOES NOT GO AWAY
----------------------------
Rosetta optimises ROSETTA energy. The competition metric is computed on a
BOLTZ-2 prediction. So every design out of here is a PROPOSAL and has to
survive a refold before it means anything (LESSONS 6: a coordinate-space
rotamer scan found 360 clash-free in-band rotamers and 0/5 Boltz seeds ever
sampled one). Use scripts/score_finalist.py on the outputs.
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ROSETTA = Path("/home/mandalorian/opt/rosetta")
BIN = ROSETTA / "source/bin/rosetta_scripts.linuxgccrelease"
DB = ROSETTA / "database"
RO = ROOT / "outputs/rosetta"
import os as _os
#: which protocol to run. enzdes_loop.py sets ENZDES_XML to the targeted
#: variant; without this the targeted XML would be silently ignored and every
#: run would redesign the whole 8 A shell instead of the named positions.
XML = RO / _os.environ.get("ENZDES_XML", "enzdes_design.xml")
#: designable positions for the targeted protocol, passed to %%des_pos%%
DES_POS = _os.environ.get("ENZDES_DES_POS", "")
CST = RO / "kemp_theozyme.cst"
PARAMS = RO / "6NT.params"

BASE_IDX, OXY_IDX, CHAIN = 102, 46, "A"
LIG_RESNAME, LIG_CHAIN, LIG_NUM = "6NT", "B", 1


ATOM_MAP = RO / "atom_map.json"


def _ours_to_rosetta() -> dict[str, str]:
    """Our ligand atom names -> the names molfile_to_params assigned.

    We do NOT rename the params. Rosetta's .params is read from fixed-width
    columns and two attempts to rewrite it in place went wrong in instructive
    ways: a split()/join() collapsed the columns ("add_bond: atom C5 and/or
    NO1 don't exist"), and a padded 4-char substitution then matched INSIDE the
    atom-TYPE field, turning `ONH2` into `ONH6` ("unrecognized atom_type_name").

    So the translation runs the other way. A PDB atom name sits in columns
    12-16 and a .cst names atoms in a plain token list - both are safe to edit,
    unlike a params file. The mapping itself is established by coordinate
    matching in make_6nt_params.py and stored in atom_map.json.

    Note the collision this protects against: our N1 (triazole, the Gln
    acceptor) is Rosetta's N2, while Rosetta's N1 is our NITRO nitrogen. An
    untranslated .cst would restrain the nitro group and report a plausible
    energy for the wrong geometry.
    """
    import json
    return json.loads(ATOM_MAP.read_text())["ours_to_rosetta"]


def translate_cst(src: Path, dst: Path, o2r: dict[str, str]) -> None:
    """Rewrite the ligand ATOM_MAP lines into Rosetta's atom names."""
    out = []
    for ln in src.read_text().splitlines():
        s = ln.strip()
        if s.startswith("TEMPLATE::") and "atom_name:" in s and ln.count(":") >= 2:
            head, names = ln.split("atom_name:", 1)
            toks = names.split()
            # only the LIGAND block is remapped; protein atom names are standard
            if all(t in o2r for t in toks):
                toks = [o2r[t] for t in toks]
                ln = f"{head}atom_name: {' '.join(toks)}"
        out.append(ln)
    dst.write_text("\n".join(out) + "\n")


def rename_ligand_atoms(lines: list[str], o2r: dict[str, str]) -> list[str]:
    """Rename 6NT atoms in-place, touching only PDB columns 12-16."""
    out = []
    for ln in lines:
        if ln[:6] in ("HETATM", "ATOM  ") and ln[17:20].strip() == LIG_RESNAME:
            name = ln[12:16].strip()
            if name in o2r:
                ln = ln[:12] + f"{o2r[name]:<4}" + ln[16:]
        out.append(ln)
    return out


def write_with_remarks(src: Path, dst: Path, o2r: dict[str, str]) -> None:
    """Prepend the REMARK 666 MATCH lines enzdes needs, then the coordinates."""
    lines = [l for l in src.read_text().splitlines()
             if not l.startswith("REMARK 666")]
    lines = rename_ligand_atoms(lines, o2r)
    remarks = [
        # block 1 -> the catalytic Asp
        f"REMARK 666 MATCH TEMPLATE {LIG_CHAIN} {LIG_RESNAME}  {LIG_NUM:3d} "
        f"MATCH MOTIF {CHAIN} ASP  {BASE_IDX:3d}  1  1",
        # block 2 -> the oxyanion Gln
        f"REMARK 666 MATCH TEMPLATE {LIG_CHAIN} {LIG_RESNAME}  {LIG_NUM:3d} "
        f"MATCH MOTIF {CHAIN} GLN  {OXY_IDX:3d}  2  1",
    ]
    dst.write_text("\n".join(remarks + lines) + "\n")


def main() -> int:
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    src = Path(sys.argv[1])
    nstruct = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    label = sys.argv[3] if len(sys.argv) > 3 else src.parent.name

    for p in (BIN, DB, XML, CST, PARAMS, ATOM_MAP):
        if not p.exists():
            raise SystemExit(f"missing prerequisite: {p}")
    if not src.exists():
        raise SystemExit(f"missing input: {src}")

    work = ROOT / "outputs/rosetta/runs" / label
    work.mkdir(parents=True, exist_ok=True)
    o2r = _ours_to_rosetta()
    cst_r = work / "theozyme_rosetta.cst"
    translate_cst(CST, cst_r, o2r)
    inp = work / "input.pdb"
    write_with_remarks(src, inp, o2r)
    print(f"  ligand atoms translated to Rosetta names "
          f"(our N1 -> {o2r['N1']}, N2 -> {o2r['N2']}, N3 -> {o2r['N3']})")
    print(f"{label}: {src}")
    print(f"  REMARK 666 -> block1=ASP{BASE_IDX}, block2=GLN{OXY_IDX}, "
          f"ligand {LIG_RESNAME} {LIG_CHAIN}{LIG_NUM}")

    cmd = [
        str(BIN),
        "-database", str(DB),
        "-s", str(inp),
        "-parser:protocol", str(XML),
        "-extra_res_fa", str(PARAMS),
        "-enzdes:cstfile", str(cst_r),
        "-run:preserve_header",          # REMARK 666 must survive parsing
        "-nstruct", str(nstruct),
        "-out:path:all", str(work),
        "-out:file:scorefile", "score.sc",
        "-overwrite",
        "-ex1", "-ex2",
        "-linmem_ig", "10",
        "-mute", "core.conformation core.scoring protocols.jd2",
    ]
    if DES_POS and ("targeted" in XML.name or "ligfix" in XML.name):
        cmd += ["-parser:script_vars", f"des_pos={DES_POS}"]
    log = work / "rosetta.log"
    print(f"  protocol={XML.name}"
          + (f"  des_pos={DES_POS}" if DES_POS and ("targeted" in XML.name or "ligfix" in XML.name) else ""))
    print(f"  nstruct={nstruct}; log -> {log}")
    t0 = time.time()
    with log.open("w") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                           timeout=10800)
    dt = (time.time() - t0) / 60
    print(f"  rc={p.returncode}  wall {dt:.1f} min")

    if p.returncode != 0:
        tail = log.read_text(errors="replace").splitlines()[-25:]
        print("  FAILED. log tail:")
        for l in tail:
            print("    " + l[:160])
        return 3

    sc = work / "score.sc"
    if sc.exists():
        rows = [l.split() for l in sc.read_text().splitlines()
                if l.startswith("SCORE:")]
        if len(rows) > 1:
            hdr = rows[0]
            print(f"\n  {len(rows)-1} design(s):")
            want = [h for h in ("total_score", "allcst", "lsasa", "description")
                    if h in hdr]
            idx = [hdr.index(h) for h in want]
            print("    " + "  ".join(f"{h:>14}" for h in want))
            for r in sorted(rows[1:], key=lambda r: float(r[hdr.index("allcst")])
                            if "allcst" in hdr else 0.0):
                print("    " + "  ".join(f"{r[i]:>14}" for i in idx))
            print("\n  LOWER allcst = constraints better satisfied. These are "
                  "PROPOSALS;\n  refold each with scripts/score_finalist.py "
                  "before believing anything.")
    pdbs = sorted(work.glob("*.pdb"))
    print(f"\n  {len(pdbs)} PDB(s) in {work}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
