"""Minimal PDB parser. No hydrogen alt-B, chain A occupancy as in the file."""

from __future__ import annotations

from pathlib import Path
from math import sqrt
from typing import Iterator


def parse_atoms(path: Path | str, chain: str | None = None, alts: tuple[str, ...] = (" ", "A")) -> list[dict]:
    atoms: list[dict] = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        alt = line[16]
        if alt not in alts:
            continue
        rec = {
            "record": line[:6].strip(),
            "an": line[12:16].strip(),
            "rn": line[17:20].strip(),
            "ch": line[21],
            "ri": int(line[22:26]),
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
            "occ": float(line[54:60]) if len(line) >= 60 else 1.0,
            "line": line[:66],
        }
        rec["het"] = rec["record"] == "HETATM"
        rec["heavy"] = not rec["an"].startswith("H") and rec["an"] != "H"
        if chain is not None and rec["ch"] != chain:
            continue
        atoms.append(rec)
    return atoms


def dist(a: dict, b: dict) -> float:
    return sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def find(atoms: list[dict], *, ri: int | None = None, rn: str | None = None, an: str | None = None, ch: str | None = None, het: bool | None = None) -> list[dict]:
    out = []
    for a in atoms:
        if ri is not None and a["ri"] != ri:
            continue
        if rn is not None and a["rn"] != rn:
            continue
        if an is not None and a["an"] != an:
            continue
        if ch is not None and a["ch"] != ch:
            continue
        if het is not None and a["het"] != het:
            continue
        out.append(a)
    return out


def one(atoms: list[dict], **kw) -> dict:
    hits = find(atoms, **kw)
    if len(hits) != 1:
        raise KeyError(f"expected 1 atom for {kw}, got {len(hits)}")
    return hits[0]


def write_pdb(path: Path | str, atoms: list[dict], remarks: list[str] | None = None) -> None:
    lines = ["REMARK   1 theozyme extract — not a full protein"]
    for r in remarks or []:
        lines.append("REMARK   2 " + r[:66])
    serial = 1
    for a in atoms:
        an = a["an"]
        lines.append(
            f"{'HETATM' if a['het'] else 'ATOM  '}"
            f"{serial:5d} "
            f"{an:>4s}"
            f" "
            f"{a['rn']:>3s} "
            f"{a['ch']}"
            f"{a['ri']:4d}"
            f"    "
            f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
            f"{a['occ']:6.2f}{0.0:6.2f}"
        )
        serial += 1
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n")
