#!/usr/bin/env python3
"""TM-score vs published structures. Does not invent a TM-score.

Uses USalign or TMalign if installed. Exit 2 if neither binary exists.
Official fail: TM > 0.7 vs a published structure.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REF_DIR = HERE.parent / "references"
MAX_TM = 0.70
REF_PDBS = ["4BS0.pdb", "3NYD.pdb", "2RKX.pdb", "9HVB.pdb", "7EPO.pdb"]


def find_bin() -> str | None:
    for name in ("USalign", "TMalign", "tmalign"):
        p = shutil.which(name)
        if p:
            return p
    return None


def run_tm(binary: str, query: Path, ref: Path) -> float | None:
    proc = subprocess.run(
        [binary, str(query), str(ref)],
        check=False,
        capture_output=True,
        text=True,
    )
    text = proc.stdout + "\n" + proc.stderr
    # TM-score=0.12345 (if normalized by length of Chain_1)
    hits = re.findall(r"TM-score\s*=\s*([0-9.]+)", text)
    if not hits:
        return None
    return max(float(x) for x in hits)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("pdb", type=Path)
    p.add_argument("--max-tm", type=float, default=MAX_TM)
    args = p.parse_args()
    binary = find_bin()
    if not binary:
        print(
            "TM-score not computed: install USalign or TMalign and put it on PATH.\n"
            "  https://github.com/pylelab/USalign\n"
            "Do not treat RMSD or Foldseek bitscore as TM-score.",
            file=sys.stderr,
        )
        return 2

    worst = 0.0
    worst_ref = ""
    for name in REF_PDBS:
        ref = REF_DIR / name
        if not ref.exists():
            print(f"skip missing {ref}", file=sys.stderr)
            continue
        tm = run_tm(binary, args.pdb, ref)
        if tm is None:
            print(f"WARN could not parse TM for {name}", file=sys.stderr)
            continue
        print(f"TM {tm:.3f}  vs {name}")
        if tm > worst:
            worst, worst_ref = tm, name
    if worst > args.max_tm:
        print(f"FAIL TM {worst:.3f} vs {worst_ref} > {args.max_tm}")
        return 1
    print(f"PASS max TM {worst:.3f} vs {worst_ref or 'n/a'} ≤ {args.max_tm}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
