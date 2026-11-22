#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

PASS_COLS = ("ame.pass", "pre_oracle_pass")
PRIORITY_ASC = (
    "refolded_motif_allatom_rmsd.rmsd",
    "design_motif_allatom_rmsd.rmsd",
    "design_unindexed_motif_allatom_rmsd.rmsd",
    "refolded_backbone_rmsd.rmsd",
    "clash_metrics.frac_atoms_clashing_1p5A",
    "bond_length_metrics.bond_length.frac_unreasonable",
    "sequence_adherence.n_violations",
)
PRIORITY_DESC = (
    "oracle.overall_plddt",
    "oracle.iptm",
    "sequence_adherence.sequence_adherence",
)
EXCLUDE_COLS = {"example_id", "designed_structure_path", "oracle.refolded_structure_path"}


def number(row, field):
    try:
        value = float(row[field])
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("metrics_dir", type=Path)
    args = parser.parse_args()
    source = args.metrics_dir / "per_example_metrics.csv"
    target = args.metrics_dir / "structures_ranked.txt"

    with source.open(newline="", encoding="utf-8") as inp:
        reader = csv.DictReader(inp)
        fieldnames = list(reader.fieldnames or ())
        if "example_id" not in fieldnames:
            parser.error("Missing required column: example_id")
        rows_raw = list(reader)

    numeric_cols = [
        c for c in fieldnames
        if c not in EXCLUDE_COLS and any(number(r, c) is not None for r in rows_raw)
    ]

    pass_cols = [c for c in PASS_COLS if c in numeric_cols]
    asc_cols = [c for c in PRIORITY_ASC if c in numeric_cols]
    desc_cols = [c for c in PRIORITY_DESC if c in numeric_cols]
    other_cols = [c for c in numeric_cols if c not in set(pass_cols + asc_cols + desc_cols)]

    def sort_key(row):
        parts = []
        for col in pass_cols:
            val = number(row, col)
            parts.append(val != 1)
        for col in asc_cols:
            val = number(row, col)
            parts.extend((val is None, val if val is not None else 0.0))
        for col in desc_cols:
            val = number(row, col)
            parts.extend((val is None, -val if val is not None else 0.0))
        parts.append(row["example_id"])
        return tuple(parts)

    rows = sorted(rows_raw, key=sort_key)
    ordered_metrics = pass_cols + asc_cols + desc_cols + other_cols
    headers = ("rank", "example_id", *ordered_metrics)

    table = []
    for index, row in enumerate(rows, 1):
        formatted = []
        for col in ordered_metrics:
            val = number(row, col)
            formatted.append(f"{val:.4f}" if val is not None and not val.is_integer() else (row[col] or "NA"))
        table.append([str(index), row["example_id"], *formatted])

    widths = [max([len(header)] + [len(r[i]) for r in table]) for i, header in enumerate(headers)]
    with target.open("w", encoding="utf-8") as out:
        for r in [headers, *table]:
            out.write("  ".join(val.ljust(widths[i]) if i == 1 else val.rjust(widths[i]) for i, val in enumerate(r)).rstrip() + "\n")
    print(f"Wrote {target}")


if __name__ == "__main__":
    main()
