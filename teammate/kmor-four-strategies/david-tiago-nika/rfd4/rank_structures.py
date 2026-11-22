#!/usr/bin/env python3
"""Rank RFD4 structures using screening metrics and motif RMSD."""
import argparse
import csv
import math
from pathlib import Path

PASS = "pre_oracle_pass"
BOND = "bond_length_metrics.bond_length.frac_unreasonable"
CLASH = "clash_metrics.frac_atoms_clashing_1p5A"
RMSD = "design_unindexed_motif_allatom_rmsd.rmsd"
ADHERENCE = "sequence_adherence.sequence_adherence"
VIOLATIONS = "sequence_adherence.n_violations"
UNMAPPED = "unindexed_insertion_mapping.n_unmapped_tokens"
MISMATCHES = "unindexed_insertion_mapping.n_sequence_mismatches"
METRICS = (PASS, BOND, CLASH, RMSD, ADHERENCE, VIOLATIONS, UNMAPPED, MISMATCHES)


def number(row, field):
    try:
        value = float(row[field])
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def key(row):
    values = {field: number(row, field) for field in METRICS}
    # Missing values sort after known values, never as perfect scores.
    def ascending(field):
        value = values[field]
        return (value is None, value if value is not None else 0)
    passed = values[PASS]
    adherence = values[ADHERENCE]
    return ((passed != 1), *ascending(BOND), *ascending(CLASH), *ascending(RMSD),
            adherence is None, -adherence if adherence is not None else 0,
            *ascending(VIOLATIONS), *ascending(UNMAPPED), *ascending(MISMATCHES),
            row["example_id"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metrics_dir", type=Path, help="Directory containing per_example_metrics.csv")
    args = parser.parse_args()
    source = args.metrics_dir / "per_example_metrics.csv"
    target = args.metrics_dir / "structures_ranked.txt"
    with source.open(newline="", encoding="utf-8") as inp:
        reader = csv.DictReader(inp)
        required = {"example_id", "designed_structure_path", *METRICS}
        missing = required - set(reader.fieldnames or ())
        if missing:
            parser.error("Missing columns: " + ", ".join(sorted(missing)))
        rows = sorted(reader, key=key)

    headers = ("rank", "example_id", "pass", "bond", "clash", "motif_rmsd", "adherence", "violations", "unmapped", "mismatches", "structure")
    table = []
    for index, row in enumerate(rows, 1):
        structure = Path(row["designed_structure_path"])
        table.append([str(index), row["example_id"], *(row[field] or "NA" for field in METRICS),
                      f"inference_samples/{structure.parent.name}/{structure.name}"])
    widths = [max([len(header)] + [len(row[i]) for row in table]) for i, header in enumerate(headers)]
    with target.open("x", encoding="utf-8") as out:
        out.write("Ranking: pass first, then lower bond, clash and motif RMSD; higher adherence; lower violations, unmapped and mismatches.\n")
        out.write("Missing values rank after known values. Structure paths are relative to the run directory.\n\n")
        for row in [headers, *table]:
            out.write("  ".join(value.ljust(widths[i]) if i in (1, 10) else value.rjust(widths[i]) for i, value in enumerate(row)).rstrip() + "\n")
    print(f"Wrote {target}")


if __name__ == "__main__":
    main()
