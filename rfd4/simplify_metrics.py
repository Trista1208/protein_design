#!/usr/bin/env python3
"""Write readable TSV summaries alongside an RFD4 inference_metrics directory."""
import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metrics_dir", type=Path, help="Directory containing the two metrics CSV files")
    args = parser.parse_args()
    source = args.metrics_dir / "aggregated_metrics.csv"
    target = args.metrics_dir / "aggregated_metrics_aligned.txt"
    with source.open(newline="", encoding="utf-8") as inp:
        rows = list(csv.DictReader(inp))
    fields = ("dataset", "metric", "reduction", "value", "n")
    table = [[Path(row[f]).name if f == "dataset" else row[f] for f in fields] for row in rows]
    widths = [max(len(f), *(len(row[i]) for row in table)) if table else len(f) for i, f in enumerate(fields)]
    with target.open("x", encoding="utf-8") as out:
        out.write("Aggregate metrics across generated structures\n\n")
        for row in [fields, *table]:
            out.write("  ".join(value.rjust(widths[i]) if i >= 3 else value.ljust(widths[i]) for i, value in enumerate(row)).rstrip() + "\n")
    print(f"Wrote {target}")

    source = args.metrics_dir / "per_example_metrics.csv"
    target = args.metrics_dir / "per_example_metrics_aligned.txt"
    prefixes = ("amino_acid_composition.", "bond_length_metrics.", "clash_metrics.", "sequence_adherence.", "sequence_probability_adherence.", "sse_distribution.", "design_motif_", "design_unindexed_motif_", "unindexed_insertion_mapping.")
    with source.open(newline="", encoding="utf-8") as inp, target.open("x", encoding="utf-8") as out:
        reader = csv.DictReader(inp)
        fields = [f for f in reader.fieldnames if f == "pre_oracle_pass" or f.startswith(prefixes)]
        rows = list(reader)
        width = max(map(len, fields), default=0)
        value_width = max((len(row[f]) for row in rows for f in fields), default=0)
        for row in rows:
            structure = Path(row["designed_structure_path"])
            out.write(f"Example: {row['example_id']}\nStructure: {structure.name}\n")
            out.write(f"Location relative to run: inference_samples/{structure.parent.name}/{structure.name}\n")
            for field in fields:
                out.write(f"  {field:<{width}}  {row[field]:>{value_width}}\n")
            out.write("\n")
    print(f"Wrote {target}")


if __name__ == "__main__":
    main()
