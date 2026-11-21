#!/usr/bin/env bash
# Re-download RCSB FASTA + PDB for identity / TM / theozyme extract.
set -euo pipefail
DIR="$(cd "$(dirname "$0")/../references" && pwd)"
mkdir -p "$DIR"
cd "$DIR"
for id in 4BS0 3NYD 2RKX 9HVB 7EPO 1MG5; do
  curl -fsSL -k "https://www.rcsb.org/fasta/entry/${id}" -o "${id}.fasta"
done
for id in 4BS0 3NYD 2RKX 9HVB 7EPO; do
  curl -fsSL -k "https://files.rcsb.org/download/${id}.pdb" -o "${id}.pdb"
done
cat 4BS0.fasta 3NYD.fasta 2RKX.fasta 9HVB.fasta 7EPO.fasta 1MG5.fasta > published_kemp.fasta
echo "references in $DIR"
