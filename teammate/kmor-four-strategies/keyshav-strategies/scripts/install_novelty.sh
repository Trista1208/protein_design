#!/bin/bash
# Install the novelty-gate toolchain and build the reference corpus.
# Dedicated env so the verified `mm` env is not disturbed.
# mmseqs2 -> sequence identity ; foldseek -> structure TM-score (both directions)
set -u
CONDA=/home/mandalorian/anaconda3/bin/conda
ENVS=/home/mandalorian/anaconda3/envs
D=/home/mandalorian/Hackathons/Enzymes-BSSE2026
CORPUS=$D/inputs/corpus
CF="--override-channels -c conda-forge -c bioconda"
H(){ echo; echo "########## $* ##########"; date '+%F %T'; }

H "1. create novelty env (mmseqs2 + foldseek)"
if [ -x "$ENVS/novelty/bin/mmseqs" ] && [ -x "$ENVS/novelty/bin/foldseek" ]; then
  echo "  already present"
else
  $CONDA create -y -n novelty $CF python=3.12 mmseqs2 foldseek 2>&1 | tail -6
fi
echo "  mmseqs:   $($ENVS/novelty/bin/mmseqs version 2>&1 | head -1)"
echo "  foldseek: $($ENVS/novelty/bin/foldseek version 2>&1 | head -1)"

H "2. fetch the novelty corpus from RCSB"
mkdir -p "$CORPUS/structures" "$CORPUS/seqs"
# handover 4.3.8: a real corpus, not just 4BS0/6TU6
IDS="6TU6 4BS0 7K4Z 7K4Q 3NYD 7EPO 7EPN 9HVB 9HVG 9HVH 1A53"
for id in $IDS; do
  f="$CORPUS/structures/${id}.cif"
  if [ -s "$f" ]; then echo "  have $id"; continue; fi
  if curl -fsSL "https://files.rcsb.org/download/${id}.cif" -o "$f"; then
    echo "  fetched $id ($(du -h "$f" | cut -f1))"
  else
    echo "  FAILED $id"; rm -f "$f"
  fi
done
# sequences for the identity gate
for id in $IDS; do
  f="$CORPUS/seqs/${id}.fasta"
  [ -s "$f" ] && continue
  curl -fsSL "https://www.rcsb.org/fasta/entry/${id}" -o "$f" 2>/dev/null || rm -f "$f"
done
cat "$CORPUS"/seqs/*.fasta > "$CORPUS/published.fasta" 2>/dev/null
echo "  corpus: $(ls "$CORPUS/structures" | wc -l) structures, \
$(grep -c '^>' "$CORPUS/published.fasta" 2>/dev/null || echo 0) sequences"

H "3. KE07/KE70 + HG-series by sequence (no single PDB id applies)"
# Rothlisberger 2008 Baker-lab designs; fetched from UniProt where available.
for acc in P0DP24; do : ; done   # placeholder: KE designs are not in UniProt
cat >> "$CORPUS/published.fasta" <<'FA'
>KE07_Rothlisberger2008_note
MNLKEKLIAPLLSGNVSAEEALEVLKELGNVDVSVLDKMIAQSQQLNLSAEELKQLLAQLQ
>KE70_Rothlisberger2008_note
MSLKEQLLASLLAGDVSAEQALDLLKQLGNVDLSVLEKMIAQSQSLNLSAEDLKALLAKLQ
FA
echo "  NOTE: the KE07/KE70 entries above are PLACEHOLDERS, not real sequences."
echo "  They exist so the gate wiring is exercised. Replace with the real"
echo "  sequences from Rothlisberger 2008 SI before trusting an identity verdict."

H "4. build mmseqs + foldseek databases"
cd "$CORPUS"
if [ -s published.fasta ]; then
  $ENVS/novelty/bin/mmseqs createdb published.fasta seqdb 2>&1 | tail -2
  echo "  seqdb built"
fi
if [ "$(ls structures/*.cif 2>/dev/null | wc -l)" -gt 0 ]; then
  $ENVS/novelty/bin/foldseek createdb structures structdb 2>&1 | tail -2
  echo "  structdb built"
fi

H "5. smoke test both gates against a real query"
Q=$D/outputs/runs/m0/theozyme_clean.pdb
if [ -s "$Q" ] && [ -f structdb ]; then
  $ENVS/novelty/bin/foldseek easy-search "$Q" structdb tm_hits.tsv tmp_fs \
     --format-output "query,target,alntmscore,qtmscore,ttmscore,fident,evalue" \
     --exhaustive-search 1 2>&1 | tail -3
  echo "  --- TM hits (both normalisations):"
  head -5 tm_hits.tsv 2>/dev/null | sed 's/^/    /'
fi
echo
echo "### INSTALL_NOVELTY DONE ###"
