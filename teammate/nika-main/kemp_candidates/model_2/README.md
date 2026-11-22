# Kemp candidate: model 2

One Boltz-2 seed of a 171-residue Kemp eliminase. This is a candidate to inspect, not a five-seed shortlist design.

## Files

- `kemp_model_2.cif` — Boltz model 2, the only seed with both catalytic contacts under 5 Å
- `kemp_model_2.fasta` — the sequence in that model
- `five_seed_contacts.json` — all five seeds of this same sequence

## What passed

- Length 171, under the 260 limit.
- Identity about 25% to the closest published sequence in our local reference set (`1MG5`), under the 70% limit.
- Asp41 and Gln65 are fixed. In this seed, Asp OD2–6NT N3 is 4.71 Å and Gln NE2–6NT N1 is 3.89 Å.

## What did not

- Seeds 0, 1, 3, and 4 of the same sequence miss at least one contact. The Asp distance ranges from 4.71 Å to 9.13 Å.
- Four new sequences designed on this structure, with only Asp41 and Gln65 fixed, also lost the Asp contact.
- TM-score against the published Kemp structures has not been run.
- There is one cysteine, at position 21.

## How it was made

RFD3 scaffolded the Theozyme_withORI motif (seed 1, length 171). LigandMPNN redesigned every residue except Asp41 and Gln65. Five Boltz-2 seeds were run on that sequence. Model 2 is the seed that kept both contacts under 5 Å.
