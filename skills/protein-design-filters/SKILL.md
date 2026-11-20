---
name: protein-design-filters
description: >-
  Applies experimentally realistic in-silico filter waterfalls for protein
  designs: pLDDT, PAE/ipAE, scRMSD, interface metrics, motif RMSD, energy,
  expression/solubility proxies, diversity, and PROTEAN G1/G2/G3. Use when
  ranking designs, setting pass criteria, or deciding what to order for wet lab.
---

# Protein design filters

Hard gates first, then rank. Do not average incompatible scores into one number until geometry and refolding pass.

## Waterfall (default)

```
- [ ] 0. File hygiene: one design per ID; keep generated CIF vs oracle refold distinct
- [ ] 1. Scope / motif geometry (fail closed)
- [ ] 2. Designability / refold
- [ ] 3. Interface or active site
- [ ] 4. Developability / expression proxies
- [ ] 5. Specificity / negative states
- [ ] 6. Diversity clustering
- [ ] 7. Pareto shortlist for DNA
```

Thresholds are **starting points**. Recalibrate on the first experimental round.

## 1. Geometry and conservation

| Check | Typical pass | Failure mode |
|---|---|---|
| Protected scaffold Cα RMSD after superposition | ≤ 0.5 Å campaign; PROTEAN contract 0.1 Å backbone | Scope leakage |
| Motif / catalytic site backbone RMSD | ≤ 1.0 Å | Dead enzyme |
| Ligand/substrate pose vs design | visual + clash-free; task PAE if used | Wrong chemistry |
| Cyclization C–N | 1.20–1.50 Å, exactly one new bond | Fake cyclization |
| Secondary-structure / distance / RSA edits | PROTEAN G1 table | Instruction not realized |

RFD4 enzyme campaigns also used protein–substrate min PAE < 1.5 Å with motif RMSD < 1 Å as a **common** in silico bar (strict; low pass fraction).

## 2. Designability (refold)

Inverse fold (ProteinMPNN or RFD4 `inverse_folding`) then predict:

| Oracle | Use |
|---|---|
| AlphaFold3 / OpenFold3 | Complexes, ligands, DNA, RFD4 search rewards |
| AlphaFold2 / ColabDesign | Classic monomer/PPI self-consistency |
| ESMFold | Fast monomers; PROTEAN G2 |
| Boltz | Co-folding / complexes when AF3 unavailable |

| Metric | Monomer (start) | Binder / complex (start) |
|---|---|---|
| pLDDT (designed chain) | ≥ 80–90 | ≥ 80 |
| PAE / ipAE (interface) | — | low ipAE; min PAE at hotspot |
| scRMSD design vs refold | < 2 Å (whole); motif < 1 Å | designed chain < 2 Å; interface residues tight |
| iPTM / ipTM | — | prefer high; don't use alone |
| Clustered pass | unique folds only | repeating one success is worthless |

PROTEAN G2: ≥1 of 4 MPNN+ESMFold sequences meets frozen RMSD cutoffs. Failure is a **sampling miss**, not proof of impossibility — but do not order G2 fails.

## 3. Interface / active site

- Hotspot contact: designed chain touches intended `C_HOT` residues.
- Buried SASA at interface; shape complementarity; H-bond / hydrophobic counts (RFD4 `interchain_interactions`).
- Unwanted target-like sequence on the binder (why PPI uses `C_CRDGRP` not full `C_CRD` + no binder sequence copy).
- Catalysis: first-shell identities fixed; alanine-scan in silico is not a substitute for experiment but flags missing geometry.
- Clash / buried unsatisfied polar / unpaired Cys at the interface.

## 4. Developability proxies (not assays)

Reject or flag:

- Unpaired cysteines (unless disulfide designed)
- N-X-S/T glycosylation sequons if not wanted
- Long hydrophobic stretches / high SAP / large exposed hydrophobic patches
- Extreme pI vs expression host; poly-R/K stretches
- Proline in the middle of a designed helix unless intended
- Rare codons only at DNA-order time, not as a fold filter

MoMST-style metrics (`hydrophobic`, `surface_expose`, `globularity`, `plddt`) **rank** after these hard rejects.

## 5. Specificity

Design the negative class:

- Paralogues / species (TrkA vs TrkB/C; human vs mouse)
- Apo vs holo (Ca-free calmodulin)
- Adjacent PTM sites (H3 pS10 vs pT11)
- Peptide–MHC off-targets (1–2 residue lookalikes)

Filter: AF3/Boltz against off-targets; require worse ipAE / no hotspot contacts. Yeast/phage still required.

## 6. Diversity

Foldseek / TM-score cluster. Keep 1–3 per cluster. RFD4 `FoldseekNoveltyMetric` + passer archive during search. Wet-lab rounds of **20–96** diverse passers beat 500 near-duplicates.

## 7. What to send to the lab

Report for each candidate:

```
id, length, cluster, G1, G3, G2_scRMSD, pLDDT, ipAE, motif_RMSD,
hotspot_contacts, developability_flags, off_target_score, rank_on_pareto
```

Expected in silico → wet-lab drop: large. Quote analogous RFD4 hit rates when the task matches; otherwise say unknown. Plan a second round that updates rewards from failures (aggregation, no expression, nonspecific binding).

## RFD4 scoring hooks

- Inline: `metrics_tags='[ppi_refolding,...]'` needs AF3/OpenFold3 env.
- Skip during smoke tests: `'~metrics'`
- Post-hoc: `rfd4-evaluate`
- Search rewards ≠ final waterfall; still run full filters on saved beams.
