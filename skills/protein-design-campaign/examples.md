# Example campaign: text-guided binder edit

Goal: locally edit an existing binder to improve on-target affinity while preserving scaffold, then Pareto-filter stability vs solubility.

## Stack

- Compile: PROTEAN two-stage (this skill pack `protean-editing`)
- Backbone: RFD4-Proteina partial diffusion + `C_CRD` on protected scaffold, `C_HOT` on target, `C_CRDGRP`+`C_SEQ` on receptor (`rfd4-proteina`)
- Sequence: RFD4 inverse folding ∪ ProteinMPNN on `M_edit` only (`inverse-folding-classics`)
- Multi-objective: MoMST `hydrophobic,surface_expose,plddt` or explicit Pareto (`momst-multiobjective`)
- Gates: G1/G3 then G2 then interface then developability (`protein-design-filters`)

Do not use Ptarmigan-1 unless the product is a small molecule.

## Conditioning

- `M_edit`: CDR-like loops / designed interface loops (not the whole chain)
- `M_relax`: 1–2 secondary-structure elements adjacent to the interface
- `M_prot`: core + distant helices; identities frozen; backbone RMSD contract 0.1–0.5 Å
- Instruction `u`: “Tighten packing against hotspots X,Y,Z; do not move the scaffold core; do not introduce unpaired Cys.”

YAML sketch: keep target as in the PD-L1 recipe; set binder contig to the existing designed chain; `generation.args.partial_trajectory_fraction_steps=0.25–0.35`; optional search with early H-bond/hydrophobic contacts and late `oracle.plddt_coi` + ipAE.

## Gates

1. G1: more hotspot atom contacts / lower ipAE vs parent; parent is the baseline, not zero.
2. G3: protected RMSD; reject even if affinity proxies look better.
3. G2: MPNN×4 + ESMFold/AF3 scRMSD < 2 Å (interface tighter).
4. Hard reject: unpaired Cys, huge hydrophobic patch, glycosylation sequon.
5. Pareto: keep high-affinity-proxy, balanced, and high-solubility representatives per Foldseek cluster.
6. Off-target AF3 against the closest paralogue.

## Failure modes

Scope leakage (PROTEAN bottleneck), refold to a different loop, aggregation from hydrophobic gain, nonspecific binding, expression loss. Expect most G1 winners to die at G2 or developability.

## Iterate

After BLI/SPR: add failures to MoMST recency memory; tighten `C_CRD` on leaking regions; raise hotspot weight in search; do not widen `M_edit` on the first failure.
