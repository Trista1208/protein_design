# RFD4 YAML recipes

Paths are relative to the spec file. Run from the RFD4-Proteina repo root unless `input` is absolute.

## Binder (PD-L1 pattern)

```yaml
pdl1_binder:
  input: docs/examples/inputs/5o45_cropped.cif.zst
  contig: "50-120, /0, A17-131"
  conditions:
    C_CRDGRP: {1: [{select: "sel('A/*/*/[N,CA,C,O]')"}]}
    C_SEQ: {true: [{select: "sel('A')"}]}
    C_HOT:
      true:
        - {select: "sel('A/*/56/[CG,OH]')"}
        - {select: "sel('A/*/115/[CG,SD]')"}
        - {select: "sel('A/*/123/[CD2,OH]')"}
```

```bash
rfd4-generate inputs=binder.yaml recycle_strategy=medium/front_loaded \
  num_replicates=8 batch_size=8 out_dir=./runs/binder \
  metrics_tags='[ppi_refolding,hotspot_conditioning]'
```

## Motif scaffold (named expandables)

```yaml
motif_scaffold:
  input: motif.cif
  contig: "A1-50, LINK:10-30, /0, BINDER:60-90"
  conditions:
    C_CRD: {true: [{select: "sel('A/*/[1, 50]')"}]}
    C_SEQ: {true: [{select: "chain_id == 'A' and res_name != 'LINK'"}]}
```

Benchmark enzymes: wrap M-CSA CIF with `input_as_annotated_cif: true`; `recycle_strategy=medium/midpeak`; `metrics_tags='[ame_refolding,motif_scaffold_refolding]'`.

## Binder sequence redesign (keep backbone)

```yaml
binder_redesign:
  input: complex.cif
  conditions:
    C_CRD: {true: [{select: "sel('B') and is_protein_backbone"}]}
    C_HOT: {true: [{select: "sel('A/*/[45, 46, 82]')"}]}
```

Or `generation=inverse_folding recycle_strategy=medium/flat`.

## Partial diffusion refine

```yaml
refine:
  input: runs/binder/inference_samples/<example>/<example>_generated.cif
  input_as_annotated_cif: true
```

```bash
rfd4-generate inputs=refine.yaml \
  generation.args.partial_trajectory_fraction_steps=0.3 \
  num_replicates=8 out_dir=./runs/refine
```

## Monomer

```yaml
monomer:
  input: data/benchmarks/monomer/monomer-short.cif.zst
  input_as_annotated_cif: true
```

```bash
rfd4-generate inputs=monomer.yaml num_replicates=16 batch_size=8 \
  out_dir=./runs/monomer metrics_tags='[unconditional_refolding]'
```

## Search sketch (early contacts, late oracle)

```yaml
search:
  beam_width: 2
  n_branch: 2
  checkpoints:
    "20":
      n_branch: 4
      metrics: [interchain_interactions]
      use_clean_prediction: true
      rewards:
        interchain_interactions.hbond_atom_pair_count_coi: {weight: 1.0}
    "40":
      metrics: [interchain_interactions]
      rewards:
        interchain_interactions.hbond_atom_pair_count_coi:
          {weight: 1.0, clip: {max: 6.0}, normalize_by: 6.0}
        interchain_interactions.hydrophobic_atom_pair_count_coi: {weight: 0.5}
    "60":
      metrics: [oracle]
      rewards:
        oracle.plddt_coi: {weight: 1.0}
generation:
  args:
    nsteps: 60
```
