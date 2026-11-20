---
name: rfd4-proteina
description: >-
  Runs and specifies RFD4-Proteina (Baker Lab × NVIDIA / RosettaCommons) for
  de novo monomers, protein/DNA binders, antibodies, enzymes, motif scaffolding,
  glycans, partial diffusion, and reward-guided search. Use when the user mentions
  RFD4, Proteina, annotated CIF, YAML condition specs, rfd4-generate, hotspots,
  or recycling/thinking presets.
---

# RFD4-Proteina

Canonical docs: https://rfd4-proteina.org/  
Code: https://github.com/RosettaCommons/RFD4-Proteina  
Paper: Corley et al., 2026. *Biomolecular design through latent reasoning with RFD4-Proteina*. doi:10.64898/2026.11.12.766037

If a live RFD4 checkout is present, prefer the repo's `skills/` (`rfd4-proteina-setup`, `generate`, `search`, `evaluate`, `campaign`) for CLI details. This skill is the scientific operator layer.

Recipes: [yaml-recipes.md](yaml-recipes.md)

## What it is

Flow-matching foundation model in a shared 32-D chemical/conformational latent space (standard + modified amino acids, sugars, DNA, RNA), trained with force-field energies.

- **Intrinsic reasoning**: recurrent trunk recycling (thinking presets).
- **Extrinsic reasoning**: reward-guided beam search (geometry early, AF3/energy late).
- Sequence–structure co-design; also inverse folding, side-chain packing, conformational sampling, refinement without extra training.

Platform: Linux x86_64, Python 3.12, NVIDIA GPU (CUDA 12/13). macOS/Windows unsupported. Browser: Hugging Face Spaces (SecureDNA screened).

## Install (default)

```bash
git clone --depth 1 https://github.com/RosettaCommons/RFD4-Proteina.git
cd RFD4-Proteina
make setup
source .venv/bin/activate
```

Or: `uv pip install "rfd4-proteina[all-cuda13]"` (CUDA 12: `all-cuda12`). Skip search oracles with `bash setup.sh --skip-search-envs`.

## Specify a design

Three equivalent inputs. Prefer YAML condition specs unless you need Python (two-body conditions, pure de novo with no anchor, programmatic selections).

Workflow: load structure → `contig` subset/insert expandables → `conditions` on the contig result. Anything omitted is **not** conditioned. Trust the printed summary.

Relative `input` paths resolve against the **spec file directory**.

### Contig

- Preserve: `A17-131`
- New chain / chain break: `/0`
- Expandable binder length: `50-120` or named `BINDER:60-90`
- At least one preserved segment is required in a spec. Pure monomers use benchmark CIFs with `input_as_annotated_cif: true`.

`res_id` in contig and selections is the **original** mmCIF `label_seq_id`. Expandables have no original `res_id`; select by `res_name`.

### Selections

```
sel('chain/res_name/res_id/atom_name')   # omitted trailing fields = *
```

`sel('B/*/N')` is wrong (parses `res_id=N`). Use `sel('B/*/*/N')`. Query syntax: `"sel('B') and is_protein_backbone"`.

### Conditions that matter

| Annotation | Typical use |
|---|---|
| `C_CRDGRP` | Relative geometry of a rigid target (preferred for PPI). One positive group ID. Do **not** also set `C_CRD` on those atoms. |
| `C_SEQ` | Reveal target (or motif) sequence. |
| `C_HOT` | Hotspot atoms on the target. |
| `C_CRD` | Absolute coordinates (motifs, protected scaffold, backbone-only redesign). |
| `C_IDX` | Index groups for motif scaffolding. |
| `C_HB` | H-bond donor/acceptor conditioning. |
| `C_SEQ_PROB_CYS` | `0.0` is “no cysteine”, distinct from unspecified. |

PPI default: `C_CRDGRP` + `C_SEQ` on target, optional `C_HOT`. Backbone-only `C_CRD` on a binder when inverse-folding: `sel('B') and is_protein_backbone`.

Inspect before GPU:

```bash
scripts/write_condition_spec_cifs.py binder.yaml --out-dir ./cifs
```

## Generate

```bash
rfd4-generate inputs=binder.yaml num_replicates=8 batch_size=8 \
  recycle_strategy=medium/front_loaded out_dir=./runs/binder \
  metrics_tags='[ppi_refolding,hotspot_conditioning]'
```

- First smoke: `inference_experiment=quickstart` and `'~metrics'` (quote `~` in zsh). Add `compile_lengths=null` to skip compile warmup.
- `recycle_strategy`: binders `medium/front_loaded`; enzymes/SM `medium/midpeak`; VHH `high/front_loaded`. Low/medium/high thinking trades throughput vs clustered AF3 pass rate.
- Outputs: `out_dir/inference_samples/`. Partial-diffuse **generated** CIFs, not oracle refolds.
- Dry-run config: `rfd4-generate inputs=binder.yaml --cfg job`

Hugging Face / Python:

```python
import rfd4_proteina as rfd4
model = rfd4.Proteina.from_pretrained("nvidia/RFD4-Proteina")
rfd4.generate(model, overrides=["inference_experiment=quickstart", "inputs=design.yaml", "out_dir=runs/hf"])
```

Diffusers takes prepared tensors from an annotated CIF, not YAML. Convert first.

## Partial diffusion and inverse folding

- Refine: `generation.args.partial_trajectory_fraction_steps=0.3` (last 30%; larger = more drift).
- Sequence only: `generation=inverse_folding recycle_strategy=medium/flat`.

## Reward-guided search

Use when single-shot pass rate is too low or constraints are geometric/energetic. Checkpoints are integration steps; each lists its own metrics and rewards (nothing carries over).

- Early (~20–40): interface area, H-bonds, hydrophobic contacts (`use_clean_prediction: true` to avoid full AF3).
- Late (~60): oracle pLDDT / ipAE / energy.
- `beam_width`, `n_branch`; empty reward keeps input order.
- Continuation criteria stop hopeless beams. Foldseek novelty against a passer archive avoids duplicate folds.
- Examples: `inference_experiment=inference_beam_example` (AF3), `sm_binder_hbond_search` (structure-only).

Search oracles need the AF3/OpenFold3 (or AF2) env from setup.

## Task routing

| Task | Recycle | Metrics tags |
|---|---|---|
| Monomer | default | `unconditional_refolding` |
| Protein binder | `medium/front_loaded` | `ppi_refolding`, `hotspot_conditioning` |
| Enzyme / motif | `medium/midpeak` | `ame_refolding`, `motif_scaffold_refolding` |
| Small-molecule binder | `medium/midpeak` | `sm_binder_refolding`, `interchain_interactions` |
| VHH | `high/front_loaded` via `examples/antibody_design` | campaign config |
| Partial refine | keep parent | usually `'~metrics'` then evaluate |

## Wet-lab context (do not overclaim)

Reported experimental campaigns include Ca-dependent calmodulin binders, peptide–MHC specificity, disordered peptide binders (YAP1/Tau), TrkA species/paralog specificity, glycan shielding, phosphopeptide binders, improved phosphotriesterases, and a lenacapavir-step ligase. In silico, clustered AF3 pass rate + target coverage beat compared methods at higher throughput; thinking longer raises pass rate and drops designs/GPU-hour.

## Limitations

- Specs cannot express two-body conditions or fully de novo designs with no preserved residue.
- Scope leakage still happens: validate protected RMSD, not just the YAML.
- Metrics that refold need the oracle environment; otherwise generate with `'~metrics'` and score later via `rfd4-evaluate`.
- Always apply the filter waterfall from `protein-design-filters` before ordering DNA.
