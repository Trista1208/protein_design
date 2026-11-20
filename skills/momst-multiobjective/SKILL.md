---
name: momst-multiobjective
description: >-
  Applies MoMST memory-aware test-time scaling for multi-objective protein
  sequence design (stability, solubility, globularity, pLDDT, hydrophobicity)
  without retraining. Use when the user mentions Pareto fronts, conflicting
  rewards, self-contrastive residue preferences, memory banks, or MoMST.
---

# MoMST multi-objective design

Paper: Yang et al., ICML 2026. *Multi-Objective Protein Design via Memory-Aware Test-Time Scaling in Diffusion Models*. https://proceedings.mlr.press/v306/yang26du.html  
Code: https://github.com/MingYangi/MoMST

## Problem

Therapeutic and industrial sequences must hit several scores at once (fold confidence, compactness, hydrophobicity, surface exposure). Retraining a diffusion model per objective mix is too expensive. Naive test-time methods repeat old errors, overfit to rare successes, and collapse to one objective.

## Method (plain language)

Iterate **noise → memory-guided denoise** on a pretrained sequence diffusion model (codebase builds on RERD / EvoDiff / OpenFold / ESMFold):

1. **Memory bank M** (hybrid, 50/50):
   - **Elite / static tier**: best samples since the start (long-term high-reward patterns).
   - **Recency / sliding window**: latest successes **and failures** so the policy sees current local landscape.
2. **Self-contrastive residue preferences**: from both winners and losers in M, extract **per-residue** relative preferences (not a single sequence-level reward). These reshape the pretrained distribution `p_pre` into a memory-aware `p_M` without destroying sequence naturalness.
3. **Inference-time Pareto alignment**: after a trajectory, select balanced proteins (rank-normalized weighted reward, Pareto filtering, worst-objective robust scores) or re-noise them for another cycle.

This is **training-free** test-time compute: more iterations / candidates instead of new weights.

## When to use

- Sequence-level trade-offs on a foldable backbone (or de novo sequence models).
- Need to change the objective mix mid-campaign (add solubility, drop Cys, raise pLDDT).
- Backbone is already conditioned (RFD4/RFdiffusion); you are ranking or refining sequences.

Do **not** use MoMST as a substitute for 3D hotspot/motif constraints. Geometry first, Pareto second.

## Practical usage

```bash
conda create -n MoMST python=3.9 && conda activate MoMST
pip install torch torchvision torchaudio
pip install -r requirements.txt
python download_model_data.py   # example PDBs into ./datasets/AlphaFoldPDB/
```

Requires PyTorch + pyrosetta. Decoding: `--decoding momst` (baselines: `SVDD`, `SVDD_edit`).

| Flag | Role |
|---|---|
| `--repeatnum` | batch size |
| `--duplicate` | candidate count |
| `--metrics_name` | comma-separated rewards |
| `--metrics_list` | weights, same length |
| `--iteration` | outer test-time cycles |
| `--seq_length` | length when not targeting a PDB |
| `--proteinname` | target PDB for `match_ss` / `crmsd` |

Single objective:

```bash
python refinement.py --decoding momst --repeatnum 20 --duplicate 20 \
  --metrics_name crmsd --metrics_list 1 --proteinname 5KPH --iteration 40
```

Therapeutic-style triple (stability / solubility / confidence):

```bash
python refinement.py --decoding momst --repeatnum 10 --duplicate 20 \
  --metrics_name hydrophobic,surface_expose,plddt --metrics_list 1,1,1 \
  --iteration 20 --seq_length 150
```

Compact scaffold:

```bash
python refinement.py --decoding momst --repeatnum 10 --duplicate 20 \
  --metrics_name globularity,plddt --metrics_list 1,1 --iteration 20 --seq_length 150
```

`match_ss` needs a named example protein in `./datasets`.

## Campaign integration

1. Generate / edit backbone (RFD4 or PROTEAN backends).
2. Inverse-fold a diverse pool (ProteinMPNN or RFD4 inverse folding).
3. Run MoMST (or a Pareto filter that mimics it) on **named, assay-relevant** metrics:
   - binders: pLDDT/ipAE + hydrophobic patch + net charge / SAP
   - enzymes: motif RMSD + pLDDT + surface exposure
4. Keep a **front**, not a winner: 3–5 representatives (high function / balanced / high developability).
5. Feed wet-lab failures into the recency tier so the same aggregation motif is not rediscovered.

## Limitations

- Rewards are proxy metrics (pLDDT ≠ activity; hydrophobicity ≠ expression).
- Memory can lock onto in silico artifacts if the oracle is biased; keep the recency tier and diversity.
- Equal weights (`1,1,1`) are a starting point; set weights from assay priorities, then inspect the Pareto set.
- Computational cost scales with `--iteration` × `--duplicate` × folding oracles.

Pair with `protein-design-filters` so MoMST does not promote sequences that fail hard gates (unpaired Cys, glycosylation motifs, poly-hydrophobic stretches).
