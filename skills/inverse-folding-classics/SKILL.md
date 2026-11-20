---
name: inverse-folding-classics
description: >-
  Uses classic computational protein design tools in the modern loop:
  ProteinMPNN, RFdiffusion, ESMFold, AlphaFold, Boltz, ColabDesign, and
  RosettaCommons sequence/structure filters. Use when inverse folding,
  self-consistency, hallucination, or legacy pipelines are needed alongside
  RFD4-Proteina.
---

# Inverse folding and classic tools

Use these when RFD4 is unavailable, for orthogonal validation, or for sequence design on a fixed backbone.

## Canonical loop (Watson/Dauparas era)

```
backbone (RFdiffusion / RFD4 / fragment) 
  → sequence (ProteinMPNN, optionally soluble-MPNN)
  → fold (AF2 / ESMFold / AF3 / Boltz)
  → keep if refold matches design (scRMSD + pLDDT + PAE)
```

RFD4 can replace both backbone and inverse folding (`generation=inverse_folding`) but MPNN remains the PROTEAN G2 standard and a useful independent check.

## ProteinMPNN (Dauparas et al., 2022)

- Input: backbone coordinates; output: sequences with native-like identities and packing.
- **Fix** motif / target / disulfide positions (`fixed_positions` / tied residues). Only redesign `M_edit` or the designed chain.
- Sample temperature ~0.1–0.2 for conservative; higher for diversity. PROTEAN uses 4 sequences per backbone as a finite G2 budget.
- Soluble-MPNN variant when cytosolic expression is the bottleneck.
- Do not MPNN the receptor sequence in binder design unless you intend to redesign the target.

## Backbone generators

| Tool | Role |
|---|---|
| RFdiffusion (Watson et al., 2023) | Pioneer diffusion backbones, motif scaffolding, binder contig syntax |
| RFdiffusion3 | PROTEAN backend; spatial bounds |
| RFD4-Proteina | Default modern generator (see `rfd4-proteina`) |
| BindCraft / hallucination (ColabDesign) | High pass rate, low throughput; RFD4 paper compares against this Pareto |
| Chroma, ESM3 | Alternative generators; still need MPNN+refold |

RFdiffusion contig strings (`50-120/0 A17-131`) map conceptually onto RFD4 YAML `contig` plus conditions.

## Structure oracles

| Tool | Strength | Weakness |
|---|---|---|
| AlphaFold2 + initial-guess / complex | Classic binder designability | Slow; protein-only |
| AlphaFold3 / OpenFold3 | Ligands, NA, complexes; RFD4 rewards | Setup/licensing; compute |
| ESMFold | Fast monomer G2 | Weak on complexes / ligands |
| Boltz / BoltzGen | Co-folding; cyclization constraints (BoltzGen) | Not a substitute for disordered-target ligand screens (see structure-free) |

Self-consistency is necessary, not sufficient. Designed binders can refold and still not bind; enzymes can pass motif RMSD and be inactive.

## ColabDesign / RosettaCommons

- ColabDesign: AF-hallucination, partial diffusion, MPNN, loss terms (pLDDT, PAE, radius of gyration, contact).
- Rosetta: FastRelax, interface energy, packstat — optional late filter, not the generator of first resort.
- Keep versions pinned; mixing AF2 vs AF3 scores across rounds invalidates thresholds.

## How this sits in a campaign

1. **Generate** with RFD4 (or RFdiffusion if that is what is installed).
2. **Sequence** with RFD4 inverse folding **and/or** ProteinMPNN; take union, then unique.
3. **Refold** with the oracle that matches the molecular type.
4. **Filter** with `protein-design-filters`.
5. **Multi-objective** with `momst-multiobjective` on the sequence pool.
6. **Orthogonal check**: if RFD4 and MPNN+ESMFold disagree on a backbone, treat it as high risk.

Cite: ProteinMPNN (Dauparas 2022); RFdiffusion (Watson 2023); AF2 (Jumper 2021); AF3 (Abramson 2024); ESMFold (Lin 2023).
