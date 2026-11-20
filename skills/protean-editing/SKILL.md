---
name: protean-editing
description: >-
  Compiles natural-language protein engineering goals into geometric constraints
  and evaluates text-guided 3D structure editing (PROTEAN). Use when the user
  wants local edits, motif preservation, scaffold conservation, LLM agents plus
  generative backends, scope leakage checks, or PROTEAN-style G1/G2/G3 gates.
---

# PROTEAN text-guided editing

Paper: Su et al., 2026. *PROTEAN: Beyond De Novo Generation — Can LLM Agents Perform Text-Guided Protein Structure Editing?* https://arxiv.org/abs/2610.17765

## Problem

Most real campaigns start from a characterized protein. The job is **local 3D editing** (affinity, immunogenicity, stability) while preserving catalytic motifs and scaffold. Generative backends need low-level geometry, not English. PROTEAN benchmarks the two-stage agent:

1. LLM compiles instruction `u` → backend-specific constraints.
2. Frozen generator remodeels the backbone.

380 sub-tasks: 350 length-preserving edits (7 types × 50) + 30 motif-scaffolding (MotifBench). Abstract also groups them as six intervention families.

## Task tuple

```
τ = ((s, X), M, u)
M = (M_edit, M_relax, M_prot)
```

- `s`: sequence; `X`: backbone N, CA, C, O.
- Length-preserving tasks: no indels, no noncanonicals; substitutions only in `M_edit`; `M_prot` identities frozen.
- After scaffold-anchored superposition, protected backbone RMSD ≤ **0.1 Å** (paper task contract). Designs that “succeed” by rebuilding the protected region **fail**.

Colors in PROTEAN figures: red = edit core, blue = relaxation shell, gray = fixed scaffold, purple = reference fragment.

## Task families and G1 bars

| Family | G1 (target attainment) |
|---|---|
| Helix enrichment | DSSP Δf = f_helix − f_sheet ≥ 0.35 |
| Sheet enrichment | Δf ≥ 0.10 |
| Inter-domain approach / separation | median distance change ≥ 1.0 Å in scaffold frame |
| Fragment mimicry | mapped backbone RMSD ≤ 2.0 Å, coverage ≥ 80%, ≥ 3 aligned residues |
| Core burial | median RSA drop ≥ 0.10 and RSA ≤ 0.40 |
| Head-to-tail cyclization | exactly one new C–N, distance 1.20–1.50 Å |
| Motif scaffolding (variable length) | motif backbone RMSD ≤ 1.0 Å + self-consistency |

## Two-tier evaluation (always run both)

**Geometric (do this on generated coordinates, before sequence):**

- **G1** target attainment (table above)
- **G3** scope compliance: protected sequence, mapping, and geometry. Off-target drift invalidates a G1 pass.

In-scaffold editing success in the paper is **G1 ∧ G3**. Best reported joint geometric compliance ~46% (Claude-Opus-5, pass@16, best backend); many models <28%.

**Biophysical:**

- **G2** ProteinMPNN on `M_edit` (protect context) → ESMFold refold (4 sequences). Need at least one below frozen whole-chain and protected RMSD cutoffs. Motif scaffolding additionally requires whole-chain Cα RMSD < 2.0 Å and motif Cα RMSD < 1.0 Å.

Strict end-to-end biophysical success starts around **7%**. Ablations: LLM compilation is already accurate. Bottlenecks are **backend scope leakage** and **refolding failure**, plus non-local contact strain and Pareto deadlocks on hard remodeling.

Cyclization: G1 can look perfect while joint G1∧G3 is 0% in the reported table — topology edits are especially leaky. Treat cyclization as high-risk.

## Agent procedure

```
- [ ] Write u as a scientific instruction (what changes, what is invariant)
- [ ] Partition residues: edit / relax / protect (include a shell around the edit)
- [ ] Choose backend (see routing)
- [ ] Compile native spec (YAML/JSON/flags). Show it to the user before GPU
- [ ] Sample K backbones (K≥16 when using PROTEAN-like pass@K)
- [ ] Superpose on M_prot; compute G1 and G3
- [ ] Reject G3 failures even if G1 looks good
- [ ] Inverse fold + refold only G1∧G3 passers (G2)
- [ ] If G3 fails: smaller edit, larger relax shell, stronger C_CRD on protect, partial diffusion
- [ ] If G2 fails: more MPNN sequences, native-biased decoding, or discard backbone
```

## Backend routing (from PROTEAN related work)

Heterogeneous, not interchangeable:

- **BoltzGen**: native head-to-tail cyclization / covalent constraints
- **ODesign**: partial diffusion, hotspot-conditioned inter-domain distances
- **RFdiffusion3**: general backbone under spatial bounds
- **RFD4-Proteina** (preferred when installed): YAML `contig` + `C_CRD`/`C_CRDGRP`/`C_SEQ`/`C_HOT`, partial diffusion, search — use `rfd4-proteina` skill

The LLM's job is to **select and compile**, not to invent Cartesian coordinates by hand.

## Implications for autonomous agents

- Decouple geometry from biophysics. Never report a single “success” score.
- Log malformed backend outputs separately from algorithmic failure.
- Do not optimize G1 by violating G3.
- Sequence-only language models are the wrong interface for 3D conservation contracts.
- Pair this skill with `protein-design-filters` after G2, and `momst-multiobjective` when affinity vs solubility vs stability conflict.

## Compile template

```
Instruction u:
Edit scope M_edit:
Relaxation shell M_relax:
Protected M_prot (identities + RMSD contract):
Backend:
Native constraints:
Negative constraints (must not happen):
G1 metric:
G3 metric:
G2 protocol:
```
