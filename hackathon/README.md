# Campaign: de novo Kemp eliminase (DBSSE × Perlin Bio)

[English](README.md) · [中文](README.zh.md)

Hackathon: **Design a Catalyst.** Basel, November 2026. Presenter: Lucas Merlicek (Perlin Bio). Contact: hello@perlin.bio

This is a **motif-scaffolded enzyme** campaign, not a binder campaign and not a local edit of a published Kemp enzyme.

**Start here (no GPU):** [basics/](basics/) — official length/identity gates + TSA vs substrate + published Kemp FASTA. Run `python3 hackathon/basics/scripts/validate_design.py --self-test`.

## Goal

- **Reaction:** Kemp elimination of 5-nitrobenzisoxazole. A catalytic **base** deprotonates carbon → isoxazole ring opens. A **H-bond donor** can stabilize the developing oxyanion.
- **Starting point:** de novo. Do not start from KE07 / KE70 / HG3 or other published Kemp folds.
- **Preserve (the catalytic contract):** Asp as base + Gln as H-bond donor, in productive geometry versus the ligand. Slide 5 uses this pair on both the DFT transition state and the TSA.
- **Change:** the entire scaffold. New sequence, new fold.
- **Scoring ligand:** **6-nitrobenzotriazole** (transition-state analogue). The true DFT TS is the chemistry; the TSA is what Boltz-2 will see.
- **Assay / success bar:** organizer scoring (below) plus ≥5 sequences sent to wet lab regardless of score. In silico pass ≠ activity. Historic de novo Kemp hit rates are low until directed evolution.

## Tool stack

| Stage | Default | Fallback if no RFD4 GPU |
|---|---|---|
| Backbone | RFD4-Proteina enzyme / motif scaffolding (`recycle_strategy=medium/midpeak`) | RFdiffusion3 / RFdiffusion |
| Conditioning | Annotated CIF or YAML: ligand + Asp/Gln motif (`C_CRD` / `C_IDX` / `C_HB`) | Contig + active-site residue spec |
| Sequence | RFD4 inverse folding **and** LigandMPNN / ProteinMPNN with Asp+Gln **fixed** | LigandMPNN / CARBonAra |
| Oracle | **Boltz-2 cofolding** of protein + TSA (hackathon scorer) | AF3 / Chai-1 / Protenix-v2 as extra checks |
| Multi-objective | Pareto on iPDE, geometry, base burial, consistency — not one scalar | MoMST only on sequence proxies after geometry passes |
| Diversity | Foldseek / TM vs each other **and** vs published enzymes | Same |

Do **not** use structure-free ligand models (Ptarmigan). The product is a protein catalyst.

## Conditioning contract

- **M_edit:** whole designed chain except catalytic identities.
- **M_prot:** Asp (base) and Gln (H-bond donor) identities + their relative geometry to the TSA. Motif backbone RMSD target ≤ 1 Å after scaffolding.
- **Relaxation:** first- and second-shell packing around the pocket; not the catalytic identities.
- **Ligand:** TSA 6-nitrobenzotriazole in the pocket during generate + score. Do not optimize against the product 3-cyano-4-nitrophenol.
- **Negative / disqualify:**
  - length > **260 aa**
  - sequence identity **> 70%** to a published enzyme
  - backbone **TM-score > 0.7** to a published structure
  - unpaired Cys unless a disulfide is intentional
  - catalytic base solvent-exposed (fails burial / pKa)

Build one annotated input (CIF) that contains: TSA coordinates, Asp, Gln, and an expandable scaffold length (aim **90–180 aa** first; 250 only if the pocket needs a barrel).

## Sampling plan (compute is scarce)

Organizers will not give tens of thousands of designs. Prefer **quality × a few scaffolds** over a huge undifferentiated library.

1. **Motif library (hours):** 2–4 Asp–Gln–TSA placements (different base approach angles, donor to the oxyanion oxygen). Reject clashing or non-productive geometries **before** backbone generation.
2. **Scaffold (GPU day 1):** ~50–200 backbones per motif. RFD4 `num_replicates` modest; `metrics_tags` enzyme/motif if the oracle env exists, else `'~metrics'` and score locally with Boltz-2.
3. **Inverse fold:** 4–8 sequences per passing backbone; **always fix Asp+Gln**.
4. **Boltz-2 screen:** 1 seed first. Keep only confident ligand placement.
5. **Consistency:** **5 seeds** on the shortlist (this matches the official scorer).
6. **Optional search:** early checkpoint = pocket H-bonds / ligand contacts; late = iPDE / pLDDT. Do not spend the whole budget on search.
7. **Partial diffusion:** 0.25–0.35 on the best generated complexes to tighten the pocket without copying a published fold.

## Evaluation waterfall

Hard gates, then rank. Do not average iPDE with burial until geometry passes.

1. **Length and novelty:** ≤260 aa; identity ≤70%; TM ≤0.7 vs published Kemp / PDB hits.
2. **Motif geometry (G1):** Asp Nε/Oδ approach to the carbon that is deprotonated; Gln H-bond to the analog of the oxyanion; distances/angles/dihedrals vs the “ideal approach” used by Perlin. Ligand clash-free.
3. **Scope (G3):** catalytic identities still Asp+Gln; TSA still in the designed pocket (not a second spurious site).
4. **iPDE:** ligand tokens confidently placed (lower is better). Fail if the TSA floats or binds the surface.
5. **Catalytic base burial:** count Cα within **12 Å** of the base; base should not be solvent-exposed.
6. **Designability:** protein pLDDT high; pocket PAE / protein–ligand error low. scRMSD of apo refold vs design < ~2 Å if you also fold without ligand.
7. **Consistency:** 5 Boltz-2 seeds reproduce pose + metrics (preorganization proxy).
8. **Developability:** no unpaired Cys, no huge hydrophobic patch, reasonable pI for *E. coli* if that is the expression host.
9. **Diversity:** Foldseek cluster; keep 1–2 per cluster.
10. **Shortlist 5+ for lab:** not five near-duplicates. Prefer a Pareto set: best geometry, best consistency, best burial, one smaller scaffold, one more polar pocket.

Official score terms (from the brief):

| Metric | Question it asks |
|---|---|
| iPDE | Is the ligand confidently placed? |
| Geometry | Is base + H-bond + ligand orientation productive? |
| Catalytic base burial | Is the base buried enough to keep a useful pKa? |
| Consistency | Does cofolding reproduce the structure across seeds? |

## Expected failure modes

- **Pretty pocket, dead enzyme:** TSA binding without productive Asp geometry.
- **Exposed base:** good iPDE, but the carboxylate is solvated → wrong pKa.
- **Copied KE/HG fold:** high TM or identity → disqualified even if scores look great.
- **Seed lottery:** one beautiful Boltz-2 pose, four garbage seeds → fails consistency.
- **Ligand hallucination:** low iPDE on a non-pocket surface.
- **Scope leakage:** Asp/Gln mutate or drift during inverse fold.
- **Aggregation / no expression:** hydrophobic core wrapping the nitroarene too aggressively.
- **Overfitting the TSA:** 6-nitrobenzotriazole ≠ live 5-nitrobenzisoxazole TS. Treat scores as a filter, not proof of k_cat.

Realistic expectation: most in silico passers will be inactive in lysate. The 5 lab slots are the real deliverable; directed evolution (if allowed later) is how Kemp enzymes historically became fast.

## Wet-lab round 1

- **N:** at least **5** designs (guaranteed). If more slots appear, add from other clusters, not more copies.
- **DNA:** *E. coli* codon optimize; His-tag; avoid glycosylation sequons.
- **Assay:** Kemp elimination of 5-nitrobenzisoxazole (UV of the phenolic product is the classic readout). Include: no-enzyme, buffer, catalytic-base knockout (Asp→Asn or Ala), and a published positive control if the organizers provide one.
- **Go / no-go:** any rate above background + knockout kills activity → iterate on that scaffold. Flat background on all five → change motif placement, not just MPNN temperature.

## Iterate after data

- Failures with no binding: loosen pocket / add a second donor / change motif rotamers.
- Binding but no activity: inspect geometry vs DFT TS; rotate the base; bury it more.
- Expression deaths: MoMST-style `hydrophobic,surface_expose,plddt` on the sequence of a **fixed** backbone.
- One weak hit: partial-diffuse that complex; do **not** jump to a published Kemp scaffold (TM rule).

## Team checklist

```
- [ ] TSA (6-nitrobenzotriazole) 3D file + Asp/Gln motif CIF
- [ ] Confirm official scoring scripts / Boltz-2 settings with organizers
- [ ] Generate diverse scaffolds (not one TIM-barrel clone)
- [ ] Fix Asp+Gln in every inverse-fold run
- [ ] Waterfall, then 5-seed consistency
- [ ] Novelty check (identity + TM)
- [ ] Pick ≥5 sequences across clusters
- [ ] Submit sequences + models; keep Asp→Ala knockouts as controls
```

## Skills in this repo

Load `skills/protein-design-campaign`, `skills/rfd4-proteina`, `skills/protein-design-filters`, `skills/inverse-folding-classics`. PROTEAN is only for later local edits of *your* hit, not for remodeling a published Kemp enzyme.
