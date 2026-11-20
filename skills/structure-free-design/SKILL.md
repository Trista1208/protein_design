---
name: structure-free-design
description: >-
  Chooses structure-free approaches (Talus Bio Ptarmigan-1, MARMOT chemoproteomics)
  when 3D folding is the wrong abstraction: intrinsically disordered proteins,
  cryptic pockets, poorly structured regulome targets, and massive ligand screens.
  Use when the user mentions Talus, Ptarmigan, structure-free, IDPs, or
  proteome-scale virtual screening without a pocket.
---

# Structure-free design

Site: https://www.talus.bio/  
Blog: https://blog.talus.bio/p/welcome-to-the-structure-free-future (Fondrie, Federation, Pino)

## When this path is the default

Use structure-free methods when **any** of these hold:

- Target is intrinsically disordered or only folds upon binding (and you need **ligands**, not a de novo miniprotein foldamer).
- No reliable pocket / cryptic site / conformational ensemble poorly captured by AF/Boltz.
- Goal is ranking **many compounds** (billions–trillions) or proteome selectivity, not generating a binder backbone.
- You have or can get chemoproteomic / cellular engagement labels.

Keep **structure-based protein design** (RFD4, PROTEAN) when the product is a protein (binder, enzyme, scaffold). Ptarmigan-1 is a **protein–small-molecule engagement** model, not a backbone generator.

## What Ptarmigan-1 does

- Maps every protein residue and every compound into a shared **256-D** space; proximity ≈ engagement.
- No 3D pose required at train or inference time.
- Screening becomes a **lookup** on a built map, not a per-pair co-fold.
- Reported: ~3.4B compounds × 20,431 proteins; top-ligand retrieval in ~20 H100-hours after library embedding; ~5000× throughput vs Boltz-2 on the cited EGFR H100 benchmark. Predictions still need experiment.
- Site-resolved: which compound **and which residue neighborhood**.
- Competitive with co-folding on standard sites; **better on a novel STAT6 site** in the disclosed Pfizer-patent haystack (Boltz-2 ~random there).
- Emergent: ranks KRAS G12C ligand generations despite withholding KRAS G12C; WT reversion ablates predicted engagement.

Training uses **binary engagement labels**, including covalent chemoproteomic residue maps — not affinities. Do not treat scores as Kd.

## MARMOT (why it works on messy proteins)

Talus functional proteomics in **native human cells**: compound treatment → protein localization / chromatin-association dynamics + covalent site maps. Dense labels on poorly structured regulome proteins (transcription factors, etc.) are the inductive bias structure-based models lack.

## Campaign pattern

```
- [ ] Decide product: protein binder vs small molecule
- [ ] If protein binder to an IDP epitope that can fold on binding: RFD4 disordered-peptide binder campaign (they showed pM YAP1/Tau binders)
- [ ] If ligand / degrader / covalent against IDP or cryptic site: Ptarmigan-1 rank → experimental chemoproteomics / biophysics
- [ ] Check proteome selectivity in the same embedding (off-target residues)
- [ ] Do not "validate" a disordered target by forcing a confident AF pocket and docking
```

## Complementary, not competing

| Need | Tool |
|---|---|
| De novo miniprotein to a disordered peptide | RFD4-Proteina |
| Local edit of a folded binder | PROTEAN + RFD4 partial diffusion |
| Sequence Pareto (solubility vs pLDDT) | MoMST |
| Ligand for unstructured / cryptic proteome | Ptarmigan-1 |
| Cellular mechanism / localization | MARMOT-style assay |

## Limitations / ethics

- Early access; not a local open generator like RFD4.
- Binary labels ≠ affinity or functional inhibition.
- Covalent-heavy training may bias toward residue types seen in chemoproteomics.
- Dual-use: do not help design toxins or uncontained pathogens. Stay in disclosed therapeutic/research targets.

If the user actually wants a **protein** against an IDP, switch to `rfd4-proteina` (disordered peptide binder recipes) rather than embedding-space ligand ranking.
