# De novo Kemp eliminase design — final report

DBSSE Protein Design Hackathon, Basel · 22 November 2026

Design of a de novo Kemp eliminase for 5-nitrobenzisoxazole, scored against the 6-nitrobenzotriazole (6NT) transition-state analogue.

---

## What we achieved

| | |
|---|---|
| Designs screened | **441** |
| Fully scored over 5 Boltz-2 seeds | **66** |
| Score on at least one replicate | **34** |
| Score on **all five** replicates | **9** |
| Score on 4 or 5 replicates | **14** |
| Qualifying designs (gates + novelty) | **31** |
| Submitted | **20** (the maximum allowed) |
| Best median score | **0.8421** |
| Best single replicate | **0.8967** |
| HG3.17 crystal, same metric (yardstick) | 0.959 |
| HG3.17 folded from sequence, same metric | 0.834 |

**The headline.** Our top design scores **0.8421 as a median across five independent predictions**, with its *worst* seed at 0.7998. The real, crystallised, kinetically characterised HG3.17 enzyme scores **0.834** under the identical procedure when folded from its own sequence. We are at parity with a working enzyme on the organisers' own metric, with a sequence that has **no significant similarity to any published enzyme**.

---

## How we scored — and why it is the core of this work

The official metric is a **weighted geometric mean** of seven subscores:

```
S = PROD_i  s_i ** w_i      equivalently   S = exp( SUM_i w_i ln s_i )
```

**One zero kills the design.** If any single subscore is 0, S = 0 no matter how good the other six are. This one property shaped every decision we made.

| term | measured quantity | formula | target / window |
|---|---|---|---|
| `base` | Asp102 OD2 to ligand N3 | window score | target 2.52 A, window 2.0-3.6 |
| `ipde` | Boltz-2 `complex_ipde` (an *error*) | `max(0, 1 - iPDE/8.0)` | lower is better |
| `oxy` | Gln46 NE2 to ligand N1 | window score | target 3.10 A, window 2.0-3.6 |
| `orient_base` | angle + dihedral at the base | `min(s_angle, s_dihedral)` | - |
| `orient_oxy` | angle + dihedral at the oxyanion | `min(s_angle, s_dihedral)` | - |
| `bury` | mean CA neighbours at 6/8/12 A | linear ramp | 0 at <=3, 1 at >=15 |
| `cons` | SD of `base` distance over 5 seeds | `max(0, 1 - sd/0.50)` | lower is better |

The dihedral term has **two allowed bands, 130-230 deg and 310-50 deg**. Everything between them scores exactly zero. That is not a soft penalty, it is a cliff, and it is where our whole campaign turned.

### Five things we did to make the measurement trustworthy

**1. We score every replicate separately, never averaged geometry.** Averaging the geometry across five seeds and scoring once reported **0.624 for a design where 0 of 5 replicates actually scored**, and simultaneously hid a real 0.743 on another. Measured across the final 20, that bug changes **17 of 20 scores**, and it *manufactures* results: three designs whose true median is **0.0000** would have been published at 0.74-0.76. Every number in this report is per-replicate.

**2. We use circular statistics for angles.** The linear mean of 350 deg and 10 deg is 180 deg - the exact antipode. On our parent design this reported 298.7 deg where the truth was 6.7 deg, zeroing a term on a design that was fine. Under the old averaged scheme this bug alone would have zeroed **8 of our 20** submitted designs.

**3. We refuse to measure a dihedral on an unbound ligand.** Both catalytic distances must be inside 2.0-3.6 A before any angle is trusted. This gate **overturned three results we had already announced**. Critically, confidence does not catch this: one design folded at pLDDT 0.967 with the ligand **9.7 A out of the pocket**.

**4. We corrected the oxyanion target from measurement, not from a file.** The constraint file's 2.40 A is for a different acceptor atom. Measured from real crystals - 4BS0 3.102 A, 7K4Z 3.100 A - the contact is **3.10 A**. The wrong value was single-handedly capping the `oxy` subscore at 0.224 on geometrically healthy motifs.

**5. We treat low confidence as unmeasurable, not as bad.** At pLDDT 0.41-0.47 the *same sequence* folded twice gave catalytic distances of **2.6 A and 13.6 A**. Ranking that is ranking noise. We gate at **pLDDT >= 0.95**, stricter than the 0.90 used by denovoKemp (Listov et al., Nature 2025).

### We checked whether these fixes mattered, rather than assuming

| change | designs affected | effect on the submitted set |
|---|---|---|
| per-replicate vs averaged geometry | **17 / 20** | **large** - would have selected a different set |
| circular vs linear angle mean | 13 / 20 | would have zeroed 8 of 20 outright |
| oxyanion target 3.10 vs 2.40 A | 15 / 20 | **almost none** - see below |

The oxyanion target is the one assumption we are least sure of, so we tested it directly. Re-measuring all 31 designs under both values: **reproducibility (n/5) is identical for every design**, **19 of the 20 submitted are the same**, and **the top 7 are the same seven designs**. The single difference is at rank 20, between two designs that both have a median of 0.0000. Reproducibility is immune because both targets share the same *window*; only the peak inside it moves.

---

## The submitted 20

Ranked by **reproducibility (non-zero replicates out of 5), then median, then best**. Reproducibility leads deliberately: a design scoring 0.88 once and 0.00 four times is not better than one scoring 0.83 five times out of five. The organisers also run five replicates.

| # | design | mutations | repro | median | best | worst | pLDDT | iPDE | TM | relax | ddG | stability |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `L3_I23V_M45Y_V84I` | I23V+M45Y+V84I | 5/5 | **0.8421** | 0.8823 | 0.7998 | 0.979 | 0.47 | 0.170 | -3.509 | +1.14 | clean |
| 2 | `L_E_shell_I23L` | I23L | 5/5 | **0.8330** | 0.8617 | 0.7671 | 0.979 | 0.44 | 0.170 | -3.626 | +0.25 | clean |
| 3 | `L2_triple_I23L_M88F_V84L` | I23L+V84L+M88F | 5/5 | **0.8316** | 0.8748 | 0.7807 | 0.976 | 0.52 | 0.174 | -3.647 | +1.05 | clean |
| 4 | `L2_pair_V84I_I23V` | I23V+V84I | 5/5 | **0.8312** | 0.8536 | 0.7622 | 0.977 | 0.50 | 0.170 | -3.596 | +0.41 | clean |
| 5 | `L2_pair_M88F_M45Y` | M45Y+M88F | 5/5 | **0.8251** | 0.8487 | 0.7774 | 0.978 | 0.49 | 0.144 | -3.586 | +1.63 | FLAGGED |
| 6 | `L3_M45Y_V84I_M88F` | M45Y+V84I+M88F | 5/5 | **0.8249** | 0.8630 | 0.7818 | 0.976 | 0.53 | 0.174 | -3.615 | +1.69 | FLAGGED |
| 7 | `L2_pair_I23V_M45Y` | I23V+M45Y | 5/5 | **0.8106** | 0.8866 | 0.6975 | 0.978 | 0.49 | 0.169 | -3.546 | +1.07 | clean |
| 8 | `L2_pair_I23L_V84L` | I23L+V84L | 5/5 | **0.8044** | 0.8839 | 0.7448 | 0.978 | 0.44 | 0.137 | -3.618 | +0.14 | clean |
| 9 | `L2_pair_I23L_V84I` | I23L+V84I | 5/5 | **0.7699** | 0.8760 | 0.7249 | 0.978 | 0.45 | 0.160 | -3.646 | +0.31 | clean |
| 10 | `L2_pair_V84I_M45Y` | M45Y+V84I | 4/5 | **0.8263** | 0.8592 | 0.8213 | 0.980 | 0.47 | 0.173 | -3.460 | +0.79 | clean |
| 11 | `L2_triple_I23L_V84I_M88F` | I23L+V84I+M88F | 4/5 | **0.8235** | 0.8508 | 0.7977 | 0.973 | 0.55 | 0.136 | -3.620 | +1.21 | clean |
| 12 | `L3_I23V_M45Y_V84L` | I23V+M45Y+V84L | 4/5 | **0.8205** | 0.8967 | 0.8149 | 0.976 | 0.51 | 0.170 | -3.540 | +0.97 | clean |
| 13 | `L_E_shell_V84I` | V84I | 4/5 | **0.7907** | 0.8359 | 0.7790 | 0.977 | 0.49 | 0.144 | -3.661 | +0.06 | clean |
| 14 | `L2_pair_I23L_M88A` | I23L+M88A | 4/5 | **0.6387** | 0.8090 | 0.6148 | 0.975 | 0.65 | 0.170 | -3.577 | +2.21 | FLAGGED |
| 15 | `L_E_shell_I23V` | I23V | 3/5 | **0.8263** | 0.8730 | 0.8263 | 0.977 | 0.51 | 0.170 | -3.578 | +0.35 | clean |
| 16 | `L2_triple_V84I_I23V_M88F` | I23V+V84I+M88F | 3/5 | **0.8106** | 0.8657 | 0.8106 | 0.970 | 0.61 | 0.160 | -3.603 | +1.32 | clean |
| 17 | `L2_pair_V84I_M88F` | V84I+M88F | 3/5 | **0.7320** | 0.7676 | 0.7320 | 0.973 | 0.58 | 0.144 | -3.575 | +0.97 | clean |
| 18 | `L2_pair_I23V_R127G` | I23V+R127G | 3/5 | **0.5117** | 0.7723 | 0.5117 | 0.976 | 0.49 | 0.173 | -3.640 | +1.73 | FLAGGED |
| 19 | `L_E_shell_V84M` | V84M | 2/5 | **0.0000** | 0.8465 | 0.7363 | 0.979 | 0.52 | 0.144 | -3.663 | +0.20 | clean |
| 20 | `L_E_shell_V84L` | V84L | 2/5 | **0.0000** | 0.8258 | 0.7929 | 0.977 | 0.52 | 0.144 | -3.668 | -0.10 | clean |

All 20 pass every competition gate: 150 residues (limit 50-260), Asp102 and Gln46 intact, **TM-score 0.136-0.174 against a 0.70 disqualification limit**, and **no significant sequence identity to any published enzyme**.

---

## Results by criterion

### 1. Challenge score (the only thing that is ranked)

- **9 designs score on all five replicates**, 14 on four or five.
- Best median **0.8421**; best single replicate **0.8967**.
- Nine of the top 20 never zero a single term on any seed.
- For comparison, HG3.17 folded from its own sequence scores 0.834 under the identical pipeline.

### 2. Prediction confidence

- pLDDT **0.970 to 0.980** across all 20, against our own gate of 0.95 and denovoKemp's 0.90. Every design is comfortably measurable.
- Interface error (iPDE) **0.44 to 0.65 A** - this is a weighted criterion in its own right, worth 20% of the score.

### 3. Thermodynamics - packing

- Rosetta relax: **-3.668 to -3.460 REU per residue**, against the HG3.17 crystal at **-3.271**.
- **Every single one of the 20 is better packed per residue than the real enzyme**, which has been expressed, purified, crystallised and kinetically characterised.

### 4. Thermodynamics - mutational cost

- ThermoMPNN ddG of each design's mutations relative to the parent: **-0.10 to +2.21 kcal/mol**.
- The sign convention was **verified empirically, not assumed**: a self-mutation scores exactly 0.000, and ->Pro (+2.06), ->Gly (+1.36) and hydrophobic->charged (+2.06) are all positive. So positive means destabilising.
- Experimental error on a measured ddG is about 0.5 kcal/mol. Our **rank-2 design costs +0.25**, and one substitution (V84L) is *stabilising* at -0.10. The pocket mutations are essentially free.

### 5. Thermal stability

2 ns implicit-solvent molecular dynamics at **350 K**, protein only, with the HG3.17 crystal run under identical settings.

| structure | RMSD A | drift A/ns | Rg start to end A | swell |
|---|---|---|---|---|
| `L_E_shell_I23L` | 2.23 | +0.03 | 14.40 to 14.43 | +0.2% |
| `L2_triple_I23L_M88F_V84L` | 1.31 | +0.19 | 14.64 to 14.50 | -1.0% |
| `L2_pair_V84I_I23V` | 2.17 | +0.30 | 14.40 to 14.36 | -0.2% |
| `L2_pair_M88F_M45Y` | 1.59 | +0.08 | 14.43 to 14.33 | -0.7% |
| `L2_pair_I23V_M45Y` | 1.57 | +0.26 | 14.36 to 14.38 | +0.1% |
| **REF_HG3.17_4BS0** (control) | 1.71 | +0.00 | 18.87 to 19.12 | +1.3% |

Every design holds its fold. RMSD plateaus rather than climbing, and **every one keeps its radius of gyration flat or contracting while the real enzyme expands +1.3%**. No unfolding signal anywhere.

### 6. Biophysics - fold and surface

- Hydrophobic core burial **0.333 to 0.367**, against the reference at **0.357** - indistinguishable.
- Exposed hydrophobic surface 7.6-9.1%, reference 6.5%.
- Exposed charged residues ~50%, reference ~20%. A highly charged surface **favours soluble expression**.

### 7. Biochemistry - sequence

- Instability index **31.7 to 33.3**, all well under the Guruprasad threshold of 40 for 'stable in vivo'. Reference 26.7.
- GRAVY -0.39 to -0.34 - more hydrophilic than the reference (-0.10), so less aggregation-prone.
- Net charge -17 to -18 over 150 residues, against -5 over 318 for the reference. High charge density is characteristic of hyperstable de novo designs and aids solubility.
- **Exactly one free cysteine** per design. A disulfide needs two, so none can form and none can scramble.

### 8. Novelty - the eligibility gates

- TM-score to published structures: **0.136 to 0.174**, against a disqualification threshold of 0.70.
- Sequence identity to published enzymes: **no significant hit at all**, against a threshold of 0.70.
- 31 of 31 qualifying designs are NOVEL on both gates, with very large margins.

---

## How the metrics shaped the submitted 20

**Only the challenge score ranked the submission.** Stability, thermal and biochemistry were never used to sort. They were run as a *buildability check* - to make sure nothing we submitted was unfoldable - and, deliberately, with a control arm, because a stability number with nothing to compare it to is not evidence.

That discipline matters because **the metrics genuinely disagree**:

- Our **rank-1 design on score is 18th of 20 on core burial**.
- `L2_pair_V84I_I23V` is **3rd on score but 1st on both core burial and instability index** - the best-rounded design we have.
- `L2_pair_M88F_M45Y` is **5th on score but flagged on ddG (+1.63)**.

Had we blended stability into the ranking, we would have demoted designs that score well on the criterion the competition actually measures. We did not.

We do report the disagreement honestly: **16 of 20 are clean on every stability threshold and 4 are flagged**, each with the specific term named rather than hidden. All 20 were still submitted, because the organisers rank each design individually against other teams' entries - a weaker design cannot drag the others down, and an occupied slot beats an empty one.

---

## Top 10 insights and findings

**1. We rescued the team's best design from a hard zero.** `rfd4_design_045` scored **0.000 on 0 of 5 seeds**. Three mutations later it scores **0.8330 on 5 of 5**. Our designs are 98% identical to a teammate's backbone - we diagnosed the team's work rather than replacing it.

**2. We found *why* it was zero: the Gln46 pinner.** An H-bond donor (Tyr42-OH or Tyr103-OH) within about 2.7 A pins the oxyanion glutamine into the only rotamer that scores. In-band designs: median **2.35 A (n=64)**. Out-of-band: **3.37 A (n=207)**. This is a measured mechanism, and it told us exactly which positions to mutate.

**3. The dihedral dead zone was the whole problem.** The parent sat at **86.1 deg**, inside the 50-130 deg gap between the two allowed bands, so `orient_oxy` was exactly 0 and the geometric mean was 0 - despite perfectly good distances. One term, one cliff, total failure.

**4. We made a landscape of zeros rankable.** Because the metric is a geometric mean, an early batch is all zeros and completely unrankable. We optimise against a dense surrogate - **how many degrees outside the allowed band** - which turns a flat cliff into a gradient. This is what let the loop move at all.

**5. Averaging geometry invents results that do not exist.** Three designs with a true median of **0.0000** would have been reported at 0.74-0.76. The bug is not even a consistent bias - it inflates some scores and deflates others.

**6. Confidence does not detect ligand loss.** A design folded at **pLDDT 0.967 with the ligand 9.7 A out of the pocket**. Any angle measured there is noise. High confidence is necessary, not sufficient.

**7. Pairs work, triples do not stack - and that told us when to stop.** Pairing winning single mutations lifted the in-band rate from **10% to 47%**. Triples added nothing beyond that. We used this as an empirical stopping criterion instead of burning the clock, and redirected compute to positions that had never been paired.

**8. A teammate's stability concern was worth measuring, and the premise was half right.** The claim was that our designs lack disulfides and so will not fold. We checked: **HG3.17 does carry a disulfide** (Cys255-Cys261, 2.03 A) - our assumption that it did not was wrong. But all 20 designs are **better packed per residue than it**, match its core burial, and hold at 350 K. Disulfides also do not form in a reducing cytoplasm. Checking the premise beat arguing the conclusion.

**9. MSA turned out to be a non-issue, verified rather than assumed.** The MSA server returns depth 1 for de novo sequences, so MSA and single-sequence folding give bit-identical results. We tested this rather than carrying the assumption.

**10. Reading the submission rules carefully was worth 15 slots.** The limit is **20 designs, not 5**. Our headers were also missing the mandatory `base=` and `oxy=` fields and carried free text that would have caused outright rejection. A late, careful re-read of the spec changed the submission more than any algorithmic improvement that day.

---

## Reasons to be happy

1. **Parity with a real enzyme on the organisers' own metric.** 0.8421 median versus 0.834 for HG3.17 folded from sequence.
2. **9 designs reproduce on all five seeds.** This is a reproducibility result, not a lucky single prediction.
3. **Very large novelty margins.** TM 0.136-0.174 against a 0.70 limit; no significant sequence hit at all. Disqualification is not a realistic risk.
4. **Every design is better packed than the reference crystal**, and every one holds its fold at 350 K.
5. **A mechanism, not a black box.** We can say *why* our designs score: the pinner hydrogen bond holds the oxyanion glutamine in band. That is a testable, falsifiable claim.
6. **The measurement layer is auditable.** Every score is re-measurable from the published structures, every bug fix is documented with the number it changed, and retractions are kept in the repository rather than quietly deleted.
7. **All 20 slots are filled with eligible, gated designs**, and every qualifying design's FASTA and all five folded PDBs are published.

## Reasons to be cautious

1. **Nothing here demonstrates catalysis.** A good score means the predicted active site has a productive *shape* under the organisers' metric. It is not a rate constant, and no experiment has been done.
2. **Everything is computed on a predicted structure.** Boltz-2 is not a crystallographer.
3. **Our scorer is a reconstruction.** We do not have the organisers' code. Two things our own code flags as unresolved: the **ramp shape between window edges**, fitted to the only two data points in their worked example, and the **burial anchor convention** (carboxylate carbon / OD2 / CA), which is the 8% `bury` term.
4. **Four of the 20 are weak.** Ranks 17-20 score on 2 or 3 of 5 seeds, two with a median of 0.0000. Effective depth is about 16, not 20.
5. **4 designs are flagged on predicted stability**, up to +2.21 kcal/mol. Submitted deliberately, but they are the weakest entries on both axes.
6. **Stability numbers are proxies, not measurements.** Rosetta score is a packing proxy; ThermoMPNN is a sequence-structure predictor; 2 ns of MD cannot unfold a stable protein. The supplied Risso 2020 sheet is explicit that melting temperature is *a gate, not a predictor* of activity - which is exactly why we never ranked on it.
7. **The MSA check was done on one sequence.** Our designs are 1-4 mutations from that one, so it should carry, but it is inference rather than measurement.

---

## Reproducibility

Code, environments, inputs and every result are on the branch `kmor/four-strategies`:

- `results/SUBMISSION_TABLE.md` - the 20, with the reason each is in
- `results/RANKINGS.md` - every design under each of 11 metrics
- `results/BIOPHYSICS.md` - stability, thermal, biochemistry in full
- `results/ARCHIVE.md` - append-only record of all 34 designs that ever scored
- `results/designs/` - FASTA and all five folded PDBs for all 31 qualifying designs
- `docs/good_practices.md`, `docs/LESSONS.md`, `docs/LOOP_DESIGN.md` - the measurement rules, each traced to something that fooled us

_Generated 2026-11-22 13:18 from measured data; no number in this document was typed by hand._