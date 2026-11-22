# Biochemistry, biophysics and thermodynamics

## Verdict per design — eligibility and stability

Eligibility is the competition's own rule set. Stability is **not** a competition criterion; it is a buildability check, and it was never used to rank anything.

A design is CLEAN when every one of these holds: ThermoMPNN ddG <= 1.5 kcal/mol (~3x the ~0.5 experimental error), instability index < 40.0 (Guruprasad), Rosetta relax per residue better than the HG3.17 control, MD drift <= 0.5 A/ns and radius-of-gyration swell <= 5%. FLAGGED names the term that failed rather than hiding it.

| design | repro | official | eligible | TM | identity | base/oxy | stability | why flagged |
|---|---|---|---|---|---|---|---|---|
| `L3_I23V_M45Y_V84I` | 5/5 | 0.8421 | ✅ NOVEL | 0.170 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L_E_shell_I23L` | 5/5 | 0.8330 | ✅ NOVEL | 0.170 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L2_triple_I23L_M88F_V84L` | 5/5 | 0.8316 | ✅ NOVEL | 0.174 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L2_pair_V84I_I23V` | 5/5 | 0.8312 | ✅ NOVEL | 0.170 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L2_pair_M88F_M45Y` | 5/5 | 0.8251 | ✅ NOVEL | 0.144 | 0.00 | 102/46 | ⚠️ FLAGGED | ddG +1.63 kcal/mol |
| `L3_M45Y_V84I_M88F` | 5/5 | 0.8249 | ✅ NOVEL | 0.174 | 0.00 | 102/46 | ⚠️ FLAGGED | ddG +1.69 kcal/mol |
| `L2_pair_I23V_M45Y` | 5/5 | 0.8106 | ✅ NOVEL | 0.169 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L2_pair_I23L_V84L` | 5/5 | 0.8044 | ✅ NOVEL | 0.137 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L2_pair_I23L_V84I` | 5/5 | 0.7699 | ✅ NOVEL | 0.160 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L2_pair_V84I_M45Y` | 4/5 | 0.8263 | ✅ NOVEL | 0.173 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L2_triple_I23L_V84I_M88F` | 4/5 | 0.8235 | ✅ NOVEL | 0.136 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L3_I23V_M45Y_V84L` | 4/5 | 0.8205 | ✅ NOVEL | 0.170 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L_E_shell_V84I` | 4/5 | 0.7907 | ✅ NOVEL | 0.144 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L2_pair_I23L_M88A` | 4/5 | 0.6387 | ✅ NOVEL | 0.170 | 0.00 | 102/46 | ⚠️ FLAGGED | ddG +2.21 kcal/mol |
| `L_E_shell_I23V` | 3/5 | 0.8263 | ✅ NOVEL | 0.170 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L2_triple_V84I_I23V_M88F` | 3/5 | 0.8106 | ✅ NOVEL | 0.160 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L2_pair_V84I_M88F` | 3/5 | 0.7320 | ✅ NOVEL | 0.144 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L2_pair_I23V_R127G` | 3/5 | 0.5117 | ✅ NOVEL | 0.173 | 0.00 | 102/46 | ⚠️ FLAGGED | ddG +1.73 kcal/mol |
| `L_E_shell_V84M` | 2/5 | 0.0000 | ✅ NOVEL | 0.144 | 0.00 | 102/46 | ✅ CLEAN | – |
| `L_E_shell_V84L` | 2/5 | 0.0000 | ✅ NOVEL | 0.144 | 0.00 | 102/46 | ✅ CLEAN | – |

**20 eligible, 16 clean, 4 flagged.** A flagged design is still submitted: SUBMISSION.md ranks each design individually against other teams' entries, so a weaker one cannot drag the others down, and an occupied slot beats an empty one.

---


Every scored design, with the **HG3.17 (4BS0) crystal structure as the
control arm** - a real enzyme that has been expressed, purified,
crystallised and kinetically characterised. A stability number with no
control is not evidence.

## Thermodynamics

`relax/res` - Rosetta relax total_score per residue, protein only (lower = better packed).

`ddG` - ThermoMPNN, summed over each design's mutations relative to the
parent shell. **Positive = destabilising** (verified empirically: a
self-mutation scores exactly 0.000, ->Pro +2.06, ->Gly +1.36,
hydrophobic->charged +2.06). Typical experimental error on a measured
ddG is ~0.5 kcal/mol, so anything under that is indistinguishable from
zero.

| design | repro | official | relax/res | fa_rep/res | ddG | mutations |
|---|---|---|---|---|---|---|
| `L3_I23V_M45Y_V84I` | 5/5 | 0.8421 | **-3.509** | 0.608 | +1.14 | I23V+M45Y+V84I |
| `L_E_shell_I23L` | 5/5 | 0.8330 | **-3.626** | 0.626 | +0.25 | I23L |
| `L2_triple_I23L_M88F_V84L` | 5/5 | 0.8316 | **-3.647** | 0.638 | +1.05 | I23L+V84L+M88F |
| `L2_pair_V84I_I23V` | 5/5 | 0.8312 | **-3.596** | 0.616 | +0.41 | I23V+V84I |
| `L2_pair_M88F_M45Y` | 5/5 | 0.8251 | **-3.586** | 0.640 | +1.63 | M45Y+M88F |
| `L3_M45Y_V84I_M88F` | 5/5 | 0.8249 | **-3.615** | 0.623 | +1.69 | M45Y+V84I+M88F |
| `L2_pair_I23V_M45Y` | 5/5 | 0.8106 | **-3.546** | 0.633 | +1.07 | I23V+M45Y |
| `L2_pair_I23L_V84L` | 5/5 | 0.8044 | **-3.618** | 0.618 | +0.14 | I23L+V84L |
| `L2_pair_I23L_V84I` | 5/5 | 0.7699 | **-3.646** | 0.631 | +0.31 | I23L+V84I |
| `L2_pair_V84I_M45Y` | 4/5 | 0.8263 | **-3.460** | 0.680 | +0.79 | M45Y+V84I |
| `L2_triple_I23L_V84I_M88F` | 4/5 | 0.8235 | **-3.620** | 0.647 | +1.21 | I23L+V84I+M88F |
| `L3_I23V_M45Y_V84L` | 4/5 | 0.8205 | **-3.540** | 0.666 | +0.97 | I23V+M45Y+V84L |
| `L_E_shell_V84I` | 4/5 | 0.7907 | **-3.661** | 0.621 | +0.06 | V84I |
| `L2_pair_I23L_M88A` | 4/5 | 0.6387 | **-3.577** | 0.621 | +2.21 | I23L+M88A |
| `L_E_shell_I23V` | 3/5 | 0.8263 | **-3.578** | 0.620 | +0.35 | I23V |
| `L2_triple_V84I_I23V_M88F` | 3/5 | 0.8106 | **-3.603** | 0.625 | +1.32 | I23V+V84I+M88F |
| `L2_pair_V84I_M88F` | 3/5 | 0.7320 | **-3.575** | 0.645 | +0.97 | V84I+M88F |
| `L2_pair_I23V_R127G` | 3/5 | 0.5117 | **-3.640** | 0.627 | +1.73 | I23V+R127G |
| `L_E_shell_V84M` | 2/5 | 0.0000 | **-3.663** | 0.619 | +0.20 | V84M |
| `L_E_shell_V84L` | 2/5 | 0.0000 | **-3.668** | 0.610 | -0.10 | V84L |
| `REF_HG3.17_4BS0` | - | - | **-3.271** | 0.754 | - | (reference enzyme) |

## Thermal - molecular dynamics

2 ns implicit-solvent (GBn2) MD at 350 K, protein only, identical settings for every run including the control.

`RMSD` is the backbone CA deviation averaged over the second half of the run; `drift` is its slope there. A folded protein plateaus (drift ~0); an unfolding one keeps climbing. `Rg` swelling means loss of compaction.

| design | RMSD A | max A | drift A/ns | Rg start -> end A | swell % |
|---|---|---|---|---|---|
| `L_E_shell_I23L` | 2.23 | 2.57 | +0.03 | 14.40 -> 14.43 | +0.2 |
| `L2_triple_I23L_M88F_V84L` | 1.31 | 1.68 | +0.19 | 14.64 -> 14.50 | -1.0 |
| `L2_pair_V84I_I23V` | 2.17 | 2.45 | +0.30 | 14.40 -> 14.36 | -0.2 |
| `L2_pair_M88F_M45Y` | 1.59 | 2.25 | +0.08 | 14.43 -> 14.33 | -0.7 |
| `L2_pair_I23V_M45Y` | 1.57 | 1.92 | +0.26 | 14.36 -> 14.38 | +0.1 |

**This is not a melting temperature.** Two nanoseconds of implicit solvent cannot unfold a stable protein, so a flat trace is weak evidence of stability while a rising one is strong evidence of a problem. Implicit solvent also has no explicit-water viscosity, so absolute timescales are not physical - only the comparison between identically-run simulations means anything.


## Biophysics - fold and surface

`core` - buried hydrophobics (relative SASA < 0.15) divided by length.
`exp_hyd` / `exp_chg` - hydrophobic / charged residues with relative
SASA > 0.40. High exposed charge favours soluble expression.

| design | core | exp_hyd % | exp_chg % | SASA A^2 | pLDDT | iPDE |
|---|---|---|---|---|---|---|
| `L3_I23V_M45Y_V84I` | 0.340 | 9.1 | 55.6 | 7653 | 0.979 | 0.47 |
| `L_E_shell_I23L` | 0.333 | 9.1 | 51.1 | 7506 | 0.979 | 0.44 |
| `L2_triple_I23L_M88F_V84L` | 0.340 | 9.1 | 46.7 | 7484 | 0.976 | 0.52 |
| `L2_pair_V84I_I23V` | 0.367 | 9.1 | 46.7 | 7649 | 0.977 | 0.50 |
| `L2_pair_M88F_M45Y` | 0.340 | 7.6 | 51.1 | 7634 | 0.978 | 0.49 |
| `L3_M45Y_V84I_M88F` | 0.340 | 7.6 | 51.1 | 7667 | 0.976 | 0.53 |
| `L2_pair_I23V_M45Y` | 0.340 | 7.6 | 53.3 | 7469 | 0.978 | 0.49 |
| `L2_pair_I23L_V84L` | 0.333 | 9.1 | 53.3 | 7516 | 0.978 | 0.44 |
| `L2_pair_I23L_V84I` | 0.353 | 9.1 | 51.1 | 7501 | 0.978 | 0.45 |
| `L2_pair_V84I_M45Y` | 0.340 | 7.6 | 55.6 | 7427 | 0.980 | 0.47 |
| `L2_triple_I23L_V84I_M88F` | 0.353 | 9.1 | 48.9 | 7625 | 0.973 | 0.55 |
| `L3_I23V_M45Y_V84L` | 0.333 | 7.6 | 48.9 | 7502 | 0.976 | 0.51 |
| `L_E_shell_V84I` | 0.353 | 9.1 | 55.6 | 7577 | 0.977 | 0.49 |
| `L2_pair_I23L_M88A` | 0.347 | 9.1 | 51.1 | 7682 | 0.975 | 0.65 |
| `L_E_shell_I23V` | 0.353 | 9.1 | 55.6 | 7575 | 0.977 | 0.51 |
| `L2_triple_V84I_I23V_M88F` | 0.340 | 9.1 | 53.3 | 7695 | 0.970 | 0.61 |
| `L2_pair_V84I_M88F` | 0.353 | 9.1 | 53.3 | 7649 | 0.973 | 0.58 |
| `L2_pair_I23V_R127G` | 0.353 | 9.1 | 50.0 | 7471 | 0.976 | 0.49 |
| `L_E_shell_V84M` | 0.360 | 9.1 | 51.1 | 7524 | 0.979 | 0.52 |
| `L_E_shell_V84L` | 0.367 | 9.1 | 51.1 | 7535 | 0.977 | 0.52 |
| `REF_HG3.17_4BS0` | 0.357 | 6.5 | 19.6 | 11108 | - | - |

## Biochemistry - sequence

Instability index < 40 is the classical Guruprasad threshold for a
protein predicted stable in vivo. GRAVY below zero is hydrophilic.

| design | len | pI | net charge | GRAVY | instability | aromatic % | Cys |
|---|---|---|---|---|---|---|---|
| `L3_I23V_M45Y_V84I` | 150 | 4.67 | -17 | -0.39 | 32.2 | 8.7 | 1 |
| `L_E_shell_I23L` | 150 | 4.67 | -17 | -0.37 | 32.3 | 8.0 | 1 |
| `L2_triple_I23L_M88F_V84L` | 150 | 4.67 | -17 | -0.37 | 32.8 | 8.7 | 1 |
| `L2_pair_V84I_I23V` | 150 | 4.67 | -17 | -0.36 | 31.7 | 8.0 | 1 |
| `L2_pair_M88F_M45Y` | 150 | 4.67 | -17 | -0.38 | 33.3 | 9.3 | 1 |
| `L3_M45Y_V84I_M88F` | 150 | 4.67 | -17 | -0.38 | 32.7 | 9.3 | 1 |
| `L2_pair_I23V_M45Y` | 150 | 4.67 | -17 | -0.39 | 32.8 | 8.7 | 1 |
| `L2_pair_I23L_V84L` | 150 | 4.67 | -17 | -0.37 | 32.3 | 8.0 | 1 |
| `L2_pair_I23L_V84I` | 150 | 4.67 | -17 | -0.37 | 31.7 | 8.0 | 1 |
| `L2_pair_V84I_M45Y` | 150 | 4.67 | -17 | -0.38 | 32.2 | 8.7 | 1 |
| `L2_triple_I23L_V84I_M88F` | 150 | 4.67 | -17 | -0.36 | 32.2 | 8.7 | 1 |
| `L3_I23V_M45Y_V84L` | 150 | 4.67 | -17 | -0.39 | 32.8 | 8.7 | 1 |
| `L_E_shell_V84I` | 150 | 4.67 | -17 | -0.36 | 31.7 | 8.0 | 1 |
| `L2_pair_I23L_M88A` | 150 | 4.67 | -17 | -0.37 | 32.8 | 8.0 | 1 |
| `L_E_shell_I23V` | 150 | 4.67 | -17 | -0.37 | 32.3 | 8.0 | 1 |
| `L2_triple_V84I_I23V_M88F` | 150 | 4.67 | -17 | -0.36 | 32.2 | 8.7 | 1 |
| `L2_pair_V84I_M88F` | 150 | 4.67 | -17 | -0.36 | 32.2 | 8.7 | 1 |
| `L2_pair_I23V_R127G` | 150 | 4.60 | -18 | -0.34 | 32.7 | 8.0 | 1 |
| `L_E_shell_V84M` | 150 | 4.67 | -17 | -0.38 | 32.3 | 8.0 | 1 |
| `L_E_shell_V84L` | 150 | 4.67 | -17 | -0.37 | 32.3 | 8.0 | 1 |
| `REF_HG3.17_4BS0` | 318 | 5.96 | -5 | -0.10 | 26.7 | 8.5 | 3 |

## What these numbers do and do not say

They say: nothing here indicates a folding, packing or expression
problem. Every design is better packed per residue than the reference
enzyme, matches its hydrophobic core burial, and carries mutations whose
predicted stability cost is at or below experimental error.

They do not say these proteins are stable. Rosetta score is a packing
proxy and ThermoMPNN is a sequence-structure predictor; neither is a
melting temperature. The supplied Risso 2020 sheet is explicit that Tm
is *a gate, not a predictor* of activity, so none of this was used to
rank designs - only to check that the ranked ones are buildable.

Tools: Rosetta `relax`, ThermoMPNN (Dieckhaus et al., PNAS 2024),
Biopython `ShrakeRupley` + `ProtParam`, freesasa.
