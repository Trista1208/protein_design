# Top 20 finalists

Ranked by **reproducibility (n/5 non-zero replicates) -> median -> best**.
Every number re-measured from the folded structures and scored **per replicate**
(never by averaging geometry).

| # | design | repro | median | best | worst | pLDDT | iPDE | len | gates | limiting term |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `L_E_shell_I23L` | 5/5 | 0.8330 | 0.8617 | 0.7671 | 0.979 | 0.44 | 150 | PASS | - |
| 2 | `L2_triple_I23L_M88F_V84L` | 5/5 | 0.8316 | 0.8748 | 0.7807 | 0.976 | 0.52 | 150 | PASS | - |
| 3 | `L2_pair_V84I_I23V` | 5/5 | 0.8312 | 0.8536 | 0.7622 | 0.977 | 0.50 | 150 | PASS | - |
| 4 | `L2_pair_M88F_M45Y` | 5/5 | 0.8251 | 0.8487 | 0.7774 | 0.978 | 0.49 | 150 | PASS | - |
| 5 | `L2_pair_I23V_M45Y` | 5/5 | 0.8106 | 0.8866 | 0.6975 | 0.978 | 0.49 | 150 | PASS | - |
| 6 | `L2_pair_I23L_V84L` | 5/5 | 0.8044 | 0.8839 | 0.7448 | 0.978 | 0.44 | 150 | PASS | - |
| 7 | `L2_pair_I23L_V84I` | 5/5 | 0.7699 | 0.8760 | 0.7249 | 0.978 | 0.45 | 150 | PASS | - |
| 8 | `L2_pair_V84I_M45Y` | 4/5 | 0.8263 | 0.8592 | 0.8213 | 0.980 | 0.47 | 150 | PASS | orient_oxy |
| 9 | `L2_triple_I23L_V84I_M88F` | 4/5 | 0.8235 | 0.8508 | 0.7977 | 0.973 | 0.55 | 150 | PASS | orient_oxy |
| 10 | `L_E_shell_V84I` | 4/5 | 0.7907 | 0.8359 | 0.7790 | 0.977 | 0.49 | 150 | PASS | orient_oxy,oxy |
| 11 | `L2_pair_I23L_M88A` | 4/5 | 0.6387 | 0.8090 | 0.6148 | 0.975 | 0.65 | 150 | PASS | oxy |
| 12 | `L_E_shell_I23V` | 3/5 | 0.8263 | 0.8730 | 0.8263 | 0.977 | 0.51 | 150 | PASS | orient_oxy |
| 13 | `L2_triple_V84I_I23V_M88F` | 3/5 | 0.8106 | 0.8657 | 0.8106 | 0.970 | 0.61 | 150 | PASS | orient_oxy |
| 14 | `L2_pair_V84I_M88F` | 3/5 | 0.7320 | 0.7676 | 0.7320 | 0.973 | 0.58 | 150 | PASS | orient_oxy,oxy |
| 15 | `L2_pair_I23V_R127G` | 3/5 | 0.5117 | 0.7723 | 0.5117 | 0.976 | 0.49 | 150 | PASS | orient_oxy |
| 16 | `L_E_shell_V84M` | 2/5 | 0.0000 | 0.8465 | 0.7363 | 0.979 | 0.52 | 150 | PASS | orient_oxy,oxy |
| 17 | `L_E_shell_V84L` | 2/5 | 0.0000 | 0.8258 | 0.7929 | 0.977 | 0.52 | 150 | PASS | orient_oxy,oxy |
| 18 | `L2_pair_M88F_V84L` | 2/5 | 0.0000 | 0.8065 | 0.7842 | 0.974 | 0.58 | 150 | PASS | orient_oxy |
| 19 | `L2_pair_I23V_V84M` | 2/5 | 0.0000 | 0.7945 | 0.7272 | 0.976 | 0.57 | 150 | PASS | orient_oxy,oxy |
| 20 | `L2_pair_M88F_M45W` | 2/5 | 0.0000 | 0.7663 | 0.7040 | 0.978 | 0.50 | 150 | PASS | orient_oxy |

A median of 0.0000 with a non-zero best means the design scores on some seeds
and is zeroed on others - a real but non-reproducible result, reported with its
reproducibility stated.

## Fold stability (Rosetta relax, protein only)

Answering: *do these fold and express, given they have no disulfide?*
The control arm is the HG3.17 crystal structure - a real, expressed,
crystallised, kinetically characterised enzyme that **does** carry a disulfide
(Cys255-Cys261, 2.03 A).

| structure | residues | total REU | REU/residue | fa_rep/residue |
|---|---|---|---|---|
| `L_E_shell_V84I` | 150 | -549.1 | **-3.661** | 0.621 |
| `L2_pair_I23L_V84I` | 150 | -546.9 | **-3.646** | 0.631 |
| `L_E_shell_I23L` | 150 | -543.9 | **-3.626** | 0.626 |
| `L2_pair_I23L_V84L` | 150 | -542.7 | **-3.618** | 0.618 |
| `L2_pair_I23V_M45Y` | 150 | -531.9 | **-3.546** | 0.633 |
| `REF_HG3.17_4BS0` | 299 | -978.2 | **-3.271** | 0.754 |

All five designs are **better packed per residue than the reference enzyme**.
Each carries exactly one free cysteine (Cys30), so no disulfide is possible and
none can scramble.
