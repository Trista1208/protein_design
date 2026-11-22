# Submitted top 20 — and why each one is in

Ranked by **reproducibility (non-zero replicates out of 5) → median → best**, every number re-measured from the folded structures and scored **per replicate**.

Reproducibility leads the sort deliberately: a design that scores 0.88 on one seed and 0.00 on four is not better than one that scores 0.83 five times out of five. The organisers run five replicates too.

All 20 pass the competition gates: ≤260 aa, Asp102 and Gln46 intact, TM-score ≤0.70 and sequence identity ≤0.70 to published enzymes. Stability is **not** a competition criterion and was never used to rank — it is shown only to flag anything unbuildable.

| # | design | mutations | repro | median | best | worst | pLDDT | iPDE | TM | stability | why it is in |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `L3_I23V_M45Y_V84I` | I23V+M45Y+V84I | 5/5 | **0.8421** | 0.8823 | 0.7998 | 0.979 | 0.47 | 0.170 | ✅ | **highest median and highest worst replicate** — scores on **every** seed with no term ever zeroed |
| 2 | `L_E_shell_I23L` | I23L | 5/5 | **0.8330** | 0.8617 | 0.7671 | 0.979 | 0.44 | 0.170 | ✅ | scores on **every** seed with no term ever zeroed |
| 3 | `L2_triple_I23L_M88F_V84L` | I23L+V84L+M88F | 5/5 | **0.8316** | 0.8748 | 0.7807 | 0.976 | 0.52 | 0.174 | ✅ | scores on **every** seed with no term ever zeroed |
| 4 | `L2_pair_V84I_I23V` | I23V+V84I | 5/5 | **0.8312** | 0.8536 | 0.7622 | 0.977 | 0.50 | 0.170 | ✅ | scores on **every** seed with no term ever zeroed |
| 5 | `L2_pair_M88F_M45Y` | M45Y+M88F | 5/5 | **0.8251** | 0.8487 | 0.7774 | 0.978 | 0.49 | 0.144 | ⚠️ | scores on **every** seed with no term ever zeroed; stability flagged (ddG +1.63 kcal/mol) |
| 6 | `L3_M45Y_V84I_M88F` | M45Y+V84I+M88F | 5/5 | **0.8249** | 0.8630 | 0.7818 | 0.976 | 0.53 | 0.174 | ⚠️ | scores on **every** seed with no term ever zeroed; stability flagged (ddG +1.69 kcal/mol) |
| 7 | `L2_pair_I23V_M45Y` | I23V+M45Y | 5/5 | **0.8106** | 0.8866 | 0.6975 | 0.978 | 0.49 | 0.169 | ✅ | scores on **every** seed with no term ever zeroed |
| 8 | `L2_pair_I23L_V84L` | I23L+V84L | 5/5 | **0.8044** | 0.8839 | 0.7448 | 0.978 | 0.44 | 0.137 | ✅ | scores on **every** seed with no term ever zeroed |
| 9 | `L2_pair_I23L_V84I` | I23L+V84I | 5/5 | **0.7699** | 0.8760 | 0.7249 | 0.978 | 0.45 | 0.160 | ✅ | scores on **every** seed with no term ever zeroed |
| 10 | `L2_pair_V84I_M45Y` | M45Y+V84I | 4/5 | **0.8263** | 0.8592 | 0.8213 | 0.980 | 0.47 | 0.173 | ✅ | scores on 4 of 5; `orient_oxy` zeroes the rest |
| 11 | `L2_triple_I23L_V84I_M88F` | I23L+V84I+M88F | 4/5 | **0.8235** | 0.8508 | 0.7977 | 0.973 | 0.55 | 0.136 | ✅ | scores on 4 of 5; `orient_oxy` zeroes the rest |
| 12 | `L3_I23V_M45Y_V84L` | I23V+M45Y+V84L | 4/5 | **0.8205** | 0.8967 | 0.8149 | 0.976 | 0.51 | 0.170 | ✅ | scores on 4 of 5; `oxy` zeroes the rest |
| 13 | `L_E_shell_V84I` | V84I | 4/5 | **0.7907** | 0.8359 | 0.7790 | 0.977 | 0.49 | 0.144 | ✅ | scores on 4 of 5; `orient_oxy,oxy` zeroes the rest |
| 14 | `L2_pair_I23L_M88A` | I23L+M88A | 4/5 | **0.6387** | 0.8090 | 0.6148 | 0.975 | 0.65 | 0.170 | ⚠️ | scores on 4 of 5; `oxy` zeroes the rest; stability flagged (ddG +2.21 kcal/mol) |
| 15 | `L_E_shell_I23V` | I23V | 3/5 | **0.8263** | 0.8730 | 0.8263 | 0.977 | 0.51 | 0.170 | ✅ | scores on 3 of 5; `orient_oxy` zeroes the rest |
| 16 | `L2_triple_V84I_I23V_M88F` | I23V+V84I+M88F | 3/5 | **0.8106** | 0.8657 | 0.8106 | 0.970 | 0.61 | 0.160 | ✅ | scores on 3 of 5; `orient_oxy` zeroes the rest |
| 17 | `L2_pair_V84I_M88F` | V84I+M88F | 3/5 | **0.7320** | 0.7676 | 0.7320 | 0.973 | 0.58 | 0.144 | ✅ | scores on 3 of 5; `orient_oxy,oxy` zeroes the rest |
| 18 | `L2_pair_I23V_R127G` | I23V+R127G | 3/5 | **0.5117** | 0.7723 | 0.5117 | 0.976 | 0.49 | 0.173 | ⚠️ | scores on 3 of 5; `orient_oxy` zeroes the rest; stability flagged (ddG +1.73 kcal/mol) |
| 19 | `L_E_shell_V84M` | V84M | 2/5 | **0.0000** | 0.8465 | 0.7363 | 0.979 | 0.52 | 0.144 | ✅ | partial pass (2/5) — kept because an occupied slot beats an empty one |
| 20 | `L_E_shell_V84L` | V84L | 2/5 | **0.0000** | 0.8258 | 0.7929 | 0.977 | 0.52 | 0.144 | ✅ | partial pass (2/5) — kept because an occupied slot beats an empty one |

### Merged duplicates

Identical sequences are one design — SUBMISSION.md says duplicates are scored once and waste a slot.

- `L2_triple_V84I_I23V_M88F` also generated as `L3_I23V_V84I_M88F`
- `L_D_seq_anchor_R127A_M45A` also generated as `L2_pair_R127A_M45A`
- `ins080_F42Y_H103Y` also generated as `ins080_r103_42_F42Y_H103Y`

### Beyond the 20

31 designs qualify in total. Every one of them has a FASTA and all five folded PDBs in [`designs/`](designs/), and keeps its measured numbers in [`ARCHIVE.md`](ARCHIVE.md) whether or not it made the cut.
