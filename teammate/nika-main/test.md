# Kemp eliminase design campaign — procedure and findings

Cluster project: `pdh-team:/mnt/work/projects/Dejvid`
Local results: `Hackathon 2026/results/`
Last updated: 2026-11-22

Target geometry, from the HG3-17 theozyme (PDB 4BS0), with 6NT
(6-nitrobenzotriazole) as the transition-state analogue:

| contact | atoms | distance | approach angle |
|---|---|---|---|
| **base** | Asp `OD2` → 6NT `N3` | 2.49 Å | `CG-OD2···N3` = 127.3° |
| **donor** | Gln `NE2` → 6NT `N1` | 3.34 Å | `CD-NE2···N1` = 158.1° |

Both are held rigid by RFdiffusion3's `select_fixed_atoms`, so all 72
protein-to-ligand distances in the motif are constraints, not just these two.

---

## Round 1 — de novo backbones

```
 RFdiffusion3          job 35        80 backbones, contig 120-180
   │   3 input defects fixed first (see Corrections)
   ▼
 structural filter                   46 of 80 pass; 7 selected
   │   no chainbreak, no ligand clash, catalytic residues >=5 from a terminus
   ▼
 LigandMPNN            jobs 48/49    7 x 256 = 1792 sequences
   │   T=0.5, omit C/M, side-chain packing ON
   ▼
 clustering (mmseqs)                 1240 distinct at 70% identity
   │
   ▼
 Boltz-2               jobs 52/56    663 predictions (cancelled early)
   │   single-sequence, 1 diffusion sample, --use_potentials
   ▼
 RESULT                              0 / 663 satisfy both contacts
```

### Counts

| stage | n |
|---|---|
| RFD3 backbones generated | 80 |
| passing the structural filter | 46 |
| selected for design | 7 |
| LigandMPNN sequences | 1792 |
| distinct at 70% identity | 1240 |
| Boltz predictions completed | 663 |

---

## Round 2 — externally supplied structures (running)

Eight structures contributed by the team, four of which already satisfy both
contacts — something round 1 never achieved.

```
 8 input structures                  136-252 residues, chains A and G
   │   motif detected geometrically; numbering differs per structure
   ▼
 LigandMPNN            T=0.05, 8 seqs each, catalytic pair fixed
   │   ligand-aware, checkpoint v_32_005_25 (lowest noise)
   ▼
 SolubleMPNN           T=0.5, ENTIRE first shell fixed (16-20 residues)
   │   surface only — it is blind to the ligand
   ▼
 mmseqs vs 180 GB DB   ONE search, TWO purposes
   │   ├─ novelty filter: drop >=70% identity to published (>=50% coverage)
   │   └─ MSA: one a3m per surviving design
   ▼
 Boltz-2               MSA-backed
```

### Input geometry

| structure | residues | Asp | base | Gln | donor |
|---|---|---|---|---|---|
| `boltz_21859852_137_seed_3000` | 252 | A43 | **2.36** | A158 | **3.42** |
| `rfd4_design_045_of3_refolded` | 150 | G102 | **2.86** | G46 | **3.09** |
| `rfd4_design_014_of3_refolded` | 200 | G178 | **3.52** | G153 | **2.93** |
| `boltz_21856114_12_seed_3002` | 237 | A185 | **3.41** | A120 | **3.14** |
| `boltz_21856114_12_seed_3001` | 237 | A185 | 3.52 | A120 | 2.89 |
| `rfd4_design_006_of3_refolded` | 150 | G43 | 3.94 | G22 | 3.40 |
| `kemp_gln12_seq2_model1` | 136 | A49 | 3.78 | A24 | 3.69 |
| `1_model_4_id62` (round 1 best) | 149 | A8 | 3.24 | A38 | 3.57 |

---

## Findings

### 1. The base contact is reachable; the donor is not

Of 663 predictions, **13** hit the base within 3.5 Å *and* its approach angle
within 25°. Only **2** hit the donor on the same criteria. **Zero** hit both.
Among the 13 good-base cases the donor sits at a median of 7.79 Å.

The best single result was `1_model_4_id62` at base 3.24 Å / donor 3.57 Å —
and its donor *angle* is 81° against a 158° reference, so the amide N–H points
away from the ligand. Proximity, not a hydrogen bond.

### 2. The failure is ligand recruitment, not folding

`3_model_1` reproduces its fold at 1.77 Å CA-RMSD with pLDDT 0.838, yet Boltz
docks 6NT **11.96 Å** from the pinned pose. `1_model_6` likewise: CA-RMSD
1.45–2.19 Å, ligand off by 13.7 Å. The scaffolds fold as designed and then fail
to hold the substrate where intended.

### 3. LigandMPNN's confidence predicts nothing

Across 311 predictions:

```
lig_conf vs base distance     r = +0.048
lig_conf vs donor distance    r = -0.058
lig_conf vs complex_plddt     r = +0.025
lig_conf vs CA-RMSD           r = +0.091
within design:  |r| <= 0.17 everywhere
```

Good-base predictions are indistinguishable from the rest (median 0.438 vs
0.403; their *maximum* is lower than the overall maximum). This invalidated the
"sample 256 deep, filter cheaply, then predict" strategy that justified round
1's sampling depth.

### 4. Input-geometry metrics never predicted outcomes

Tip bond angles, junction bond lengths and χ strain were each tested as
predictors and each failed. `1_model_17` had the cleanest input chemistry of all
80 designs (0.018 Å bond deviation, 14° worst angle) and produced among the
worst catalytic geometry. `3_model_1` had some of the worst input chemistry and
the smallest catalytic displacement under packing.

### 5. Sequence diversity is concentrated, not uniform

Zero cross-backbone pairs reach 70% identity (ceiling 45.8%). Within backbones,
redundancy sits almost entirely in the two highest-alanine designs: `1_model_6`
(38.3% Ala) and `3_model_1` (31.7%) account for 98.6% of all ≥70% pairs. Their
256 samples compress to 29 and 35 distinct sequences respectively.

---

## Corrections made during the campaign

Recorded because several changed conclusions:

| issue | consequence |
|---|---|
| `--no_kernels` assumed optional | **Required** — `cuequivariance_torch` is absent; a pilot failed on all 10 without it |
| `msa:` field omitted | Does **not** mean single-sequence; the input is silently dropped with no output. Must be explicit `msa: empty` |
| Donor angle via `CG` | Wrong atom — `CD` is bonded to `NE2`. Reference is 158°, not the 121° first reported |
| CIF coordinates read positionally | Broke on RFdiffusion4's 62-column format, where the last column is `id`, not `model_num`. Produced plausible but meaningless distances. Fixed by `src/cifread.py`, which reads the loop header |
| 4-way GPU parallelism assumed | AWS vCPU quota is 8; one g6.2xlarge consumes it entirely |
| mmseqs at 12 GB / 4 CPUs | Database is ~180 GB; the working kemp-4 config uses 30 GB / 8 CPUs |

---

## Key scripts

All under `pdh-team:/mnt/work/projects/Dejvid/src/`:

| script | purpose |
|---|---|
| `cifread.py` | header-aware mmCIF reader; the three CIF flavours put `Cartn_x` at index 10, 17 and 57 |
| `sanitize_rfd3_pdb.py` | renumbers atom serials so biotite can parse the theozyme |
| `prepare_ligandmpnn.py` | RFD3 → LigandMPNN bridge, motif from `diffused_index_map` |
| `prepare_external.py` | same for externally supplied structures, motif detected geometrically |
| `detect_motif.py` | standalone motif detection with ambiguity reporting |
| `check_packed_geometry.py` | catalytic geometry in packed structures + theozyme-exact variants |
| `measure_boltz.py` | distances, confidences, CA-RMSD, ligand displacement |
| `scan_angle.py` | approach angles against the 127.3° / 158.1° references |
| `iter2_filter_and_msa.py` | novelty filter + MSA-backed Boltz inputs |

---

## Open questions

- Is the donor geometry physically achievable in these scaffolds at all? A Boltz
  run with an explicit `constraints: contact` on the donor would answer this
  directly, rather than inferring it from cases where it never arose.
- Does the round-2 chain (lower-temperature LigandMPNN, SolubleMPNN surface
  pass, MSA-backed prediction) preserve the four input structures that already
  satisfy both contacts?
- With `ligand_confidence` ruled out, is there any cheap filter that predicts
  Boltz success — or must every design be predicted?
