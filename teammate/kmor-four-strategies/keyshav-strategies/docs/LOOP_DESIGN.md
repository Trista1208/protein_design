# Design notes for the loop-based pipeline — 22 Nov 2026, 03:00

Written for whoever builds the next pipeline, including a future session with
no memory of tonight. Supplements `LESSONS.md`; **where the two disagree, this
file wins**, because §1 below corrects two things in LESSONS that were measured
wrong.

Everything is marked **[measured]**, **[published]** (with its source) or
**[open]**. Nothing here is inferred and presented as fact. Read §1 and §6
first.

---

## 1. Corrections to LESSONS.md — read before trusting its §6

### 1.1 `rfd4_design_045` has FOUR zeroed terms, not one **[measured]**

LESSONS §6 says *"Everything is solved except the oxyanion donor orientation"*
and prints `oxy 0.994`, `orient_base 0.364`, `sd_base 0.0759`. Those numbers
are **not** what the organisers' protocol produces. That table mixes two
predictors: the geometry came from the OpenFold3 model, only `ipde` came from
the Boltz refold — the exact incomparability LESSONS §3.2 warns about, landing
inside the headline design's own scorecard.

The actual 5-seed Boltz refold (`outputs/refold_d045/verdict.json`,
`msa: empty`, the organisers' "5 replicates per sequence"):

| term | LESSONS §6 claims | actually measured |
|---|---|---|
| base | 0.687 | **0.182** |
| oxy | **0.994** | **0.000** |
| ipde | 0.914 | 0.914 |
| orient_base | 0.364 | 0.000 → **see §1.2, an averaging bug** |
| orient_oxy | 0.000 | 0.000 |
| bury | 1.000 | 1.000 |
| cons | (sd 0.0759) | **0.000** (sd_base **1.047**) |

**The oxyanion contact never forms under Boltz.** Per-seed `d_oxy`:

```
seed 1  3.964   OUT        window is 2.0 - 3.6 A
seed 2  3.654   OUT
seed 3  3.816   OUT
seed 4  9.654   OUT        ligand largely lost: d_base 5.272 too
seed 5  4.327   OUT
```

Not one seed lands inside the window. So the real blockers are **three**, not
one: oxyanion *distance*, oxyanion *orientation*, and *consistency*. The
`sd_base 0.0759` in LESSONS was computed "excluding one outlier seed" — the
official `cons` term permits no such exclusion, and with all five seeds it is
1.047, which zeroes `cons` on its own.

**Consequence for planning:** "we are 3.4° away" is false. That was the most
load-bearing claim of the previous session and it does not survive measurement.

### 1.2 Dihedrals were averaged linearly across replicates — FIXED **[measured]**

`official.combine_replicates()` took an arithmetic mean of a circular quantity.
Base dihedral across the five seeds:

```
354.8   350.2   349.6   83.4   355.6
linear mean   = 298.7 deg  -> 230-310 DEAD ZONE -> s_dih 0.000
circular mean =   6.7 deg  -> 310-50 band, near ideal -> s_dih 0.889
```

Four of five seeds are in band and close to ideal; the linear mean put the
aggregate in the dead zone and zeroed `orient_base` (weight 0.12) as a pure
artefact. Fixed in `official.py` with `_circular_mean()`; `d_base`, `d_oxy`,
the two angles and burial stay linear. **Averaging 350 and 10 must give 0,
never 180.** This is LESSONS §3.11's "dihedral conventions bite", one layer up.

After the fix the oxy dihedral aggregate is unchanged at 237.2° (all five seeds
were already in the same dead zone, so there was nothing for the bug to hide).

### 1.3 Better still: aggregate the SUBSCORES, not the angles **[open]**

Averaging any geometry and then scoring once is a second-best. The handover's
own rule (4.3.7, echoed in LESSONS §5.2) is *rank on median plus dispersion or
worst case, never the best seed*. Scoring each replicate and aggregating the
**scores** would need no circular-mean reasoning at all and would expose the
seed-4 blow-up instead of smearing it. Not done tonight; worth doing.

---

## 2. What is actually installed — the constraint that kills three pipelines

Verified tonight, not assumed:

| Tool | Reality | Consequence |
|---|---|---|
| **Rosetta** | **BUILT AND WORKING.** `rosetta_scripts`, `relax` and `match` all compiled 21 Nov 13:50-13:51, all run, database present. No PyRosetta. | RosettaScripts / FastDesign / RosettaRelax / enzdes constraints are **available**. pSUFER additionally needs snakemake (pip, minutes). |
| **APNovo** | Code + env install clean, `alphafold3.structure` imports, **no weights** (gated request) | Generator with no model. Same dead end as RFD4 locally. |
| **RFD4** | env exists, no binary, no checkpoint | Team runs it elsewhere |
| **OpenMM stack** | openmm 8.6.1, openff 0.19.0, openmmforcefields 0.16.0, parmed 4.3.1, rdkit 2026.09.2, **CUDA live** | Works. Usable for restrained minimisation / preorganisation checks |
| **Foldseek + MMseqs** | installed, wired into the novelty gate | Works |
| **Boltz-2** | `envs/her2_boltz2/bin/boltz` | The only evaluator that counts |

**CORRECTION, 22 Nov ~04:45.** An earlier version of this section claimed
Rosetta was not built and used that to rule out Pipelines 1, 2 and 4. That was
wrong, and the error is worth recording because of how it happened: the check
was `ls -la source/bin/ | head -5`, and `head -5` cut the listing after the two
alphabetically-first `match` entries. `relax` and `rosetta_scripts` were on the
next lines. A truncated command was then reported as a verified absence and
written into this document as fact.

LESSONS §2.1's rule generalises: *if a check would condemn something you
expect to work, suspect the check.* Here it condemned three pipelines.
**Verify a negative before building a plan on it.**

---

## 3. From denovoKemp (Listov et al., Nature 2025) **[published]**

`PMC12310539` *is* this paper (Nature 643, doi 10.1038/s41586-025-09136-2).
Both links in the brief point at the same work. nature.com and PMC block
automated fetches; the full text is readable at
`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12310539/fullTextXML`.

### 3.1 They deliberately have NO oxyanion donor

Verbatim from the paper:

> *"Typical design studies also introduced a polar interaction with the
> isoxazole oxygen to stabilize the developing negative charge in the
> transition state. We excluded this requirement from our theozyme because a
> water molecule can satisfy it, and a misplaced polar group could reduce
> reactivity by lowering the pKa of the catalytic base."*

Their theozyme is **two** elements: Asp/Glu base + aromatic (Trp/Phe/Tyr)
π-stack onto the benzene centroid. No Gln. The best design reached
k_cat/K_M ≈ 12,700, and **Phe113Leu — deleting the designed aromatic — raised
it to ≈123,000**. Both designed elements turned out dispensable or absent.

**What this means for us:** the official metric ties 0.28 of the score
(`oxy` 0.20 + `orient_oxy` 0.08) to a Gln donor that the highest-performing
published Kemp eliminase deliberately omitted — and whose coordinates in our
own input were *modelled back in*, not observed (handover §1.1: B-factor 0.00
on all 17 atoms). We cannot opt out, because the geometric mean zeroes us. But
treat the oxyanion as a **scoring obligation, not a catalytic one**: spend
exactly enough to clear zero, and no more.

### 3.2 Their only enforced dihedral is 180° ± 7° — and it does NOT transfer

The entire repo has one constraint file, `2.MatchDesign/1.Match/theozime.cst`
(sic). In the base block only `torsion_A` has a non-zero force constant:

```
CONSTRAINT:: distanceAB:   2.93  0.2   100.  1    1     <- "1" = COVALENT flag
CONSTRAINT::    angle_A:  127.6  5.0   80.0  360. 1
CONSTRAINT::    angle_B:  120.0  5.0   80.0  360. 1
CONSTRAINT::  torsion_A:  180.0  7.0   70.0  360. 1     <- the only live dihedral
CONSTRAINT:: torsion_AB:    0.0  0.0   0.00    5. 1     <- k = 0, inert
CONSTRAINT::  torsion_B:  180.0 10.0   00.0  360. 1     <- k = 0, inert
```

Columns are `optimal / tolerance / penalty k / periodicity / n samples`
(confirmed by the paper's Methods and Supplementary Table S6).

This **independently corroborates `DIHEDRAL_IDEAL = 180.0`** in `official.py`.

**Do not port the ±7°.** Their `torsion_A` is `N1–C1–C7–OD` on the *substrate*,
where Asp contacts ring **carbon** C3; ours is `CG–OD2–N3–N2` on **6NT**, where
Asp contacts ring **nitrogen** N3. Different ligand, different contact atom,
different quadruple. Porting the number would repeat the cross-frame error that
the ORI trap and LESSONS §2.2 both warn about. It confirms the centre of the
band; it says nothing about our edges. **§1.8 Q3 stays open.**

Also note the base contact there is 2.93 Å declared *covalent* to fake TS
geometry — another reason the numbers are not comparable to our TSA distances.

### 3.3 Numeric gates worth adopting **[published]**

From the Methods (none of these are in the repo files):

- **Discard pLDDT < 90.** → *we are using **≥ 0.95***, tightened by the team.
  d045's seeds are 0.963–0.967 and survive; yesterday's local designs at
  0.41–0.47 correctly vanish before geometry is ever scored.
- **Preorganisation: discard if catalytic-base RMSD > 1.2 Å** between
  ligand-bound and ligand-free models. The slide calls `cons` *"a proxy for
  preorganization"*; this is the direct measurement. One extra apo fold per
  design → finalists only.
- **Design shell: 8 Å** around ligand + catalytic residues. Use this instead of
  a hand-picked radius.
- Catalytic residues **fixed for the first two design rounds, then released**.
- AF2-vs-model RMSD < 1 Å for their final series.
- ΔΔG thresholds (+1 / +3 R.e.u.) — Rosetta-only, not available to us.

### 3.4 Their selection objective **[published]**

`3.FilterDesigns/design_choose.ipynb`: a **sigmoid per term with μ and σ taken
from the surviving population** — data-adaptive, not fixed thresholds — then

```
if_vdw · total_score/len · ligand_fa_sol · totalS_rank
        · (dihedral_cst · angle_cst · atom_pair_cst)^(1/2)
```

keeping the top 1000. Two things to steal: population-relative normalisation
(which our percentile ranking already does), and the **square root on the
geometry-constraint block** — their explicit judgement that constraint scores
are noisier than energies and should be down-weighted.

---

## 4. From AI.zymes (Merlicek et al., 2025) **[published]**

bioRxiv 10.1101/2025.01.18.633707. Optimises Kemp elimination in KSI D38N.
Code not yet released. The **selection machinery** is the transferable part and
needs no Rosetta:

- **Boltzmann selection over the entire pool**, not discrete rounds —
  *"operating on the entire design pool rather than using rounds-based
  approaches"*, continuously running a fixed number of designs.
- **Fitness = "potential" = mean of a variant's own score and its
  first-generation descendants'** — forward-looking, keeps good lineages alive.
- **Standardise each objective, then average.** Never a product.
- **Annealed temperature**, exploration → exploitation. Values are in Table S1,
  which was not readable — **pick your own schedule and say so**.
- **Mutation load capped**, with an explicit bias toward the parent sequence
  (Rosetta: a mutation penalty weight; ProteinMPNN: a bias to the input
  sequence). For us: bias LigandMPNN to the parent. Same instruction as
  LESSONS §7.6.
- **Geometry is a hard filter *before* selection, not a term in the score.**
  Matches LESSONS §4.5 and how our gates already split.

**Calibration, and it is sobering:** 104 designs made → 7 tested → **1
improved**, for 7.7× — starting from an enzyme that already worked. Do not plan
a sprint around an expected breakthrough.

### 4.1 Why a product is fatal for us but fine for them

LESSONS §4.3 says "percentile-rank scalarisation, never a product". denovoKemp
*does* use a product — but of **sigmoids, which never reach 0**, so nothing
annihilates. The sharper rule:

> A product is safe when no term can be exactly zero. The official metric is
> lethal precisely because its terms can be, and are.

Our dense objective (`kemp/objective.py`) is on the right side of that line.

---

## 5. From Theo's RFD4 protocol **[measured, his report]**

- Backbones: `tsr_k` 8/16/24, 2 seeds × 25 = 150 backbones.
- Inverse folding: local-latent `tsr_k` 1/2/4, 20 sequences per setting.
- **One sequence per setting selected by clashes, then motif RMSD.** Adopt
  this — it is the concrete fix for `pipeline.py:410`'s `seqs[0]`, which throws
  away 7 of 8 sequences spanning `d_base` 2.557–9.329 Å.
- **Refolded "with 6NT and ColabFold MSAs"** — a different protocol from our
  `msa: empty`. See §6.

---

## 6. The MSA question — the one live experiment

Two protocols were in circulation for the same design lineage and **nothing
recorded which produced which number**:

- Theo: ColabFold MSAs
- `scripts/refold_d045.py:50`: `msa: empty`

This matters because LESSONS §3.1 shows the MSA setting moved HG3.17's `d_oxy`
from **12.698 → 2.952 Å** — the same term that is zeroed on every d045 seed
(§1.1). LESSONS §3.2 counters that ColabFold returns depth 1 for a de novo
design, making the two identical in effect; that was never checked on this
sequence. The organisers' slide 7 lists *"MSA depth for de novo design can
impact folding confidence"* as a known limitation.

`scripts/msa_test_d045.py` runs it: same sequence, same seeds, same settings,
**MSA the only variable**, and it records the actual MSA depth. Note that
single-sequence pLDDT here is already 0.963–0.967, so this is *not* a
foldability test — it asks only whether an MSA changes the Gln rotamer and the
oxyanion distance.

### 6.1 RESULT — settled **[measured, 22 Nov 03:15]**

All 5 seeds, `--use_msa_server` vs `msa: empty`:

```
MSA depth = 1
  (outputs/msa_test_d045/seed1/.../msa/input_unpaired_tmp_env/
   bfd.mgnify30.metaeuk30.smag30.a3m)

with MSA   dihedral_oxy: [239.4, 234.4, 239.0, 239.8, 233.4]
no MSA     dihedral_oxy: [239.4, 234.4, 239.0, 239.8, 233.4]
band gaps:               [  9.4,   4.4,   9.0,   9.8,   3.4]
seeds in band: 0/5
```

**The MSA server returns the query and nothing else.** With depth 1,
`--use_msa_server` and `msa: empty` are literally the same input, and every
seed is identical to three decimals on `d_base`, `d_oxy` and the dihedral.
**LESSONS §3.2 is confirmed on this sequence.**

Three things follow:

1. The 237° rotamer and the 3.65–9.65 Å oxyanion distance are properties of
   the **sequence**, not of the fold protocol. Only a sequence change can move
   them. The loop's premise survives.
2. Theo's "refolded with ColabFold MSAs" and our `msa: empty` are, for de novo
   designs, **the same protocol wearing two labels**. The scores are
   comparable after all — but only because depth came back 1. Keep stamping
   `msa_mode` and `msa_depth`: for any sequence with real homologues (a
   natural control such as HG3.17) they are emphatically *not* the same, which
   is exactly how the two became incomparable before.
3. There is no cheap protocol fix waiting. Do not spend more time here.

### 6.2 d045's true state, after the circular-mean fix **[measured]**

The same run, scored through the corrected `combine_replicates`:

```
base=0.182  oxy=0.000  ipde=0.914  orient_base=0.407
orient_oxy=0.000  bury=1.000  cons=0.000
official_score = 0.0   limiting = oxy
```

`orient_base` recovered from 0.000 to **0.407** — that is §1.2's bug, confirmed
on a real run. Three zeros remain: **oxy, orient_oxy, cons**.

**The three zeros are not three problems.** Look at where they come from:

| seed | d_base | d_oxy | dihedral_oxy |
|---|---|---|---|
| 1 | 3.023 | 3.964 | 239.4 |
| 2 | 2.864 | **3.654** | **234.4** |
| 3 | 2.975 | 3.816 | 239.0 |
| 4 | **5.272** | 9.654 | 239.8 |
| 5 | 2.882 | 4.327 | 233.4 |

- `d_base` over the four good seeds is 2.864–3.023, **sd ≈ 0.075** — essentially
  HG3.17's 0.068. Consistency is excellent four times in five and catastrophic
  once. `cons` is zeroed by **seed 4 alone**, not by general sloppiness.
- `oxy` and `orient_oxy` are both properties of the **same Gln46 side-chain
  tip**: `d_oxy` is NE2···N1, `orient_oxy` is the dihedral CD–NE2···N1–N2. One
  chi2/chi3 change moves both.
- The misses are small and systematic, not noisy: best seed needs
  **0.054 Å** on distance (3.654 vs the 3.6 edge) and **3.4°** on the dihedral.

So the design objective is **one rotamer**, and fixing it addresses two of the
three zeros at once. The previous session was right in kind and wrong in
magnitude: not "3.4° on one term", but "one Gln rotamer, worth 0.28 of the
score", plus a separate seed-4 stability problem.

---

## 6.3 THE SUPPLIED THEOZYME'S Gln ROTAMER IS WRONG **[measured, 05:10]**

The most consequential measurement of the night. Oxyanion torsion
(`CD-NE2...N1-N2`, the scored `orient_oxy` quantity):

| structure | Gln50 | torsion | `s_dih` |
|---|---|---|---|
| **4BS0** crystal | real, refined | **2.74°** | **0.954** |
| **7K4Z** crystal | real, refined | **355.58°** | **0.926** |
| **supplied theozyme** | **modelled, B=0.00** | **50.90°** | **0.000** |

Both crystals sit in the **310-50 band, near 0/360**, and score ~0.94. The
supplied theozyme's modelled Gln50 sits **~48° away from the crystallographic
rotamer**, lands 0.9° outside that band, and is the only one of the three that
scores zero.

The band definition is therefore NOT the problem - it is the one thing here
that is validated. LESSONS §2.1 says a check that condemns the reference is
wrong; this check *exonerates* both references and condemns only the modelled
input.

**The handover's §1.1 claim that the modelled Gln50 "DOES match the
crystallographic Gln50 in 4BS0/7K4Z - confirmed" is wrong on the torsion.** It
may match on distance (3.33 vs 3.10 Å); it does not match on the rotamer, and
the rotamer is what `orient_oxy` scores.

### Why this reframes the whole oxyanion effort

Every design scaffolded on this theozyme inherits a donor torsion that cannot
score. d045 sits at 237° because it was built around a motif whose own donor
torsion is out of band.

Worse, the "we are 3.4° from the band edge" framing was aiming at a nearly
worthless target. The taper runs from an ideal of 180° with half-width 60:

```
reaching the NEAREST edge, 230 deg  ->  s_dih = 1 - 50/60 = 0.167
reaching the CRYSTAL rotamer, ~0 deg ->  s_dih ~ 0.99
```

So the entire "nudge it 3.4° into band" campaign was worth **at most 0.167** on
`orient_oxy`, while the crystallographically correct answer - the *other* band,
near 0° - is worth ~1.0. We were optimising toward the wrong one of the two
allowed bands because it happened to be nearer.

### What follows

1. **Target the crystal rotamer (~0°/359°), not the nearest band edge.**
   `outputs/rosetta/kemp_theozyme.cst` now does this: torsion targets 339.9°
   (base) and 359.2° (oxy), the circular means of 4BS0 and 7K4Z, rather than
   the 180° an earlier draft assumed.
2. **Consider re-scaffolding from a corrected theozyme.** Replacing the
   modelled Gln50 rotamer with 4BS0's would give every future design a donor
   geometry that can reach ~1.0 instead of ~0.167. `scripts/fix_oxyanion.py`
   already exists but re-optimises toward a 2.4 Å distance target, not toward
   the crystal torsion - it needs redirecting.
3. The base torsion is in band everywhere (328-341°, s 0.47-0.68) and is not
   the problem.

---

## 6.4 THE SUPPLIED DATA FILES — read at 05:20, should have been read first **[measured]**

Three files sat in `/home/mandalorian/Hackathons/` all night and were never
opened. They contain literature-derived filters for *this reaction*, and one
of them directly refutes a route I spent hours building.

### 6.4.1 Arg127 salt-bridges the catalytic base in EVERY design

Privett 2012 / Alexandrova 2008 hard filter **H3**: no charged residue within
4 Å of the catalytic base. Measured with `kemp/litfilters.py`:

| structure | charged residues within 4 Å of the base |
|---|---|
| **4BS0 crystal** | **0** |
| **7K4Z crystal** | **0** |
| d045 parent | **2** — ARG127.NE @3.49 Å, ARG127.NH2 @3.68 Å |
| ins080 | **2** — ARG127.NH1 @3.26 Å, ARG127.NH2 @3.98 Å |

Both positive controls pass; every design of ours fails. This is the KE16
failure mode verbatim: *"Lys201–Asp48 salt bridge in the apo form competes
with the substrate."* And it has a fix with a measured payoff — in KE07 the
equivalent Lys222–base salt bridge was removed and **kcat/KM improved
2.5-fold**.

**Independent corroboration:** the untargeted Rosetta enzdes run, which knew
nothing about this filter, chose **R127I, R127L and R127M** unprompted in all
three of its designs. Two unrelated methods converging on one residue.

`R127X` is therefore the highest-priority proposal class in the campaign.

### 6.4.2 Rosetta score does NOT predict activity — on this exact reaction

`geometrie_kemp_risso2020_v2.xlsx`, sheet "Dati paper", section C,
*"Parametri che NON hanno funzionato come filtro di attività"*:

| parameter | result |
|---|---|
| **Rosetta score (FuncLib)** | **No correlation with activity.** Good for stability ranking, not for catalysis |
| pKa of catalytic base (PROPKA) | no significant correlation |
| Reorganization energy | no significant correlation |
| Stability (Tm) | all 20 variants stable, activity spans 2 orders of magnitude — *"use as a gate, not as a predictor"* |
| Per-residue electrostatics | small, <2 kcal/mol each |

Risso et al. 2020 (Chem. Sci. 11, 6134) tested this on Kemp elimination of
5-nitrobenzisoxazole — our reaction. It does not invalidate using Rosetta to
*propose* geometry, but it rules out ranking candidates by Rosetta energy, and
it retroactively explains why the enzdes designs refolded badly. It also means
the FieldTools/electrostatics branch was correctly deprioritised.

### 6.4.3 The geometry that DOES predict activity

Risso G1/G2, from EVB trajectories:

```
G1  D–A distance   WT 2.64 A   range 2.64-4.25   soft <=3.0, hard >4.0
                   R^2 0.84 vs calculated dG‡, -0.71 vs log(kcat/KM)
G2  D–H...A angle  WT 166.5    range 129.8-167.1  soft >=150, hard <135
                   R^2 0.71 vs log(kcat/KM)
```

**We were not measuring G2 at all.** It is the donor–hydrogen···acceptor
angle, near-linear when good — a different triple from the organisers'
`angle2` (centred on the ligand acceptor, target 120°). Both now computed;
`litfilters.py` keeps them separate and flags G2 as derived, since Boltz
output has no hydrogens and the H must be reconstructed.

Validation against the controls is reassuring:

| | G1 | G2 | base rel. SASA | S2 apolar |
|---|---|---|---|---|
| 4BS0 | 2.532 | **167.16** | 2.22% | 0.727 |
| 7K4Z | 2.515 | **164.09** | 3.96% | 0.727 |
| Risso WT | 2.64 | 166.5 | — | — |
| d045 | 2.864 | 139.44 | 3.40% | 0.556 |
| **ins080** | **2.644** | **150.93** | 3.33% | 0.500 |

ins080's G1 of 2.644 is essentially the Risso WT value, and its G2 clears the
soft threshold. Its remaining literature defects are H3 (Arg127) and S2
(apolar fraction 0.500 against the controls' 0.727).

### 6.4.4 Other filters now implementable

From `filtri_backbone_v2.xlsx` — hard: H1 base–substrate H ≤3 Å, H2 base
relative SASA ≤10%, H3 above, H4 no clashes, H5 catalytic residues not on
flexible loops. Soft: S2 apolar fraction ≥0.60, S4 apo-vs-holo site RMSD
≤1 Å, S7 donor contacts ≥6 within 4.5 Å, S9/S10 solvation terms.

`filters hackaton.docx` adds Bhattacharya 2022, including the point that
matches our own result: **"the triple mutant was 3× more active than predicted
from the three single mutations"** — combinations beat additivity, which is
exactly what F42Y+H103Y did.

**Caveat the sheets state themselves:** the Risso thresholds are
*"mie ipotesi dalle letture di Fig. 9: da tarare"* — hypotheses read off a
figure, to be calibrated, averaged over EVB trajectories rather than measured
on one predicted model. Treat a literature FAIL as a flag, never as an
override of the organisers' metric.

---

## 6.5 THE MECHANISM: a Tyr103–Gln46 hydrogen bond gates the score **[measured]**

`ins080` scores 0.743 on seed 1 and 0.000 on seeds 2–5. Comparing the Gln46
environment across those five structures isolates a single discriminator:

```
distance to the Gln46 side-chain tip
              GOOD(48.2 deg)   bad    bad    bad    bad   mean|diff|
  TYR103           2.36        7.43   7.01   6.70   6.89      4.65
  PHE77            9.41       10.05   9.64  10.79  13.60      1.62
```

Tyr103 is three times more discriminating than anything else, and the contact
is specific:

```
seed1 (scores 0.743)   Tyr103.OH ... Gln46.NE2   2.36 A
seed2 (scores 0.000)   no Tyr103-Gln46 contact under 5 A
seed4 (scores 0.000)   no Tyr103-Gln46 contact under 5 A
```

**The Tyr103 hydroxyl hydrogen-bonds the Gln46 amide and pins the rotamer in
the scoring geometry.** When that bond is absent the Gln46 side chain is free
and flips into a dead zone.

> **CORRECTED 06:25 — the residue was too specific, the principle holds.**
> Designs from library 1 reach the band with **Tyr103 6.7-6.9 A away**; there
> the pinner is **Tyr42**, the tyrosine introduced by the F42Y we already
> carried. Measured over 50 bound+measurable designs:
>
> | | nearest H-bond donor to the Gln46 amide |
> |---|---|
> | **in band (10)** | median **2.34 A**, range 2.29-2.72 |
> | out of band (40) | median **3.34 A**, range 2.38-3.56 |
>
> All ten in-band designs are pinned by **Tyr42.OH or Tyr103.OH**. So the rule
> is *some* H-bond-capable atom reaching the Gln46 amide, not one particular
> residue supplying it. It is **necessary but not sufficient**:
> `L_B_anchor_M45F` has a pinner at 2.38 A and is still out of band, so rank
> on it and explain with it, never gate on it.
>
> Implemented as `litfilters.gln_pinner()`, computed on every 1-seed screen -
> a cheap proxy for something that otherwise costs five folds to establish.

### This explains the F/Y pair we already had

| design | residue 103 | hydroxyl | replicates scoring |
|---|---|---|---|
| `ins079` | **Phe** | no | **0/5** |
| `ins080` | **Tyr** | **yes** | **1/5, best 0.743** |

Those two differ at exactly one atom — the phenol oxygen. It is not a
correlation across a noisy batch; it is a single-atom causal difference with a
measured 2.36 A contact behind it.

### What lets Tyr103 swing away

Every residue beyond its immediate neighbours shifts by a near-uniform ~6.5 A
between the good and bad states, which is the signature of the whole Tyr103
side chain rotating as a unit rather than the pocket rearranging. Its only
persistent contacts across all five states are **Met45 and Leu49**.

So the design target is now specific: **make the Tyr103-OH···Gln46-NE2 bond
permanent**, by biasing the Tyr103 chi angles through the bulk of Met45 and
Leu49, its anchor residues. Both are already in the campaign's scan set.

This also predicts which substitutions at 103 should work: donors that can
reach Gln46 (Tyr, Trp, and possibly Thr/Asn if the geometry allows) should
score; Phe, Leu, Ile and Ala should not. That is a falsifiable prediction the
campaign is testing.

---

## 6.6 TWO FAILURE MODES, not one **[measured 06:35]**

`L_B_anchor_M45Y` scored **0.787 on its screening seed** and **0/5 on five**.
That single-seed number was a candidate, not a result, and it did not survive.
But the way it failed separates two problems that had been conflated:

```
M45Y, 5 seeds        d_base   d_oxy   dihedral
  seed1               2.802   3.159     40.4   <- in window AND in band
  seed2               3.127   3.278     62.1   <- dead zone
  seed3               2.947   3.149     35.3   <- in window AND in band
  seed4               4.822   7.154    115.1   <- LIGAND LOST
  seed5               2.688   3.150     60.7   <- dead zone
```

M45Y reaches full scoring geometry on **2 of 5** seeds; `ins080` manages
**1 of 5**. On the rotamer problem M45Y is the better design. It scores zero
anyway because **seed 4 loses the ligand**, taking `sd_base` to 0.879 and
zeroing `cons` — and `cons` is a whole-batch term, so one bad replicate
zeroes *every* replicate.

Excluding that seed (diagnostic only; the official term permits no
exclusion): per-replicate `[0.761, 0.000, 0.720, 0.000]`, sd_base 0.19.

### The two modes

| mode | signature | seen in |
|---|---|---|
| **A. Rotamer flip** | distances fine, oxyanion dihedral in a dead zone | ins080 (4/5 seeds), M45Y (2/5) |
| **B. Ligand loss** | one seed with d_base > 4 Å; `cons` collapses | M45Y, ins081 (sd 1.45), d045 (seed 4) |

They need different fixes. Mode A is the Gln46 rotamer and responds to the
pinner. Mode B is whether the pocket holds the ligand at all, and no amount of
rotamer engineering touches it.

**Consequence for selection:** rank on *how many seeds reach full scoring
geometry*, not on the best single seed. By that measure M45Y (2/5) is ahead of
ins080 (1/5) even though its official score is lower — and a design that fixes
mode B on top of M45Y's geometry would score on 2-3 replicates rather than 0.

---

## 6.7 BLOCK A: INCONCLUSIVE, and a verdict rule that overstated it **[measured 07:05]**

The falsifiable test of the donor hypothesis, on 225 screened designs:

```
BLOCK A, bound + measurable only
  donors  (Y W T S N Q H E D K R at 103)  n=9  median gap 27.6 deg  in band 0
  control (F L I V A M G at 103)          n=6  median gap 30.0 deg  in band 0
  donor gaps  23.4 25.5 26.7 27.5 27.6 29.3 29.3 31.6 33.1
  control     28.6 28.9 29.3 30.6 32.0 32.4
```

**Neither arm put a single design in band.** The distributions overlap almost
entirely. Substituting residue 103 alone - donor or not - does not decide the
oxyanion dihedral.

### The verdict rule was wrong before the verdict was

`supervisor.py` reported **"SUPPORTED (by gap, not band)"**. It fired on
`donor_best_gap < control_best_gap - 5`, i.e. 23.4 against 28.6: a best-of-9
compared with a best-of-6. Taking the extreme of the larger sample and calling
the difference a result is the same error as quoting a one-seed score. The
rule now requires a **median** separation of 10 deg, and reports
INCONCLUSIVE when neither arm reaches the band.

### What this does NOT overturn

The pinner correlation, now on all 146 bound+measurable designs:

| | nearest H-bond donor to the Gln46 amide |
|---|---|
| **in band (15)** | median **2.35 A** |
| out of band (131) | median **3.38 A** |

In-band pinners: **Tyr42.OH x7, Tyr103.OH x8**. The correlation is intact and
stronger than when first measured on 50.

### Reconciling the two

A close pinner accompanies an in-band dihedral, but **putting a donor at
position 103 does not create one**. The contact is a consequence of the whole
pocket's packing, not something installed by one substitution - which is why
the 15 in-band designs came from the anchor and shell positions (M45, L49,
F26, V84), not from the donor scan.

**Design consequence:** stop treating "has a donor at 103" as the lever. The
lever is whatever packing brings an existing tyrosine's hydroxyl to 2.35 A -
and that is a property of several residues at once, which is why pairs and
triples (library 2) are the right next move rather than more single
substitutions.

---

## 6.8 PAIRS WORK: 0.8308, and mode B is solved **[measured 08:00]**

Library 2 screened PAIRS of the single mutations that library 1 found. The
effect is the largest anything has shown:

| library | in band / bound+measurable | rate |
|---|---|---|
| 1 (single substitutions) | 15 / 146 | **10%** |
| 2 (pairs) | 23 / 49 | **47%** |

### The best design so far

`L2_pair_I23L_V84M` = ins080 + I23L + V84M, re-measured from its five
structures:

```
per-replicate totals  [0.8308, 0, 0, 0, 0]      best 0.8308 on 1/5
mean pLDDT 0.979      sd_base 0.0266            (HG3.17 reference 0.068)

seed1  d_base 2.540  d_oxy 3.041  dih  33.35   total 0.8437  limiting orient_base
seed2  d_base 2.504  d_oxy 3.046  dih  61.89   total 0.0000  limiting orient_oxy
seed3  d_base 2.577  d_oxy 2.998  dih  53.26   total 0.0000  limiting orient_oxy
seed4  d_base 2.550  d_oxy 3.032  dih 100.47   total 0.0000  limiting orient_oxy
seed5  d_base 2.556  d_oxy 4.129  dih 251.58   total 0.0000  limiting oxy
```

For scale: HG3.17's crystal structure scores 0.959 and the same enzyme folded
with an MSA scores 0.834. **0.8308 is level with the folded reference
enzyme.**

### Mode B is solved; only mode A remains

`sd_base` is **0.0266**, four times tighter than HG3.17's 0.068, and `d_base`
spans only 2.504-2.577 across all five seeds. The ligand-loss failure (6.6
mode B) that zeroed `cons` on M45Y and d045 is gone.

Everything that varies is the oxyanion dihedral: 33.4, 61.9, 53.3, 100.5,
251.6. And two of the four failures are NEAR misses - seed 3 at 53.26 misses
the 310-50 band edge by **3.26 deg**, seed 2 by 11.9. Three of five seeds sit
within 12 deg of scoring.

### What this says about the search

`I23L` appears in every one of library 2's best pairs and was unremarkable on
its own in library 1. A residue with no single-mutant signal became the
strongest component of the best pairs, which is the concrete form of
"combinations beat additivity" (Bhattacharya 2022 report a triple 3x better
than its singles predicted).

**Consequence:** triples of the library-2 winners are the next move, not more
pairs - and the target is now a 3.26 deg rotamer nudge on a design whose every
other term is already strong.

---

## 6.9 A DESIGN THAT SCORES ON MOST REPLICATES **[measured 08:27]**

`L_E_shell_V84I` = ins080 + **V84I**, a single substitution. Re-measured from
its five structures:

```
per-replicate  [0.8359, 0.7790, 0.8205, 0.7907, 0.0000]
best 0.8359    MEDIAN 0.7907    worst 0.0000    4 of 5 replicates score
mean pLDDT 0.9773    sd_base 0.0421   (HG3.17 reference 0.068)

seed1  d_base 2.611  d_oxy 3.086  dih  38.02   0.8517
seed2  d_base 2.577  d_oxy 3.193  dih  43.23   0.7937
seed3  d_base 2.609  d_oxy 3.108  dih  36.53   0.8360
seed4  d_base 2.507  d_oxy 3.143  dih  48.64   0.8057
seed5  d_base 2.582  d_oxy 4.098  dih 262.18   0.0000
```

**This is the first design with a non-zero MEDIAN.** Four of five replicates
put the oxyanion dihedral inside the 310-50 crystal band (36.5-48.6 deg), and
the fifth is the familiar single-seed ligand slip (d_oxy 4.098).

For scale: HG3.17's crystal structure scores 0.959 and the same enzyme folded
with an MSA scores 0.834.

| | best | median | reproducibility |
|---|---|---|---|
| `L_E_shell_V84I` | **0.8359** | **0.7907** | **4/5** |
| `L2_pair_I23L_V84M` | 0.8308 | 0.0000 | 1/5 |
| `ins080` | 0.7435 | 0.0000 | 1/5 |
| HG3.17 folded with MSA | 0.834 | — | — |

### A correction to 6.8's framing

6.8 called pairs "the effect". On 1-seed in-band RATE they clearly are
(47% against 10%). But the best CONFIRMED design is a **single** substitution,
and the pair that was 5-seeded scored on only 1 of 5. A high in-band rate at
one seed is not the same property as reproducing across five, and the two have
now visibly diverged.

Both V84I and V84M act at position 84, which did not appear in the pinner
analysis at all - its effect is on the ligand's own seating (`d_oxy` 3.09-3.19
on the four good seeds), not on the tyrosine contact.

### What is left

One failure mode, one seed in five: the ligand slips to d_oxy ~4.1 A and the
dihedral falls to ~262 deg. `cons` survives it (0.9159 worst-case) because
`d_base` barely moves, so this costs the median rather than zeroing
everything.

---

## 6.10 THREE DESIGNS SCORE ON EVERY REPLICATE **[measured 08:55]**

Re-measured from their structures, scored per replicate:

| design | per-replicate | best | **median** | worst | repro | sd_base |
|---|---|---|---|---|---|---|
| `L_E_shell_I23L` | 0.833 0.767 0.843 0.817 0.862 | 0.8617 | **0.8330** | 0.7671 | **5/5** | 0.0605 |
| `L2_pair_I23L_V84L` | 0.856 0.745 0.804 0.763 0.884 | 0.8839 | 0.8044 | **0.7448** | **5/5** | **0.0140** |
| `L2_pair_I23V_M45Y` | 0.811 0.698 0.873 0.789 0.887 | 0.8866 | 0.8106 | 0.6975 | **5/5** | 0.0719 |
| `L_E_shell_V84I` | 0.836 0.779 0.821 0.791 0.000 | 0.8359 | 0.7907 | 0.0000 | 4/5 | 0.0421 |

**No term is zeroed on any replicate** for the top three. Every one passes the
length, identity and TM gates and screens as NOVEL.

For scale: HG3.17's crystal structure scores **0.959**; the same enzyme folded
with an MSA scores **0.834**. `L_E_shell_I23L`'s median of 0.8330 is level
with the folded reference enzyme, and `L2_pair_I23L_V84L`'s worst replicate is
0.7448 - its *floor* is higher than anything this run produced before 08:00.

`L2_pair_I23L_V84L` also has `sd_base` **0.0140**, five times tighter than
HG3.17's 0.068.

### Which number to quote

By **median**, `L_E_shell_I23L` (0.8330) leads. By **worst case**,
`L2_pair_I23L_V84L` (0.7448) leads. By **best single replicate**,
`L2_pair_I23V_M45Y` (0.8866). All three are defensible; the median is the
honest headline and the worst case is the conservative one. Quote them with
their reproducibility, never alone.

### Position 23 is the lever

`I23L` or `I23V` appears in three of the four. Position 23 never appeared in
the pinner analysis and had no notable single-mutant signal in library 1's
donor scan - it emerged only once the shell block and the pairs were screened.
Position 84 and `M45Y` are the other recurring components.

This is the third mechanism found tonight, after the Gln46 pinner (6.5) and
the Arg127 salt bridge (6.4.1), and the only one that produced designs scoring
on all five replicates.

---

## 6.11 TRIPLES DID NOT STACK **[measured 09:48]**

Library 3 screened every 3- and 4-way combination of the components that had
produced a scoring design. The prediction was that if pairs beat singles,
triples should beat pairs. They did not.

Best of each order, by median, all re-measured from structures:

| order | design | repro | median | best | **worst** |
|---|---|---|---|---|---|
| single | `L_E_shell_I23L` | 5/5 | **0.8330** | 0.8617 | 0.7671 |
| **triple** | `L2_triple_I23L_M88F_V84L` | 5/5 | 0.8316 | 0.8748 | **0.7807** |
| pair | `L2_pair_V84I_I23V` | 5/5 | 0.8312 | 0.8536 | 0.7622 |

The best single, pair and triple are separated by **0.0018 in median** - a
dead heat. The triple does have the highest floor of anything measured
(worst replicate 0.7807), which is the one place the extra mutations show.

### What the progression actually looks like

```
singles  ->  pairs        10% -> 47% in-band at one seed, and the first
                          designs scoring on 5/5
pairs    ->  triples      no further gain in median; a slightly better floor
```

The large effect was getting from one mutation to two. Beyond that the pocket
appears to be at a local optimum that more substitutions of the same kind do
not improve. That is a useful stopping signal: further combinatorial depth on
these positions is not where the next gain lies.

### Where the remaining loss is

Every one of the seven 5/5 designs is limited by `orient_base` or
`orient_oxy` on its worst replicate, never by distance, burial, iPDE or
consistency. `d_base`, `bury` and `cons` are effectively solved. What is left
is the same rotamer question the run opened with - only now it costs a design
0.07 of median rather than zeroing it.

---

## 7. The loop, as specified

**Objective** — `kemp/objective.py`, built and tested tonight (36 tests):

- `dihedral_band_gap()`: degrees to the nearest sp2 band edge. The official
  sub-score is a flat 0 across the whole dead zone, so 237° and 290° tie at
  zero; the gap separates them (7.2 vs 20). This is what makes a batch of
  zeros rankable.
- `badness()` → continuous defect per official term, defined where the
  sub-scores are flat 0.
- `rank_batch()` → percentile per term across the batch, weighted by the
  official weights. Validated: on the real designs it orders
  HG3.17 crystal 0.860 > d045 0.500 > d00000 0.405 > HG3.17 no-MSA 0.235,
  where the official score ties the last three at 0.000.
- **`worst_term` ≠ the zero.** It is batch-relative and blind to whether a
  value is acceptable; on d045 it names `ipde` (sub-score 0.914, needs no
  work). Read `official.limiting_term` for the zero. Both are in the CSV.

**Gates, in order** (hard filters before any ranking — AI.zymes' ordering):

1. length ≤ 260, novelty (identity ≤ 0.70, TM ≤ 0.70) — existing
2. **pLDDT ≥ 0.95** (team's threshold; denovoKemp publishes 0.90)
3. motif preserved, chain integrity — existing `kemp/strain.py`
4. geometry hard filter against the theozyme's own geometry (trusted to
   ≤0.005 Å, LESSONS §5)
5. preorganisation: apo-vs-holo base RMSD < 1.2 Å — finalists only

**Search:**

- parent = a design that already folds (d045 or a successor), never a fresh
  sequence — LESSONS §7.6, HER2/Catalyst §4.1
- mutation shell 8 Å around ligand + catalytic residues, catalytic core frozen
- ≤ 2 edits per variant, LigandMPNN biased to the parent
- select among the 8 MPNN sequences by clashes → motif RMSD (Theo), never
  `seqs[0]`
- Boltzmann selection over the whole pool, annealed, descendant-averaged
  potential (AI.zymes)
- **a refold is the only evaluator.** Coordinate repacks do not survive
  (LESSONS §6: 0/5 seeds sampled any of the 360 clash-free in-band rotamers)
- cluster before picking the final 5, so the pool cannot collapse to one lineage

**Instrumentation, non-negotiable before a campaign launches:**

- degeneracy alarm `objective.audit()` / `assert_ok()` — fails on all-missing,
  >50% missing, ≤1 unique value, std < 1e-8. Note `skipped` is **not** `passed`
  (`status.json` now carries `degeneracy_judged` separately).
- **campaign early-abort checkpoint** — `objective.campaign_checkpoint()`,
  wired into `pipeline.scale()`. See §7.1.
- `protocol_stamp()` on every record before folding, so a score always carries
  its MSA mode; `comparable()` refuses a cross-protocol comparison.
- dihedrals and band gaps in `designs.csv` for every design and every replicate.

---

### 7.1 The degeneracy gate — decision, and why it is not where you'd expect

**It is not a launch gate.** A check in `launch_campaign.sh` runs when there
are zero scored designs, so it can only ever pass. That is theatre. The LESSONS
§2.5 failure is *only* visible once several designs have been scored and their
geometry turns out identical.

So it is an **early-abort checkpoint inside `pipeline.scale()`**, fired once,
after **10 scored designs** (`CAMPAIGN_CHECKPOINT_AFTER`). At the measured
~34 s/design that costs about six minutes of a 4h41m/500-design campaign.
`checkpoint_after=0` disables it — a legitimate choice for a 3-design pilot,
where no distribution exists to judge.

**Fatal vs advisory**, because conflating the two is how an alarm gets switched
off:

| finding | verdict | why |
|---|---|---|
| term constant across designs | **FATAL, abort** | the §2.5 signature: the pipeline is reporting its input |
| term has std < 1e-8 | **FATAL, abort** | same thing, through rounding |
| term entirely missing | **FATAL, abort** | that term can never be scored |
| term >50% missing | **advisory, continue** | designs legitimately fail gates and fail to fold; expected early in a sweep |

Verified against the real §2.5 numbers (d_base 2.486 / d_oxy 3.338 inherited,
dihedrals genuinely varying at 167/114/128/138/58/13/...): it aborts and names
exactly `['base', 'oxy']`, leaving the five terms that do vary alone. The
verdict is written to `outputs/runs/m4/checkpoint.json` either way.

---

## 8. Open questions, in value order

1. **[open]** Exact dihedral band edges and `RAMP_MODE`. Under `flat`, anything
   inside the window scores 1.0 and the ranking **reorders**. denovoKemp
   confirms the 180° centre but not our edges. One message to the organisers;
   cheapest high-value action available.
2. **[open]** Does the organisers' scoring fold with an MSA? Slide 7
   acknowledges MSA depth matters but does not say what they do.
3. **[open]** Burial anchor convention (8% weight).
4. **[open]** Does TM ≤ 0.7 apply to all published structures or only to
   published Kemp eliminases? Decides whether scaffold-based routes are legal.
5. **[open]** Whether to aggregate subscores rather than geometry (§1.3).

---

## 9. Things not to re-litigate

- Theozyme reference geometry — round-trips to ≤0.005 Å.
- The scoring implementation in `official.py` — geometric mean, fails closed.
  The circular-mean bug was in `combine_replicates`, now fixed, with the
  reasoning in a comment at the fix site.
- The novelty gate — Foldseek + MMseqs, working.
- `kemp/strain.py` — validated against crystal structures.
- `max_distance` pocket constraint — proven near-useless and *worse* with an
  MSA. Delete the knob, do not tune it.
- ProteinMPNN instead of LigandMPNN — the composition pathology is in the
  *invocation* (temperature / bias / `omit_AAs`); the teammate hits pLDDT
  0.80–0.90 on the identical protocol. A ligand-blind model for a catalytic
  pocket would misdiagnose it.

## 6.12 The disulfide question — measured, with a control arm

A teammate raised: *"they are overall quite unstable because they lack any
sort of Cysteine bond, so probably they won't express or fold like predicted."*

Checked rather than argued. **The premise is partly right and the conclusion
is wrong.**

Right: the reference enzyme HG3.17 (4BS0) *does* carry a disulfide —
Cys255–Cys261 at 2.03 Å, confirmed from the mmCIF `disulf1` record and
re-measured as an SG–SG distance. I had expected none; checking corrected me.

Wrong: that disulfide is not what makes it stable, and its absence is not what
would make ours unstable.

| evidence | designs | HG3.17 control |
|---|---|---|
| Rosetta relax, REU/residue | **−3.55 to −3.66** | −3.27 |
| fa_rep per residue | 0.63 | 0.75 |
| buried hydrophobics / length | 0.333–0.353 | 0.357 |
| exposed hydrophobic % | 7.6–9.1 | 6.5 |
| instability index (Guruprasad) | 31.7–32.8 | 26.7 |
| GRAVY | −0.36 to −0.39 | −0.10 |
| net charge | −17 / 150 aa | −5 / 318 aa |
| Boltz-2 pLDDT (5 seeds) | 0.977–0.979 | — |

All five designs are **better packed per residue than the enzyme that has
actually been expressed, purified and crystallised.** Core burial matches it.
They are more hydrophilic and far more charged, which favours soluble
expression, not against it. HG3's stability comes from its scaffold — a
*Thermoascus aurantiacus* (thermophile) IGP synthase barrel — not from a
surface cysteine pair 18–20 Å from the active site.

Two further points:

- Each design has **exactly one** free cysteine (Cys30). A disulfide needs two,
  so none is possible — and a lone Cys is the safe case, since scrambling also
  needs two. This is not luck: the generators ran with `OMIT_EVERYWHERE = "C"`.
- Disulfides do not form in the reducing *E. coli* cytoplasm, so a designed
  one would not be present in the expressed protein anyway without periplasmic
  export or an oxidising strain.

**If a disulfide is wanted next round**, 13 residue pairs in `L_E_shell_I23L`
already sit in disulfide-compatible geometry (CB–CB 3.5–4.5 Å, CA–CA 4.0–6.5 Å,
|i−j| ≥ 5), all clear of Asp102/Gln46/Tyr42/Tyr103. The ones worth trying are
the **long-range** pairs, which is where disulfides actually buy stability by
cutting unfolded-state entropy: `Ala14–Thr67` (Δ53), `Met98–Ala130` (Δ32),
`Leu6–Ala29` (Δ23). Local pairs such as `Ala13–Ala22` (Δ9) buy little.

Script: `scripts/stability_relax.py`. Note it strips the ligand — the question
is fold stability, and scoring a ligand without its params would silently drop
it regardless.

**Caveat, stated plainly:** none of this is a melting temperature. Rosetta
score is a packing proxy, and the supplied Risso 2020 sheet is explicit that Tm
is *"a gate, not a predictor"* of activity. The honest claim is "nothing here
indicates a folding or expression problem", not "these are stable".
