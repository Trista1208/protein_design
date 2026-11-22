# Lessons from the Kemp eliminase hackathon — 21–22 Nov 2026

Written after a full day of debugging. Everything here is **measured on this
machine** unless marked *(unverified)*. Where a claim was made and later
retracted, both the claim and the refutation are recorded, because the
retraction is usually the more useful half.

Read §1 and §6 first if you are short of time.

---

## 1. The scoring function — what it actually does

### 1.1 It is a weighted GEOMETRIC mean. One zero kills everything.

```
S = Π sᵢ^wᵢ
base 0.25 · ipde 0.20 · oxy 0.20 · orient_base 0.12 · orient_oxy 0.08 · bury 0.08 · cons 0.07
```

A single zero term collapses the total regardless of the other six. `d00000`
had ipde 0.918 and bury 1.0 and still totalled 0.000 because `oxy` was 0.
`official.py` emits a `zero_warning` saying exactly this.

**The arithmetic weighted sum in `official.py` is a DIAGNOSTIC only.** It gives
0.611 where the geometric mean gives 0.000. Do not quote it.

### 1.2 The windows are TAPERED, not thresholds

| term | target | window | headroom above target |
|---|---|---|---|
| `d_base` | **2.52** | 2.0–3.6 | 1.08 Å |
| `d_oxy` | **3.10** | 2.0–3.6 | **0.50 Å** |

Score decays linearly from 1.0 at target to 0.0 at the nearer window edge.
**Being "just inside" is worth almost nothing.** `d_oxy` 3.52 scores 0.159;
`d_base` 3.52 scores 0.074. Every 0.1 Å above the `d_oxy` target costs 0.2 of
the subscore.

Corollary: a design sitting at the window edge is not "nearly scoring". Early in
the day I called a design "0.2 Å from scoring" when it was ~0.6 Å and ~1.3 Å
away in the terms that mattered. Compute the subscore, never eyeball the
distance.

### 1.3 There are TWO orthogonal geometric requirements

`orient = min(angle_score, dihedral_score)` — **both must be good.**

- angle window 80–160°, target 120°
- dihedral must fall in an sp2 band: **130–230 or 310–50**
- therefore the **dead zones are 50–130 and 230–310**

Distance and orientation are independent. Fixing a distance never unzeros a
design whose dihedral is out of band. This was the single most expensive thing
to learn: hours went into tuning distance constraints while the dihedral — which
nothing in the pipeline controls or even reports — was what zeroed six of eight
designs.

### 1.4 `cons` is sd_base across replicates, and it is multiplicative

High replicate variance zeroes a design on its own, independent of geometry.
Reference: **HG3.17 sd_base = 0.068 Å**. Team designs measured 0.891 and 1.270 Å
→ `cons` = 0 for both, before geometry was even considered.

### 1.5 An ABSENT quantity must score 0, not "unmeasurable"

`official.score(strict=False)` maps a missing `d_oxy` to subscore `None` and
then **excludes that term**, flattering the Q50M knockout to 0.825 against the
positive's 0.834 — nonsense. `scripts/controls_full.py:135-157` has the
`donor_absent` override; use it. Production (`pipeline.py:262`) calls
`official.score(raw)` with `strict=True` and fails closed, which is correct.

### 1.6 Three different "base angle" conventions are in circulation

| convention | theozyme | 4BS0 | a bad design |
|---|---|---|---|
| CG–OD2···N3 (teammate's table) | 126.95 | 130.87 | 105.00 |
| **OD2···N3–N2 (what the score uses)** | 124.64 | 130.60 | **97.33** |
| OD2–N3–N1 (`config.REF_GEOMETRY`) | 158.35 | 165.47 | 131.52 |

They agree within 0.3° on good structures and diverge by ~8° on bad ones —
exactly where ranking decisions get made. The official one is defined by
`ORIENT_ATOMS` in `official.py`.

### 1.7 The realistic ceiling is ~0.78–0.87, not 1.0

HG3.17's own orientation terms score 0.649 and 0.589. Measured on crystal
coordinates it totals ~0.87; folded through Boltz with an MSA, 0.834. **A
perfect enzyme does not score 1.0 under this metric.** Know when to stop
optimising.

### 1.8 Open questions for the organisers

1. **Ramp shape** — `RAMP_MODE` is unresolved in our own code (taper vs flat).
   Under `flat`, anything inside the window scores 1.0, which **reorders the
   designs**.
2. **Angle target** — 120.0 in `TARGETS` is provenance-flagged ("5TS cst"). The
   theozyme measures 124.6, the crystal 130.6.
3. **Exact dihedral band edges.** Our best design misses by 3.4°. Worth
   confirming before redesigning around an inferred boundary.
4. **Do they fold with an MSA?** See §3.1 — this changes everything.
5. **Does the TM ≤ 0.7 novelty limit apply to all published structures, or to
   published Kemp eliminases?** Decides whether scaffold-based routes are legal.

---

## 2. Scientific analysis — how to not fool yourself

### 2.1 Run the control BEFORE believing any measurement

Three findings died on controls today. All three were plausible, internally
consistent, and wrong.

| claim | control that killed it |
|---|---|
| "All catalytic residues sit at disallowed φ" | 4BS0/7K4Z/6TU6 crystals fail the same check (Asp127 φ=+88, Gln50 φ=+71). My Ramachandran boundaries were too strict — they condemned the reference enzyme. |
| "The ideal catalytic pose clashes in every design" | 4BS0 and 7K4Z report 21 clashes at 0.39 Å, yet they bind 6NT fine. The method was invalid (see 2.2). |
| "CA–CA span predicts dihedral failure" | Across 27 structures the correlation **inverted**: 2.047 Å mean deviation in the "in band" group vs 1.540 Å in the "out" group. |

The control that *validated* a finding is equally important: the peptide-bond
junction result survived because Boltz folds gave mean 1.327 Å, sd 0.007, **zero
outliers in 1000 bonds**, against RFD3 junctions 0.15–2.94 Å off ideal. Clean
tenfold separation, so the finding stands.

**Rule: if a check would condemn the reference structure, the check is wrong.**

### 2.2 Two H-bond distances do not determine a ligand pose

There is free rotation about the axis through the contact atoms, plus a possible
ring flip. A rigid transplant of the theozyme ligand places the contact end
correctly and the far end anywhere — up to 3.93 Å off at the nitro group.

Consequences:
- Single-point scoring of a hand-placed pose **cannot work** for this problem.
  The remaining degrees of freedom must be *optimised* under restraints.
- This is the same degree of freedom the dihedral term measures. Two independent
  lines of evidence converging on one uncontrolled variable.

### 2.3 A dihedral on an unbound ligand is not an angle

It is noise. `foldtest none` had the ligand 13 Å away and scored "both dihedrals
in band". Restrict all geometric analysis to structures where both distances are
in window; otherwise you are computing statistics on garbage.

### 2.4 Measurement noise scales with fold confidence

| | pLDDT | d_base spread |
|---|---|---|
| HG3.17, 5 seeds | 0.983 | **sd 0.068 Å** |
| RFD3 **seq4**, same sequence twice | 0.41–0.47 | **2.56 → 13.64 Å** |

For a low-confidence fold, `s2_true = max(0, s2_between − s2_within/n)` floors
at **zero**: replication cannot help and best-of-N selection selects noise.
(Contrast HER2, ICC 0.504, where "spend on seeds not designs" is correct.)

**CORRECTION — do not over-generalise this.** I quoted "our designs fold at
pLDDT 0.41–0.47" throughout 21–22 Nov. That is **seq4 and s0 only**. Measured
across all eight LigandMPNN sequences on the same backbone:

| seq | s0 | s1 | s2 | s3 | **s4** | s5 | s6 | s7 |
|---|---|---|---|---|---|---|---|---|
| pLDDT | 0.349 | **0.834** | **0.838** | **0.849** | 0.411 | **0.859** | 0.605 | 0.552 |
| scRMSD Å | 9.47 | 2.60 | 1.25 | 2.76 | 10.51 | **0.97** | 8.31 | 5.57 |

**Four of eight fold at 0.83–0.86 with scRMSD 0.97–2.76 Å.** The set is bimodal,
not uniformly broken. I generalised from seq4 because the previous session had
selected it for the foldtest — on its *geometry* (best `d_base` 2.557), not its
foldability. **Never characterise a population from a member that was selected
on a different axis.**

### 2.4b Foldability and catalytic geometry are in TENSION here

Across those eight sequences, `r(pLDDT, d_base) = +0.539` (n=8). `d_base` is
better when lower, so **higher confidence goes with worse geometry**:

- best `d_base` (2.557 Å) belongs to the sequence with pLDDT **0.411**
- highest pLDDT (0.859) has `d_base` **7.304 Å**

This is the same shape as HER2's finding that sequence-plausibility scores were
*negatively* correlated with interface quality. It is also a direct instance of
§2.6: I proposed composition as a proxy for design quality, and on this backbone
it has the **wrong sign** for the thing that scores.

n=8, so this is a hypothesis, not a law. But it means "make it fold better" and
"make it catalyse" may be competing objectives on a diffused backbone, and
optimising one blindly can cost the other.

### 2.5 Identical numbers across different designs are inherited, not earned

RFD3 reported `d_base` 2.486 / `d_oxy` 3.338 / CA–CA 10.606 **identical to three
decimals** across designs whose motifs landed at 167/114, 128/138 and 58/13. It
was a rigid copy of the input. It carried zero information about the design.
Any metric that does not vary between designs is measuring your input.

### 2.6 Verify the sign of a cheap proxy before ranking by it

HER2 measured every sequence-plausibility score as **negatively** correlated with
interface quality (−0.427 raw, −0.602 disattenuated) — a language model reverts
an affinity-matured parent toward germline. They changed the score's *role* (a
floor, not a maximand) rather than retuning weights.

I proposed three cheap proxies today (fg-RMSD, composition, CA–CA span) and
verified the sign of none of them. The one I did test inverted.

### 2.7 Separate "the measurement is broken" from "the design is bad"

Most of this day was the former. See §3.

---

## 3. Tooling — what works, what lies, what is missing

### 3.1 Boltz in single-sequence mode cannot place a ligand

**The most important empirical result of the day.** Same sequence, same seed, no
ligand constraint, MSA the only variable:

| | d_base | d_oxy | score | pLDDT | iPDE |
|---|---|---|---|---|---|
| HG3.17, no MSA | 2.157 | **12.698** | **0.000** | 0.738 | 3.00 |
| HG3.17, with MSA | 2.760 | **2.952** | **0.775** | 0.983 | 0.53 |
| crystal reference | 2.532 | 3.102 | 0.959 | — | — |

**A crystallographically perfect Kemp eliminase scores 0.000 under
single-sequence folding.** The "oxyanion never lands" signature — base roughly
OK, oxyanion catastrophically far, geometric mean zeroed — is reproduced exactly
by the evaluator on a known-good enzyme. It was never purely a design failure.

Also: MSA folds reproduce the crystal CA–CA span to 0.02–0.03 Å; no-MSA folds
are off by 0.84–1.02 Å. Single-sequence mode gets the catalytic residue
*spacing* wrong, not just the ligand.

### 3.2 But de novo designs have no MSA

The ColabFold server returns **depth 1** for an RFD3 design — the query and
nothing else. So:
- `--use_msa_server` and `msa: empty` are the same thing for our designs.
- Controls (natural, MSA) and designs (de novo, no MSA) are on **incomparable
  protocols**, and nothing in `designs.csv` or the ledger records which.
- **Record the MSA setting alongside every score.** This is how the control and
  the designs silently became incomparable.

### 3.3 The pocket constraint is distance-only and near-useless

`max_distance` has no orientation term. With an MSA it makes things slightly
*worse* (0.775 → 0.756). Without one it cannot rescue the oxyanion
(12.70 → 11.45 Å, still ~8 Å outside). Every hour spent tuning it from 6 Å to
5 Å to 3.5 Å was optimising one of two orthogonal axes — which is why tightening
always traded one contact against the other.

### 3.4 RFD3's self-reported metrics are unreliable

`n_chainbreaks` reported **0** for a design containing a **2.68 Å** C–N bond,
and 1 where 2 cuts were measured. The pre-MPNN gate trusted that number.
**Measure the coordinates.**

### 3.5 RFD3 welds the motif in with non-physical peptide bonds

All 10 outputs: junction C–N **0.15–2.94 Å** off ideal, while the rest of the
chain is within **0.141 Å**. Perfect separation by location. Since a structure
predictor builds physically valid chains, it cannot reproduce that junction, so
the catalytic residues must move during folding.

### 3.6 `pipeline.py:410` throws away 7 of 8 sequences

`header, seq = seqs[0]` — LigandMPNN generates 8 and the pipeline uses the
first, with no selection. Those 8 spanned `d_base` **2.557–9.329 Å**; the one
used scored 4.427. (But see §2.4 — selecting among them is only meaningful once
ICC > 0.)

### 3.7 The composition pathology is NOT tool-specific

| | Ala | Gly | aromatic | Trp | entropy | pLDDT |
|---|---|---|---|---|---|---|
| our LigandMPNN seqs | **24–27%** | 13% | **1.5–2%** | **0** | 0.73–0.81 | 0.35–0.86 *(see §2.4)* |
| **RFD4 inverse folding** | **26.7%** *(its own reported metric)* | — | — | — | — | — |
| teammate designs | 7.6–19.5% | 5–8% | 4.7–7% | 0–2 | 0.84–0.91 | 0.80–0.90 |
| **rfd4_design_045** | **8.7%** | 2.0% | **7.3%** | 0 | 0.899 | **0.96** |
| HG3.17 | 11% | 8% | **9%** | **7** | 0.95 | 0.98 |

**Both inverse-folding tools produce ~25% alanine.** RFD4-Proteina reports
`amino_acid_composition.ala_fraction` **0.2674** over 10 sequences in its own
`inverse_folding_n10_ll2.0` run — the same pathology as LigandMPNN's 24–27%.
RFD4 tracks it as a named metric, which suggests its authors consider it worth
watching.

So this is **not** a LigandMPNN defect. It is what inverse folding does on
diffused backbones generally, and any fix transfers across generators.

**But the design that folds well does not have it.** `rfd4_design_045` sits at
8.7% Ala / 7.3% aromatic, close to HG3.17. Whether that came from different
sampling settings or from downstream selection is not established here.

Reading: composition looks weaker as a *single cause* and stronger as a
*selection criterion*. Filtering or biasing on it is cheap and
generator-independent.

Note this also contradicts a comment in `BACKBONE_GATE`: *"LigandMPNN REPLACES
the sequence wholesale, so the backbone's alanine fraction says little"* — the
problem survives inverse folding rather than being fixed by it.

### 3.7b The composition bias test — RUN, and the result is mixed

`scripts/mpnn_bias_test.py`. Same backbone, same fixed catalytic residues, same
seed 37 and temperature 0.1. **Only `--bias_AA` changed**
(`A:-2.0,G:-1.0,W:1.5,F:1.5,Y:1.5,L:0.5,I:0.5,M:0.3`).

**It did exactly what it was asked to, on composition:**

| | Ala | aromatic | Trp | entropy |
|---|---|---|---|---|
| baseline (8 seqs) | 22.8% | 1.6% | 0.0 | 0.807 |
| **biased (8 seqs)** | **8.8%** | **4.6%** | **2.0** | **0.897** |
| rfd4_design_045 | 8.7% | 7.3% | 0 | 0.899 |

Catalytic residues preserved in all eight. Gly barely moved (13%) — the −1.0
bias was too weak.

**Foldability: mean improved, best case did not.**

| | scRMSD mean | scRMSD best | pLDDT |
|---|---|---|---|
| baseline | 5.180 Å | **0.971 Å** | 0.349–0.859 (bimodal) |
| biased | **2.832 Å** | 2.373 Å | 0.554–0.830 |

The bias **compressed the distribution** — it removed the catastrophic folds
(9–10 Å scRMSD) and also removed the best one. If you want reliability it helps;
if you want a single champion it does not.

**Catalytic geometry: worse, in all four folded.**

`d_base` 18.195, 9.431, 7.385, 4.381 — every one outside the 2.0–3.6 window,
against a baseline that had one at 2.557. See §2.4b: on this backbone, folding
better and catalysing are pulling in opposite directions.

**Verdict:** composition bias is a real, one-parameter lever on *reliability*,
and it transfers across generators (RFD4 has the same pathology). It is **not**
a fix for catalytic geometry, and on this backbone it cost geometry. Use it to
reduce variance, not to chase a champion — and re-check the sign on your own
backbone before trusting it.

Raw results: `outputs/mpnn_bias/results.json`, `logs/mpnn_bias.log`.

**Trap:** `--bias_AA` takes `A:-2.0,W:1.5`, **not** JSON. Passing JSON gives
`ValueError: could not convert string to float: ' 0.3}'`.

### 3.8 Predictors disagree about rotamers, and neither is reliable

On `rfd4_design_045`, OpenFold3 chose an oxyanion dihedral of **86.1°** and
Boltz-2 chose **~237°** (mean 237.2, sd 3.06 over 5 seeds). Both out of band, in
*different* dead zones. A coordinate repack will not survive re-folding.

### 3.9 Tool notes

| tool | verdict |
|---|---|
| **OpenMM + openmmforcefields + OpenFF + RDKit + PDBFixer + ParmEd** | all installed. The right stack for restrained minimisation. |
| **gnina** | not installed. No distance restraints; CNN optimises binding affinity, not catalytic geometry. Wrong tool. |
| **PLIP** | not installed (needs OpenBabel). Useful for *diagnosis* — what holds the ligand, π-stacking — not scoring. |
| **Rosetta** | built, `match.linuxgccrelease` available. No scaffold library though. |
| **Foldseek + MMseqs** | installed and already wired into the novelty gate. Works — rejected `d00002` at TM 0.705. |
| **RFD4** | env exists but **no binary and no checkpoint**. The team runs it elsewhere. |
| **OpenFold3** | emits `of3_min_ipae`, **not** Boltz `complex_ipde`. RFD4 designs folded with OF3 cannot be scored on the 20%-weight ipde term without a Boltz refold. |

### 3.10 Organiser-recommended tools

> Enzyme design: **RFD4-Proteina**, RFDiffusion3, EnzyGen2, Disco
> Optimization: ProteinMPNN and **RFOptimization**
> Binders: RFD4-Proteina, BoltzGen, BindCraft2

We were on RFdiffusion3; **RFD4-Proteina is listed first**. The handover said
"do not put RFD4 on the critical path" — that call looks wrong: the team's best
designs are RFD4, and RFD4 reports 12.1% vs 5.7% cluster pass rate on enzyme
scaffolding with an explicit `enzyme_design` task. **RFOptimization** is the
organisers' own pointer at the optimisation loop.

### 3.11 Small traps that cost time

- `official.measure_raw` reads **PDB only**. Convert CIF first
  (`pipeline.cif_to_pdb`). Cost two crashes.
- **Dihedral sign conventions bite.** A `b0 = a−b` instead of `b−a` gives a 180°
  offset that makes trans peptides read as cis. Validate against known values.
- mmCIF `_atom_site` column order is **not** the PDB's. Parse by header name or
  you will read `B_iso` as a coordinate.
- Auto-locating the catalytic pair by "nearest Asp to N3" **invents** a pair
  where none exists (it found a spurious Gln in the Q50M knockout) and picks the
  wrong residue when the true indices are known. Use the index map.

---

## 4. Method — what the sibling projects do that this one doesn't

From HER2 (`BreastCancerHackathon2026`), Catalyst, and FGFR2.

1. **Neither HER2 nor Catalyst lets a model write a sequence.** Both perturb a
   known-good parent: substitution-only, fixed length, 5–8 edits, a
   geometry-derived frozen core (HER2 freezes anything within **3.0 Å** of the
   target, derived at runtime), and a per-position **top-5 whitelist** from a
   one-time masked scan. The composition pathology is structurally impossible
   there, not filtered after the fact.
2. **Pre-compute a per-position score table once; candidate scoring is a
   lookup.** HER2 scored **1,000,000 candidates for 515 s of CPU**. This
   decouples proposal count from fold count — the real answer to "scale".
3. **Percentile-rank scalarisation, never a product**, for internal ranking. A
   missing term contributes 0 to that axis only. Gives a usable ordering even
   when every official score is 0.000. (`robust_rank_scores`,
   `Boltz-2/run_boltz2_pipeline.py:1145-1177`.)
4. **Degeneracy audit** — FAIL the run when a score column is all-NaN, >50%
   missing, has ≤1 unique value, or std < 1e-8. ~20 lines
   (`ESM/audit_esm_pipeline_run.py:428-452`). **Would have fired on day 1 here.**
5. **Hard gates separate from soft ranking**, with named rejection reasons
   counted into a histogram. You get *why* the campaign died, not just a count.
6. **Per-objective band tables** (strong/acceptable/borderline/red per metric,
   decoupled from the rank) — the "which term is the zero?" instrument. Our two
   dihedral bands map onto this directly.
7. **Statistics** (HER2 `stage6_select_audit.py:285-406`): ICC variance
   decomposition → seeds-for-reliability → Welch vs baseline → BH FDR at
   q<0.05. Lift wholesale. **Fix their bug on the way in:** compare Δ against
   the SE of the difference, not the single-seed SD (their two criteria
   disagree because of this).
8. **Budget the expensive stage.** 500 designs = **4h41m on this exact GPU**
   (RTX 3090), folding being 58% of it. Scale was never our blocker — the
   campaign stopped at 5 by choice, at ~49 s/design.
9. **Operational**: resume by per-design artifact completeness; non-fatal job
   failures with per-design status; ranked snapshot of an in-flight campaign;
   GPU memory preflight; stop-loss artifact written at the earliest stage;
   smoke-test abort.
10. **Counter-screen.** FGFR2 counter-screens binders against FGFR1. Our unused
    analogue sits in `config.DECOY_LIGANDS`: score against **5NT** (wrong nitro
    position) as well as 6NT. A generic aromatic pocket holds both; a real TSA
    site discriminates — **and it still discriminates when everything scores
    0.000**.
11. **Neither project has an optimisation loop.** Both are one-pass cost-tiered
    cascades. Don't go looking for a GA or BO to copy.
12. **Neither has angular constraints.** Distance only, everywhere. The dihedral
    problem is ours alone.

---

## 5. What is trustworthy

Don't re-litigate these:

- **Theozyme reference geometry** — round-trips to ≤0.005 Å.
- **The scoring implementation** — geometric mean, fails closed. My *scripts*
  had the bugs, not `official.py`.
- **The novelty gate** — Foldseek + MMseqs, working.
- **`kemp/strain.py`** — chain integrity, dihedrals, catalytic distances, all
  validated against crystal structures. 77 tests pass.
- **The control pair, folded with an MSA** — 0.834 vs 0.000, ranks correctly.
  The old "STOP" verdict was a protocol artefact of single-sequence folding.

---

## 6. State of play, and the one open problem

**Best design: `rfd4_design_045`** (150 aa, RFD4 + OpenFold3).

| term | value |
|---|---|
| base | 0.687 |
| oxy | **0.994** |
| ipde | **0.914** (from the Boltz refold) |
| orient_base | 0.364 |
| **orient_oxy** | **0.000** ← the only zero |
| bury | 1.000 |
| sd_base | **0.0759** excluding one outlier seed (HG3.17: 0.068) |

Everything is solved except the **oxyanion donor orientation**.

- Boltz refold, 5 seeds: dihedral **233.4–239.8°**, mean 237.2, **sd 3.06°**.
  Band edge is 230. **It misses by 3.4°** — systematic and reproducible, not
  noise.
- A rotamer scan found **360 clash-free Gln rotamers** that satisfy both the
  distance window and the dihedral band; the best gives `orient_oxy` 0.745 and
  a 5-term geometric mean of **0.721**.
- **But 0/5 seeds sampled one.** So a coordinate repack will not survive the
  organisers re-folding the sequence.
- Minimal nudges (Δchi ≤ 30°) do reach the band but score ~0.17 — they land at
  the *edge* of both window and band, where the taper makes them worthless.

**The fix is a design change near Gln46** that disfavours the ~237° rotamer and
favours one near 180°. The target is specified, the viable rotamers are
enumerated, and each candidate mutation is testable with a 4-minute refold.

---

## 7. If starting over

1. **Get a trustworthy measurement first.** Nothing else means anything until
   folds are confident. A design at pLDDT 0.41 is unmeasurable, not bad.
2. **Fold a known enzyme through your exact protocol** before trusting any
   number. That one test would have found the MSA problem on day 1.
3. **Report every term, every design, every replicate** — especially the
   dihedrals. Six of eight designs died on a term no table showed.
4. **Add the degeneracy alarm.** All-identical score columns must fail loudly.
5. **Rank by percentile sums internally**, keep the geometric mean for
   reporting.
6. **Don't generate whole sequences.** Perturb a foldable parent with a frozen
   catalytic core.
7. **Ask the organisers §1.8 early.** Several of those answers change the
   strategy, and all are cheap to ask.
