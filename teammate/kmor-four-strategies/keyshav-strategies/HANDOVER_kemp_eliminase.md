# Kemp Eliminase Hackathon — Verified Handover for Claude Code

**Prepared:** 21 November 2026
**Supersedes:** `Claude_Kemp_Eliminase_Implementation_Brief.md` (the ChatGPT brief). That document is
largely sound and its compliance thinking is good. Section 1 below lists the places where it is
wrong or where it missed something, because several of those errors are silent and expensive.
**Hardware assumed:** 1 × RTX 3090 (24 GB), optionally one more cluster GPU.
**Inputs:** `DBSSE_Hackathon_Perlin.pdf`, `Theozyme_withORI.pdb`, `Theozyme_noORI.pdb`.

> **Reading rule for the implementer.** Everything in §1–§3 was verified against primary sources
> on 21 Nov 2026 and the verification is cited. Everything in §4 onward is a plan, not a fact.
> Where this document and the slides disagree, the slides win. Where this document and a pinned
> repository's actual CLI disagree, the repository wins — read the installed version's schema
> before writing any config.

---

## 1. Corrections to the previous brief

### 1.1 The theozyme provenance is solved: it is 6TU6, with Gln50 modelled back in

The previous brief flagged a "critical provenance mismatch" and left it unresolved. It is
resolvable in one line. The `CRYST1` record of both supplied files is:

```
CRYST1   98.370  108.690  108.720  88.53  90.15  89.82 P 1          16
```

PDB **6TU6** has unit cell 98.370, 108.690, 108.720 / 88.53, 90.16, 89.82, space group **P 1**,
**16 chains** (A–P). The supplied atoms are chain **G**. That is a match on every cell length to
three decimals, on the space group, and on the chain count. The files are cropped from 6TU6.

6TU6 is *Kemp Eliminase HG3.17 mutant Q50M, E47N, N300D complexed with 6-nitrobenzotriazole*,
1.05 Å, Bloch & Hilvert, from Kries, Bloch, Bunzel, Pinkas & Hilvert, *ACS Catalysis* 2020,
10(8):4460–4464, "Contribution of Oxyanion Stabilization to Kemp Eliminase Efficiency".

The brief was right that 6TU6 has **Met**, not Gln, at position 50 — and that is exactly why the
file looks the way it does. Evidence from the file itself:

| Residue | Occupancy | B-factor | ANISOU | Hydrogens |
|---|---|---|---|---|
| ASP 127 | 0.98 | 8.1–13.2 | yes (real refined) | yes |
| 6NT 401 | 0.90 | 10.7–13.7 | yes (real refined) | yes |
| GLN 50  | 1.00 | **0.00 for all 17 atoms** | **none** | yes, idealised |

Asp127 and the ligand are genuine crystallographic coordinates from 6TU6. **Gln50 was modelled
in** — the Q50M mutation was computationally reverted. That is the whole story, and it has
consequences:

1. The **3.33 Å Gln NE2···N1 distance is a model's opinion, not an observation.** Do not treat it
   as ground truth for the Geometry metric.
2. The **6NT pose comes from a structure that had no Gln50.** 6TU6 is the oxyanion-stabiliser
   knockout. The ligand may sit slightly differently in a Gln-containing structure.
3. **Action (M0, 30 min, high value):** cross-check the supplied theozyme against **4BS0**
   (HG3.17 E47N/N300D + 6NT, 1.09 Å, Blomberg et al. *Nature* 2013) and **7K4Z** (HG3.17 + 6NT,
   1.08 Å, Kern lab). Both have a real Gln50. If the organisers' modelled Gln50 disagrees with the
   crystallographic one, you want to know before you scaffold 500 backbones around it.

### 1.2 The two supplied files are in different coordinate frames — ORI is not transferable

Verified by Kabsch superposition: the 45 shared atoms superpose at **RMSD 0.0005 Å**. Identical
motif, different frame (6TU6 chain G and a symmetry copy).

The `ORI` pseudo-atom exists only in `Theozyme_withORI.pdb`, at `[25.568, 0.048, -23.462]`, which
is **0.99 Å from ligand H7** and 3.8 Å from the ligand centroid.

**The trap:** if you use `Theozyme_noORI.pdb` as the RFD3 `input` but paste the ORI coordinates
from the withORI file into `ori_token`, the origin lands **70.5 Å from the ligand**. RFD3 will
not error. It will quietly try to build a protein whose centre of mass is 70 Å from your active
site, and every design in that batch is wasted.

Mapped into the noORI frame, ORI would be `[-33.913, 30.683, -33.510]`.

**Rule: pick ONE frame and never mix.** Use `Theozyme_withORI.pdb`'s frame, since the ORI
coordinate is meaningful there. Assert in code that `||ori_token − ligand_centroid|| < 12 Å`.

### 1.3 ORI is a JSON field, not an atom

The brief implies building an ORI pseudo-atom into the PDB. The RFD3 enzyme tutorial supplies it
as a JSON field, `"ori_token": [x, y, z]`. The `HETATM ... ORI` line in the supplied file is a
**visualisation** of the intended token, not an input format. You read its coordinates into the
JSON and you **strip the atom** before anything touches RDKit, Boltz-2, or the scoring geometry.
Note the supplied ORI atom also has a 3-character element field, which is invalid PDB and will
break RDKit and OpenBabel outright.

### 1.4 Verified RFD3 enzyme-design schema (the brief guessed; these are the real field names)

From the RFD3 enzyme design tutorial (`rosettacommons.github.io/foundry/models/rfd3/tutorials/enzyme_design_tutorial.html`),
settings live in a named calculation block:

| Field | Meaning | Tutorial example |
|---|---|---|
| `input` | path to the theozyme structure | `"path/to/1euv_lig.pdb"` |
| `ligand` | ligand identifier, `chain:resname` style | `"l:g"` (colon prevents CCD matching) |
| `unindex` | motif residues with unspecified sequence position | `"A514,A531,A574,A579-581"` |
| `length` | designed protein length range | `"100-200"` |
| `ori_token` | centre of mass of the designed portion | `[0,1,0]` |
| `select_fixed_atoms` | residue → atoms held fixed | `"A580":"SG,CB,CA,N,C,O"` |
| `select_buried` / `select_exposed` | RASA conditioning on ligand atoms | — |
| `select_unfixed_sequence` | residues whose identity may change | `"A579,A581"` |

The ligand is held in place automatically — you do not fix it yourself. Do **not** fix side-chain
atoms of residues listed in `select_unfixed_sequence`. Batch size is controlled by
`diffusion_batch_size` on the `rfd3 design` command; default produces 8 designs numbered 0–7, and
**all designs in a batch share one length**, so use batch size 1 if you want length variety.

There is no `contig` or `motif` field in the enzyme tutorial — the brief's mental model was the
RFdiffusion1 contig string, which does not apply here. (A `contig` does appear in the RFD3
protein–protein interface tutorial, so don't copy across tutorials either.)

The docs give no VRAM figure; they say only that RFD3 "runs best on GPUs" and suggest an
interactive GPU node. **Measure it yourself in M1.**

### 1.5 Attribution corrections

- **HG3 is not a Baker-lab enzyme.** HG3 comes from Privett, Kiss, Lee, Blomberg, Chica, Thomas,
  Hilvert, Houk & Mayo, *PNAS* 2012, 109(10):3790–3795 — Mayo/Houk/Hilvert. HG3.17 is Hilvert-lab
  directed evolution (Blomberg et al., *Nature* 2013). **KE07/KE70** are the Baker-lab designs
  (Röthlisberger et al., *Nature* 2008). The novelty filter must cover all of them regardless, so
  this changes nothing operationally — but do not write "Baker-lab HG3" in a submission.
- **3NYD** is HG-2 with **5-nitrobenzotriazole**, not 6NT. Different ligand. Check before use.
- **Listov numbers:** the *Nature* abstract gives efficiency **12,700 M⁻¹s⁻¹** and rate **2.8 s⁻¹**
  (brief said 2.85). The improved variant reaches **>10⁵ M⁻¹s⁻¹ and 30 s⁻¹**, described in the
  abstract as the result of *"designing a residue considered essential in all previous Kemp
  eliminase designs"* — see §3.3, this is mechanistically important and the brief missed it.

### 1.6 The hardware claim that matters most

The brief's compute table is vague where it needs to be precise. The single most important fact:

**Boltz-2 runs fine on a 24 GB RTX 3090 via the open-source `boltz` package.** The 48 GB
requirement the brief implies is for the **NVIDIA NIM container**, not the open-source release.
ChimeraX's benchmarks report Boltz 2.1.1 on Linux with a 24 GB NVIDIA GPU (3090/4090) predicting
~1300 residues in ~4 minutes. Your designs are ≤260 aa plus one small ligand — comfortably inside
that. **Use the open-source `boltz` package (github.com/jwohlwend/boltz), not the NIM.**

---

## 2. What the challenge actually requires (re-read from the slides)

| Rule | Value | Where |
|---|---|---|
| Reaction | Kemp elimination of 5-nitrobenzisoxazole | p.4 |
| Scoring ligand | 6-nitrobenzotriazole (TSA) | p.5 |
| Max length | 260 AA (motif residues included) | p.7 |
| Sequence novelty | >70% identity to published enzymes disqualifies | p.7 |
| Structure novelty | TM-score >0.7 to published structures disqualifies | p.7 |
| Replicates | 5 per sequence | p.7 |
| Metrics | iPDE, Geometry, Catalytic base burial (CA within 12 Å), Consistency | p.6 |
| Lab slots | ≥5 designs per team, independent of score | p.7 |

Reference geometry measured from the supplied theozyme (verified, round-trips to deviation 0.0):

| Quantity | Value |
|---|---|
| Asp127 OD2 ··· 6NT N3 | **2.48 Å** |
| Asp127 OD2 ··· H3 | 1.53 Å |
| Gln50 NE2 ··· 6NT N1 | **3.33 Å** *(modelled, see §1.1)* |
| Asp127 CA ··· Gln50 CA | **10.64 Å** |
| OD2–N3–N1 angle (heavy atoms) | **158.4°** |
| OD1 ··· N3 (non-engaging oxygen) | 4.56 Å |

The 10.64 Å CA–CA span is the real scaffolding constraint.

---

## 3. Strategy assessment

### 3.1 Verdict on the brief's four strategies

| Brief's strategy | Verdict | Why |
|---|---|---|
| 1. Rosetta + FuncLib on a known enzyme | **Cut as a design path.** Keep sequences as controls. | Output is disqualified by construction. PROSS and FuncLib are queued Weizmann web services — turnaround is unpredictable and unsuited to a hackathon clock. |
| 2. AI.zymes full evolutionary loop | **Cut the loop. Steal FieldTools.** | Needs a Rosetta licence and substantial setup; MIT-NC licence; and it needs a compliant scaffold to start from, which you don't have until Pipeline C runs. See §3.4 for why FieldTools alone is worth keeping. |
| 3. RFdiffusion3 motif scaffolding | **This is the plan.** 80% of effort. | Only path that is novelty-compliant, fits a 3090, and matches the organisers' own stated pipeline (slide 3). |
| 4. denovoKemp / Rosetta Matcher | **Cut as a pipeline. Steal two insights.** | Listov et al. is a multi-month CPU campaign with Rosetta licensing. Not reproducible on this clock. The insights are in §3.3. |

The brief's suggested effort split (C 45% / D 25% / B 20% / A 10%) spreads one GPU across four
workstreams and will finish none of them. **Recommended split: C 70% / controls+scoring 20% /
timeboxed lottery tickets 10%.**

### 3.2 The control pair the organisers accidentally handed you

This is the highest-value idea in this document.

The theozyme was cropped from **6TU6 = HG3.17 Q50M**, the oxyanion-stabiliser knockout from a
paper *specifically about how much oxyanion stabilisation contributes to Kemp eliminase
efficiency*. So you have, free and published:

- **Positive control:** HG3.17 (Gln50 present) — the most active Kemp eliminase of its lineage.
- **Negative control:** HG3.17 Q50M — same fold, same base, oxyanion stabiliser deleted,
  experimentally known to be worse.

Run both through your full Boltz-2 ×5 scoring pipeline before you generate a single design.

**If your pipeline cannot rank HG3.17 above HG3.17-Q50M, your Geometry and Consistency metrics are
not detecting oxyanion stabilisation, and every number you produce afterwards is noise.** This is a
one-hour test that validates or invalidates the entire weekend. Nobody else will run it.

It also calibrates: it gives you the seed-to-seed variance of a *real, crystallographically
verified* enzyme, which is your noise floor. Any design scoring within that spread of another is
indistinguishable from chance, and ranking inside it is superstition.

### 3.3 Two insights worth stealing from Listov et al. 2025

1. **TIM barrels, and you are right at the size ceiling.** Listov's workflow designs in TIM-barrel
   folds from natural backbone fragments. Their deposited designs (**9HVB, 9HVG, 9HVH**, DOI
   10.1038/s41586-025-09136-2) are **263 residues** — *above* the hackathon's 260 AA cap. The
   state of the art for this exact reaction sits just over your limit. Set RFD3 `length` to
   roughly `"190-250"` to leave headroom, and count motif residues against the cap.
2. **The "essential" residue may not be essential.** The abstract reports that *designing a
   residue considered essential in all previous Kemp eliminase designs* pushed efficiency to
   >10⁵ M⁻¹s⁻¹ and 30 s⁻¹. Read the paper before committing all five lab slots to the canonical
   Asp+Gln arrangement. At minimum, make one of your five designs a motif variant.

### 3.4 Know who is marking your work

**Lucas Merlicek — the presenter on slide 1, CTO and co-founder of Perlin Bio — is the first
author of AI.zymes** (Merlicek, Neumann, Lear, Degiorgi, de Waal, Cotet, Mulholland & Bunzel,
*Angew. Chem. Int. Ed.* 2025, 64:e202507031, DOI 10.1002/anie.202507031). AI.zymes was benchmarked
on exactly this chemistry — promiscuous Kemp eliminase activity in ketosteroid isomerase — giving
a 7.7-fold activity increase after zero-shot experimental validation of only 7 variants.

And **H. Adrian Bunzel**, its corresponding author, is a co-author of the ACS Catalysis paper
behind **6TU6**, the structure your theozyme was cut from.

Two practical consequences, not flattery:

- The central thesis of AI.zymes is that **electrostatic preorganisation** — the catalytic
  electric field projected on the breaking bond — captures catalytic quality that structural
  confidence metrics miss. Running **FieldTools** (github.com/bunzela/FieldTools, CPU, cheap) on
  your five finalists is a mechanistically meaningful tiebreaker that iPDE cannot give you, and it
  speaks the organisers' language. Keep it. Report it as a proxy, explicitly not as ΔG‡.
- Do not submit anything KSI-derived. They will recognise it instantly, and it fails the novelty
  rules anyway.

### 3.5 What actually wins these competitions

From the one well-documented post-mortem of a comparable event — *Crowdsourced Protein Design:
Lessons From the Adaptyv EGFR Binder Competition* (bioRxiv 2025.04.17.648362):

- The **overall winner (Cradle Bio) optimised an existing framework**, reaching 1.21 nM, ~8×
  better than the commercial antibody. **That route is closed to you by the novelty rules** — this
  hackathon has deliberately removed the strategy that wins when it is allowed.
- The best *de novo* binder came from **BindCraft** at 82 nM, and it had no close structural match
  in either the PDB or the AFDB. The lesson is not "use BindCraft" (wrong problem) but that **a
  tight integrated design-and-filter loop beat large undirected generation sweeps.**
- Of 400 proteins tested, **378 expressed (95%)** but only **53 bound (14%)**. Hit rates in this
  field are low. **Diversity across your five slots is worth more than depth on one pocket.**
- Ranking evolved from AlphaFold-Multimer **iPAE** in round 1 to adding **ESM2 pseudo-log-likelihood**
  in round 2, specifically to favour "natural-looking" sequences more likely to express.

**That last point is the cheapest edge available to you.** The official scoring rewards nothing
about expressibility — but your five lab-validated designs are worthless if they don't express.
Add an ESM2 pseudo-log-likelihood filter and basic solubility heuristics as a *final* tiebreak
among compliant, well-scoring designs. It costs almost nothing and it is the difference between
five data points and five empty wells.

### 3.6 Timeboxed lottery tickets (10% of effort, hard stop)

- **AlphaProtein Novo has a built-in `kemp_eliminase` evaluation suite.** Its design manifest
  takes `evaluation` with an enzyme suite chosen from `"kemp_eliminase"`, `"serine_esterase"`,
  `"dehp_esterase"`, `"carbene_transfer"`, `"nitrene_transfer"`
  (github.com/google-deepmind/alphaprotein-novo, released Oct 2026, code Apache-2.0,
  non-commercial generator weights). DeepMind built a scorer *and* generator for your exact
  reaction. **Worth one person-hour to check whether it installs and runs on a 3090.** If it
  works, it is a free second generator and an independent scorer. If it doesn't install in an
  hour, drop it and never look back.
- **RFD4-Proteina** (NVIDIA + IPD) reports a 12.1% cluster pass rate on enzyme scaffolding versus
  5.7% for AlphaProtein-Novo co-design, with 5.9× throughput, and has an explicit `enzyme_design`
  task for scaffolding an unindexed catalytic site. **But it was released on build.nvidia.com and
  NGC on 16 Nov 2026 — five days ago** — and is distributed as a NIM container likely wanting far
  more than 24 GB. Do not put it on the critical path. If someone has a big cluster GPU and a
  spare hour, try the hosted endpoint.

---

## 4. Recommended pipeline

```
M0  AUDIT & CONTROLS            (CPU, ~2h)
    verify theozyme vs 4BS0/7K4Z · fix frames · strip ORI · build ligand chemistry
    HG3.17 and HG3.17-Q50M through full scoring  ->  CALIBRATION + NOISE FLOOR
                         |
M1  RFD3 SMOKE TEST            (GPU, ~2h)
    one config, batch 1, measure VRAM and wall-clock  ->  THROUGHPUT NUMBER
                         |
M2  LIGANDMPNN + BOLTZ SMOKE   (GPU, ~2h)
    1-3 sequences, catalytic residues FIXED, one Boltz run, extract raw iPDE
                         |
M3  SCORING + HARD GATES       (CPU, ~3h)
    4 metrics as separate functions · length/identity/TM gates · unit tests
                         |
M4  SCALE                      (GPU, overnight)
    ORI sweep x motif variants x lengths -> backbones -> 8 seqs each
    prefilter hard, Boltz x5 ONLY on finalists
                         |
M5  SELECT                     (CPU, ~2h)
    cluster for diversity · ESM2 pLL expressibility tiebreak · FieldTools on final 5
```

Each milestone gates the next. Do not start M4 until M0's control pair ranks correctly.

### 4.1 Budget arithmetic — do this before generating anything

With one 3090, Boltz-2 at ~1–2 min for a 260 aa + ligand complex, and 5 seeds per finalist, a
finalist costs ~5–10 GPU-minutes. One overnight run (12 h) is roughly 70–140 finalists *if Boltz
is the only GPU consumer*, which it isn't — RFD3 and LigandMPNN need the same card.

Realistic plan: **Boltz ×1 as a cheap prefilter on a few hundred designs, Boltz ×5 on 20–40
finalists only.** Measure the real numbers in M1/M2 and re-plan; do not trust this paragraph.

### 4.2 Motif variants to scaffold (diversity across the five lab slots)

1. **Canonical** — Asp base + Gln donor, as supplied.
2. **Crystallographic Gln** — same, but with the real Gln50 rotamer from 4BS0/7K4Z if it differs.
3. **Backbone-NH donor** — Asp base, oxyanion stabilised by a backbone amide instead of a side chain.
4. **Alternative donor** — Ser/Thr/Tyr hydroxyl in place of Gln.
5. **Base variant** — Glu instead of Asp (longer reach, different rotamer space).

Scaffold at least three. Five correlated designs from one pocket is the single most likely way to
waste all five lab slots.

### 4.3 Non-negotiable implementation rules

1. **Fix the catalytic residues in LigandMPNN.** It will otherwise mutate your base away. Assert
   their identity in the output sequence and fail loudly.
2. **One coordinate frame.** Assert `||ori_token − ligand_centroid|| < 12 Å` before every run.
3. **No ORI atom downstream.** Strip it before RDKit, Boltz-2, and all geometry.
4. **Altloc handling.** The supplied files mix altlocs (Asp `A`, ligand `B`, Gln blank). Naive
   parsers that keep only altloc `A` or blank will silently drop the ligand. Normalise to
   single-conformer first and assert the ligand survived.
5. **6NT chemistry.** 6-nitrobenzotriazole, C6H4N4O2, CCD `6NT`. Validate bond orders and tautomer
   against the CCD. A test must fail if 6NT is confused with H5J (5-nitro-1,2-benzoxazole, the
   real substrate from 7EPO) or with 5-nitrobenzotriazole (the 3NYD ligand).
6. **Identical MSA settings across all Boltz runs.** Slide 7 warns MSA depth affects de novo
   folding confidence. If settings drift between designs, confidence numbers are not comparable.
   De novo designs have no meaningful MSA — pick single-sequence mode and keep it.
7. **Seeds recorded, never cherry-picked.** Rank on median plus dispersion or worst-case, never on
   the best-looking seed. Consistency is an official metric; gaming it defeats its purpose.
8. **Novelty screened against a real corpus**, not just 4BS0/6TU6. Include KE07/KE70, HG-series,
   Des27 family (9HVB/9HVG/9HVH), KSI, and 1A53-2. Use Foldseek or US-align for TM and MMseqs2 for
   identity. Record both TM normalisation directions.
9. **Fail closed.** Missing score, missing reference, unresolved atom mapping — reject the
   candidate, do not default it through.

---

## 5. Verified source list

**Challenge chemistry and structures**
- 6TU6 — HG3.17 Q50M/E47N/N300D + 6NT, 1.05 Å, P1, 16 chains. *Source of the supplied theozyme.*
  https://www.rcsb.org/structure/6TU6
- Kries H., Bloch J.S., Bunzel H.A., Pinkas D.M., Hilvert D. "Contribution of Oxyanion Stabilization
  to Kemp Eliminase Efficiency." *ACS Catal.* 2020, 10(8):4460–4464. DOI 10.1021/acscatal.0c00575
- 4BS0 — HG3.17 E47N/N300D + 6NT, 1.09 Å. Blomberg et al. "Precision is essential for efficient
  catalysis in an evolved Kemp eliminase." *Nature* 2013. https://www.rcsb.org/structure/4BS0
- 7K4Z / 7K4Q — HG3.17 and HG3 + 6NT, 1.08 Å, Kern lab (Brandeis/HHMI). From Otten et al. "How
  directed evolution reshapes energy landscapes in enzymes to boost catalysis." *Science* 2020,
  370(6523):1442–1446. Confirms Asp127 as base and a *newly acquired* H-bond donor in HG3.17.
- 7EPO — *M. smegmatis* KSI + H5J (5-nitro-1,2-benzoxazole = the real substrate), 2.35 Å, Glu127
  proposed as base. Apo: 7EPN. https://www.rcsb.org/structure/7EPO
- 3NYD — HG-2 + **5**-nitrobenzotriazole (different ligand). Privett et al. *PNAS* 2012.
- HG3 origin: Privett, Kiss, Lee, Blomberg, Chica, Thomas, Hilvert, Houk, Mayo. "Iterative approach
  to computational enzyme design." *PNAS* 2012, 109(10):3790–3795. DOI 10.1073/pnas.1118082108

**Design methodology**
- Listov D. et al. "Complete computational design of high-efficiency Kemp elimination enzymes."
  *Nature* 2025, 643(8074):1421–1427. DOI 10.1038/s41586-025-09136-2. PMID 40533551.
  Structures 9HVB/9HVG/9HVH (263 aa). Code: https://github.com/Fleishman-Lab/denovoKemp
- Merlicek L.P. et al. "AI.zymes: A Modular Platform for Evolutionary Enzyme Design."
  *Angew. Chem. Int. Ed.* 2025, 64:e202507031. DOI 10.1002/anie.202507031.
  https://github.com/bunzela/AIzymes · https://github.com/bunzela/FieldTools
- Butcher J. et al. "De novo Design of All-atom Biomolecular Interactions with RFdiffusion3."
  bioRxiv 2025.09.18.676967. SI reports 100 designs per motif, 8 LigandMPNN sequences per backbone.
- Ahern W. et al. "Atom level enzyme active site scaffolding using RFdiffusion2." *Nature Methods*
  2025. bioRxiv 2025.04.09.648075.
- RFD3 docs and enzyme tutorials: https://rosettacommons.github.io/foundry/models/rfd3/index.html
- Baker-lab enzyme-design worked examples (RFD2/RFD3 tutorials, demo pipeline, design data):
  https://github.com/baker-laboratory/Metallohydrolase_Enzyme_Design
- AlphaProtein Novo (has a `kemp_eliminase` evaluation suite):
  https://github.com/google-deepmind/alphaprotein-novo
- RFD4-Proteina model card (released 16 Nov 2026): https://build.nvidia.com/ipd/rfd4-proteina/modelcard
- Boltz-2: https://github.com/jwohlwend/boltz · bioRxiv 10.1101/2025.06.14.659707
- LigandMPNN: https://github.com/dauparas/LigandMPNN

**Competition evidence**
- Chen A. et al. "Crowdsourced Protein Design: Lessons From the Adaptyv EGFR Binder Competition."
  bioRxiv 2025.04.17.648362
- Cradle Bio write-up of the winning entry:
  https://www.cradle.bio/blog/adaptyv-protein-design-competition

**Hardware**
- Boltz on 24 GB consumer GPUs (~1300 residues, ~4 min):
  https://www.rbvi.ucsf.edu/chimerax/data/boltz-apr2025/boltz_help.html
- The 48 GB minimum applies to the **NVIDIA NIM container only**:
  https://docs.nvidia.com/nim/bionemo/boltz2/support-matrix.md

---

## 6. Open questions to resolve before scaling

1. Does the organisers' modelled Gln50 match the crystallographic Gln50 in 4BS0/7K4Z? *(M0)*
2. What is the exact anchor convention for the 12 Å burial count — carboxylate carbon, OD2, CA?
   Does it include the catalytic residue's own CA? Report alternatives and ask the organisers.
3. What is the exact Boltz-2 field name behind "iPDE", and is smaller better? Do not invent a
   formula. Read it from the organisers' scoring code if released.
4. Is the >70% identity threshold whole-chain or local? Use whole-chain and a conservative buffer
   until told otherwise.
5. Does RFD3 fit in 24 GB at length 250 with a ligand? *(M1 — measure, don't assume)*
