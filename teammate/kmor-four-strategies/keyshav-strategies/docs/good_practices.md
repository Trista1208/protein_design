# good_practices.md — read this BEFORE starting or confirming any task

Standing instructions for this project. Check this list before initiating work,
and again before reporting anything as done.

---

## 0. The pre-flight checklist

Before writing or running anything:

- [ ] Read `docs/LESSONS.md` and `docs/LOOP_DESIGN.md`. Where they disagree,
      LOOP_DESIGN wins (it corrects two measured errors in LESSONS §6).
- [ ] Check `DBSSE_Hackathon_Perlin.pdf` when anything about scoring,
      eligibility or winning criteria is in doubt. The slides are the
      authority — above this file, above the brief, above any paper.
- [ ] Check the supplied data files (§4) before inventing a threshold.
- [ ] Check `Claude_Kemp_Eliminase_Implementation_Brief.md` for an approach
      that already covers the problem.
- [ ] **Triple-check the assumption you are about to build on.** See §1.

---

## 1. Triple-check assumptions — the specific failures this prevents

Every one of these cost real time tonight. They are not hypothetical.

| What was assumed | What was true | How it was caught |
|---|---|---|
| "Rosetta is not built" | `rosetta_scripts` + `relax` built 21 Nov | `ls \| head -5` truncated the listing. The *command* was wrong, not the facts |
| "d045 has one zeroed term" | It had four | LESSONS §6 mixed OpenFold3 geometry with Boltz ipde |
| "F26Y fixes the dihedral" | Its ligand was 9.4 Å away — the angle was noise | Printing `d_base` alongside `d_oxy` |
| "V27I hits both terms" | `d_base` 7.2 Å on all 5 seeds | 5-seed score |
| "The params has N1/N2/N3" | Its `N1` was our NITRO nitrogen | Checking **bonded connectivity**, not name presence |
| "The dense rank includes ipde" | `ipde` was `None` on all 84 designs | The degeneracy alarm fired |
| "Theozyme Gln rotamer matches the crystals" | ~48° away; it is the only one scoring 0 | Measuring 4BS0 and 7K4Z directly |

**The rule distilled:** a negative result from a command you wrote is a claim
about your command, not about the world. Verify a negative before building a
plan on it. If a check would condemn a reference structure, the check is wrong
(LESSONS §2.1) — but if it *exonerates* the references and condemns only your
input, believe it.

---

## 2. Measurement discipline (non-negotiable)

- **Bound-ligand filter first.** Geometry is meaningful only when BOTH
  `d_base` and `d_oxy` are inside 2.0–3.6 Å. A dihedral on a drifted ligand is
  noise, not an angle (LESSONS §2.3). This has overturned three apparent
  breakthroughs. Never quote geometry before this filter.
- **pLDDT ≥ 0.95** (`config.PLDDT_MIN`). Below it a design is *unmeasurable*,
  not bad — keep it in a reserve tier, never rank on it.
- **pLDDT does not detect ligand loss.** d045 seed 4: pLDDT 0.967 with the
  ligand 9.7 Å away. Protein fold ≠ ligand placement.
- **Report every term, every design, every replicate** — especially the
  dihedrals (LESSONS §7.3). A term that is not in the table is a term that
  kills designs silently.
- **Five seeds or it is not a score.** One seed cannot compute `cons` at all.
- **SCORE EACH SEED, THEN AGGREGATE THE SCORES.** Never average the geometry
  and score once. Averaging a *bimodal* angular distribution produces a value
  no replicate holds. Measured both failure directions in one hour:
    - `ins079`: dihedral is ~99 deg or ~245 deg, never between. Circular mean
      138.8 deg lands inside the 130-230 band -> reported **0.624** while
      **0/5** replicates scored. A false positive, one command from being
      pushed to the team branch and sent to collaborators.
    - `ins080`: circular mean landed in a dead zone -> reported **0.000**
      while **1/5** replicates genuinely scored **0.743**. A real result,
      erased.
  Use `official.score_replicates()`. Quote `worst`/`median` with
  `n_nonzero_replicates`, e.g. "0.743 achievable, 1/5 reproducibility" -
  never a bare "scores 0.74".
- **Re-measure from the structure, not from stored replicate records.** Those
  records do not carry `angle2_oxy`; rebuilding a RawMetrics from them leaves
  the term MISSING rather than zero, and `score(strict=False)` then drops it
  and inflates the total. That is the LESSONS 1.5 failure, and it produced
  phantom 0.92 scores during the very check meant to catch the averaging bug.
- **Record the protocol with every number** (`objective.protocol_stamp`).
  MSA mode changed HG3.17 from 0.000 to 0.775.
- **Run the degeneracy alarm.** It has already caught one real bug.

## 3. Only a refold counts

Rosetta, LigandMPNN, a rotamer scan and OpenMM all produce **proposals**. The
competition metric is computed on a **Boltz-2 prediction**, so nothing is a
result until the sequence has been refolded.

Evidence: a coordinate scan found 360 clash-free in-band rotamers that 0/5
Boltz seeds sampled (LESSONS §6); Rosetta satisfied the enzdes constraints
(allcst 0.006) by sliding the ligand 1.35 Å, and the refold came back unbound.

**Corollary:** when using Rosetta, hold the ligand fixed
(`enzdes_ligfix.xml`). Geometry bought by moving the ligand is not
attributable to the sequence and cannot transfer.

## 4. The supplied data files — check before inventing a number

Live in `/home/mandalorian/Hackathons/`:

| File | Why it matters |
|---|---|
| `DBSSE_Hackathon_Perlin.pdf` | **The authority** on scoring, the 260 aa cap, the 70% identity and TM 0.7 limits, 5 replicates, and that 5 designs are accepted regardless of score |
| `filters hackaton.docx` | Filter definitions — check before setting any threshold |
| `filtri_backbone_v2.xlsx` | Backbone filters |
| `geometrie_kemp_risso2020_v2.xlsx` | Kemp catalytic geometry references |

Our open questions (dihedral band edges, `RAMP_MODE`, the burial anchor) must
be checked against these **before** any further inference.

## 5. Learn from the designs that worked

Keep the lineage visible; do not restart from scratch.

| design | what it taught |
|---|---|
| `rfd4_design_045` | the only parent that folds above the gate; 4 zeros |
| `ins002 H103F` | first bound-ligand variant with strong distances |
| `ins080 F42Y+H103Y` | **six of seven terms non-zero**, `cons` 0.853 at the HG3.17 level; only `orient_oxy` left. Found by COMBINING two singles that each failed alone |

The combination principle is the one that worked: single mutants fix one term
and break another; pairs did what neither could.

## 5b. A mechanism is a hypothesis until its control arm reports

The Tyr103 reading was right that a hydrogen bond pins the Gln46 rotamer and
wrong about which residue makes it. The correction came from designs where
Tyr103 was 6.9 A away and the dihedral was in band anyway.

Two habits that caught it:

- **Name the residue AND the general principle separately.** "Tyr103.OH
  H-bonds Gln46" is an observation on one structure; "an H-bond donor
  reaching the Gln46 amide pins the rotamer" is the claim worth testing.
- **Build the control arm into the library.** Block A deliberately includes
  non-donors at 103 so the mechanism can be falsified rather than confirmed by
  assumption. A library containing only the designs you expect to work cannot
  tell you anything.

When a mechanism survives correction, integrate it into the CODE, not just the
prose: `litfilters.gln_pinner()` measures it on every screen, so the next
version of this belief is testable against numbers rather than memory.

## 5c. A one-seed score is a candidate; say so when you report it

`L_B_anchor_M45Y` screened at 0.787 and was reported here as "beats ins080".
Five seeds returned 0/5. The claim was wrong the moment it was made, because
a single replicate cannot establish `cons` and cannot show reproducibility.

Report a screening hit as "0.787 on one seed, 5-seed confirmation running" —
never as a comparison against a design that HAS been confirmed. The two
numbers are not the same kind of thing.

## 5d. Never conclude from a best-of-N comparison

`supervisor.py` compared the BEST gap in a 9-design arm against the best in a
6-design arm, found 23.4 against 28.6, and reported the hypothesis
"SUPPORTED". The medians were 27.6 and 30.0, the distributions overlapped
almost completely, and neither arm had a single design in band.

The extreme of a larger sample is not a result. Compare **medians or full
distributions**, state n for each arm, and when neither arm reaches the
outcome being tested, the answer is **inconclusive** - the variable did not
decide anything, which is itself worth knowing.

This is the same error as quoting a one-seed score (§5c), in a different
costume: one favourable draw presented as the effect.

## 5e. Restart a long-running process after editing its code

Three times today a running job used a stale module: the campaign and
supervisor were blind to `litfilters` edits, and the publisher reported
"107 designs, best S=0" while 298 had been folded and one scored 0.83 - it
had been started before both the fix and the new screen paths.

Python loads modules at process start. Editing a file changes nothing for a
process already running it. After any edit that a running job depends on,
either restart that job or add a post-hoc pass that recovers the measurement
from data on disk (`annotate_screens.py` exists for exactly this, and saved
~70 folds).

## 5f. Never kill by a pattern your own shell matches

Four times today a `pgrep -f X` or `pkill -f X` matched the very shell running
it, because that shell's command line contains X. Twice it killed the launcher
before the job started; once it deadlocked a wrapper waiting for itself.

Get the PID first, then kill the PID:

    PID=$(ps -eo pid,cmd --no-headers | grep "[s]cript.py" | awk '{print $1}')
    [ -n "$PID" ] && kill "$PID"

The bracket trick (`[s]cript`) keeps grep's own line out of the match, and
killing a number cannot kill the wrong thing. Equally, count running jobs with
`ps -eo cmd | grep -oE ...`, never `pgrep -c`, which counts your shell too.

## 6. Tools — what is installed and what each is for

Use the right tool; do not substitute. FuncLib and PROSS are **not** options.

| Tool | Use it for |
|---|---|
| **Boltz-2** | the ONLY evaluator that counts |
| **LigandMPNN** | fold-compatible proposals. Forbid the native residue or it returns the parent (seq_rec 1.0) |
| **Rosetta** (`rosetta_scripts`, `relax`, `match`) | enzdes design scoring distance+angle+torsion jointly; ligand MUST be fixed |
| **OpenMM** | rotamer-stability questions; which neighbour pins a side chain |
| **Foldseek + MMseqs** | the novelty gate |
| **APNovo metrics** | the authoritative geometry measurement |
| **PyMOL / RDKit / OpenBabel** | inspection and chemistry prep |
| Chai-1, ESMFold, RoseTTAFold3 | secondary predictors — note the organisers score with Boltz |
| FieldTools, Q6/EVB, ORCA | mechanistic proxies; not in the metric |
| ProteinMPNN | ligand-blind — wrong tool for a catalytic pocket |

## 7. Keep the record current

- Add new findings to `docs/LESSONS.md` / `docs/LOOP_DESIGN.md` as they are
  measured, including retractions — a retraction is usually the more useful
  half.
- When a claim here is falsified, correct it here rather than silently
  working around it.

## 8. The floor

The organisers accept **at least five designs regardless of score**
(slide 7). Having zero submittable designs at the deadline is the only
unrecoverable failure. Bank five before chasing a better score.

### 5g. A claim about your designs deserves a measurement, not a rebuttal

When someone says the designs have a problem, the useful reply is a number with
a **control arm**, not an argument. The disulfide question (LOOP_DESIGN §6.12)
is the worked example: the teammate's premise about the reference enzyme turned
out to be *right* and my assumption wrong — the control is what caught it — and
the conclusion was still wrong, which only the per-residue comparison showed.

Check the premise before you dispute the conclusion.
