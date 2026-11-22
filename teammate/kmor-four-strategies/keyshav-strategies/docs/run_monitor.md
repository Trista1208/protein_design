# run_monitor.md — check this every 15 minutes, without exception

Standing rules for the autonomous run. Read alongside `good_practices.md` and
`LESSONS.md`. Where this file and a convenience disagree, this file wins.

---

## 1. The authority is the PDF, not this pipeline

`DBSSE_Hackathon_Perlin.pdf` defines the challenge. Re-check against it
**every time a new design produces a new score**. The binding facts:

| | |
|---|---|
| Catalytic base | **ASP** |
| H-bond donor | **GLN** |
| Transition-state analogue | **6NT** (6-nitrobenzotriazole) — the scoring ligand, an approximation of the DFT transition state |
| Replicates | **exactly 5 Boltz-2 seeds per sequence** |
| Max length | **260 aa** |
| Sequence identity to published enzymes | **>70% disqualifies** |
| Backbone TM-score to published structures | **>0.7 disqualifies** |

### The four scoring criteria, in the organisers' own terms

1. **iPDE** — a Boltz-2 metric estimating error in protein–ligand token
   distances. Confidence in ligand placement. **Always taken from Boltz-2
   output, never synthesised.** Field: `complex_ipde`; it is an ERROR, so
   smaller is better.
2. **Geometry** — distances, angles and dihedrals between the catalytic base,
   the H-bond donor and the ligand, against ideal catalytic geometry.
3. **Catalytic base burial** — Cα atoms within 12 Å of the catalytic base.
   Buried rather than solvent-exposed, to preserve its pKa.
4. **Consistency** — reproducibility of the metrics across the Boltz-2 seeds.
   A proxy for active-site preorganization.

**Never substitute a proxy for one of these four.** The literature filters
(Risso/Privett/Alexandrova) and the dense objective are *internal* instruments
for ranking and diagnosis. They do not replace the official metric and must
never be reported as if they did.

---

## 2. Every 15 minutes, in order

- [ ] Re-read `good_practices.md`, `LESSONS.md`, and this file.
- [ ] Check the run: `logs/library_screen.log`, `logs/supervisor.log`,
      `logs/live_publish.log`, `logs/annotate.log`,
      `outputs/supervisor_state.json`.
- [ ] Confirm the publisher is alive and **pushing**; restart it if dead.
- [ ] Report the current best score **per replicate, with n/5
      reproducibility**. Never a bare average, never a single-seed number
      presented as a result.
- [ ] Triple-check any new claim by **re-measuring from the structures**, not
      from stored numbers.
- [ ] Record any verified finding in `LESSONS.md` or `LOOP_DESIGN.md`; record
      any new practice here or in `good_practices.md`.

## 3. On a serious bug mid-run

1. **Do not discard completed folds.** They are on disk; a measurement can
   almost always be recovered post-hoc (`annotate_screens.py` exists because
   of exactly this).
2. Consult `good_practices.md` §1 for the known failure classes before
   diagnosing: truncated commands reported as negatives, falsy zeros,
   stale modules in running processes, self-matching `pgrep`/`pkill`.
3. Fix the code, then **pause gently and resume** — kill by PID, never by a
   pattern that matches the killing shell's own command line.
4. Re-annotate rather than re-fold wherever the data allows it.

## 4. What gets published, and for whom

To `kmor/four-strategies`, every cycle:

- `designs.csv` — every measured design, all metrics, good and bad
- `README.md` — plots embedded, ✅ scoring and ❌ zeroed sections, the
  limiting term for each failure
- `structures/` — PDBs
- `gifs/` — **rotating 3D GIFs for designs that meet the criteria**
- `plots/`, `status.json`

**Publish partial passes too.** A design that scores on fewer than 5
replicates is a real result and goes up with its PDB, its GIF and its
`n/5` reproducibility stated. Suppressing it would hide the only genuine
result this run has produced (`ins080`, 0.7435 on 1/5).

## 5. Eligibility is checked, not assumed

Before anything is called submittable:

- length ≤ 260 (`gates.check_length`)
- sequence identity ≤ 0.70 to published enzymes (`gates.check_identity`,
  MMseqs2)
- backbone TM ≤ 0.7 to published structures (`novelty.screen`, Foldseek)

A design failing any of these is **disqualified**, regardless of its score.
Re-run the gate rather than trusting an earlier verdict after a sequence
changes.

## 6. Two failure modes, kept separate

| mode | signature | fix |
|---|---|---|
| **A. Rotamer flip** | distances fine, oxyanion dihedral in a dead zone | the Gln46 pinner — an H-bond donor within ~2.7 Å of the amide |
| **B. Ligand loss** | one seed with d_base > 4 Å; `cons` collapses to 0 for every replicate | pocket packing; rotamer work does not touch it |

Rank on **how many seeds reach full scoring geometry**, not on the best
single seed.
