# User guide

Prerequisite: `docs/SETUP.md` is done and `make doctor` is clean.

```bash
cd /home/mandalorian/Hackathons/Enzymes-BSSE2026
make help         # all targets
make doctor       # envs, ~/opt installs, weights
```

Everything is driven by `make`, so the shortest possible path is:

```bash
make m0           # audit + controls
make m1           # RFD3 smoke, measures VRAM
make m3           # scoring + gates self-test
make m4           # refuses until M0 controls pass - by design
```

or `make all` for the gated chain. Logs land in `logs/<milestone>_<timestamp>.log`;
machine-readable reports in `outputs/runs/<milestone>/report.json`; cross-milestone
state in `outputs/state.json`.

---

## M0 — audit and controls

```bash
make m0
```

Does six things:

1. Loads both theozyme files **safely** — altloc normalisation that keeps the
   highest-occupancy conformer *per atom* rather than preferring altloc `A`.
   Preferring `A`/blank silently deletes the ligand, which carries altloc `B`.
2. Proves the two-frames story. Expect `superposed RMSD 0.0005 Å`, `raw RMSD 68.7 Å`.
3. Checks the `ORI` marker is within 12 Å of the ligand **in each frame**, and
   asserts that deliberately mixing frames is *detected*.
4. Re-measures the active-site geometry and diffs it against the verified table.
   All six values should agree to ≤0.05.
5. Writes `outputs/runs/m0/theozyme_clean.pdb` — single conformer, **no ORI atom**.
   Everything downstream uses this, never the raw file.
6. Reports catalytic burial under **all three** anchor conventions.

### The control pair — do not skip this

`make m0` reports controls as `PENDING`, and **`make m4` refuses to run** while
that is true. To complete it:

1. Fetch the HG3.17 sequence (`4BS0`) and HG3.17-Q50M (`6TU6`).
2. Run both through `boltz_replicates` (5 seeds) and the full scorer.
3. Feed both into `kemp.scoring.NoiseFloor` and check `.verdict()`.

If the positive does not outrank the negative, **stop**. Your Geometry and
Consistency metrics are not detecting oxyanion stabilisation and nothing produced
afterwards means anything. The spread across the five seeds is your **noise
floor**: two designs closer than that are indistinguishable.

```python
from kemp.scoring import NoiseFloor, aggregate
nf = NoiseFloor(positive=aggregate(hg317_scores), negative=aggregate(q50m_scores))
print(nf.verdict())
```

---

## M1 — RFD3 smoke test

```bash
make m1                       # default 20 timesteps
```

Builds a correct enzyme-design spec via `kemp.rfd3.build_spec` and runs one
batch. Use it to get **your** VRAM and wall-clock numbers — the RFD3 docs state
neither.

Writing a spec by hand:

```python
from kemp import rfd3
from kemp.theozyme import load
tz   = load("outputs/runs/m0/theozyme_clean.pdb")
spec = rfd3.build_spec("canonical", tz, length="190-250")
js   = rfd3.write_inputs_json([spec], "outputs/runs/m1/inputs.json")
rfd3.run(js, "outputs/runs/m1/designs", batch_size=1, timesteps=200)
```

**Schema notes** (read from the installed `rfd3/inference/input_parsing.py`, not
from prose):

- A composition field is **mandatory**. `input` alone raises
  `Input provided but unused in composition specification`.
- `ori_token` is a JSON field of three floats — *not* an atom in the PDB.
- `select_hbond_donor` / `select_hbond_acceptor` do atom-wise H-bond
  conditioning. Oxyanion stabilisation *is* an H-bond, so the `canonical`
  variant conditions `GLN:NE2` as a donor.
- `diffusion_batch_size=1` because all designs in a batch share one length.
- Never fix side-chain atoms of a residue listed in `select_unfixed_sequence` —
  `build_spec` raises if you try.

Available motif variants (`kemp.config.MOTIF_VARIANTS`): `canonical`,
`crystal_gln`, `backbone_nh`, `alt_donor`, `base_variant`. **Build at least three.**

---

## M2 — sequence design and folding

```python
from kemp.runners import ligandmpnn, assert_catalytic_preserved, boltz_replicates

res = ligandmpnn("design.pdb", "outputs/runs/m2/lmpnn",
                 fixed_residues="A127 A50", n_seqs=8)
assert_catalytic_preserved(res.sequences, {126: "D", 49: "Q"})   # 0-based

for header, seq in res.sequences[1:]:
    boltz_replicates(seq, f"outputs/runs/m2/boltz/{header}")
```

`--fixed_residues` is not optional in practice: without it LigandMPNN mutates the
catalytic aspartate away. `assert_catalytic_preserved` then fails loudly rather
than letting a dead design through.

Boltz settings are **pinned** in `kemp.config.BOLTZ_SETTINGS` and the runner will
not parameterise them. Drifting MSA settings between designs makes confidences
incomparable, and MSA depth is called out on slide 7.

---

## M3 — scoring and gates

```bash
make m3
```

Self-test: the reference theozyme must score ~1.0 on geometry (it scores 0.999),
a 240-aa sequence must be `ELIGIBLE`, a 300-aa one `DISQUALIFIED`.

```python
from kemp import scoring, gates

geo = scoring.geometry_score(measured)          # dict + 0-1 aggregate
bur = scoring.catalytic_burial("design.pdb")    # all three anchors
key, ipde = scoring.extract_ipde(conf)          # raises if absent - never invents
agg = scoring.aggregate(values)                 # median minus spread
con = scoring.consistency(values)

g = gates.run_all(seq, structure=pdb, corpus_dir=refs,
                  catalytic={126: "D", 49: "Q"})
assert g["passed"]
```

Three deliberate behaviours:

- `extract_ipde` **raises** when no recognised confidence field exists. The exact
  iPDE field name is unconfirmed; guessing a formula would silently fabricate
  scores. Pass `key=` once the organisers publish theirs.
- `aggregate` has **no `best` mode**. Ranking on the best seed defeats
  Consistency, which is an official metric.
- Missing tools or empty corpora make novelty gates **fail closed** — unproven
  novelty counts as failure, not as a pass.

---

## M4 — scale (overnight)

Gated on M0. Per batch: sweep ORI offsets × motif variants × lengths → backbones
→ 8 LigandMPNN sequences each → hard filter → **Boltz ×1** as a cheap prefilter →
**Boltz ×5** on 20–40 finalists only.

Budget from the measured numbers in `docs/CONCURRENCY.md` before launching. Do
not run five-replicate folding on hundreds of designs; you will not finish.

---

## M5 — select the final five

Cluster for **diversity** first, then tie-break:

1. **ESM2 pseudo-log-likelihood** — proxy for whether it will express. The
   official scoring rewards nothing here, and five designs that don't express are
   five empty wells.
2. **FieldTools** — electrostatic preorganisation. Report as a **proxy**, never
   as ΔG‡.

```python
from kemp.runners import fieldtools
fieldtools("final_design.pdb", "outputs/runs/m5/field")   # runs under the mm env
```

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `ORI is 68.7 Å from the ligand centroid` | mixed coordinate frames | use one file's frame throughout; `ori_for()` refuses to guess |
| `6NT vanished during parsing` | altloc filter dropped the ligand | use `kemp.theozyme.load`, not a naive parser |
| `Input provided but unused in composition specification` | RFD3 spec has no `length`/`unindex`/`contig` | use `build_spec` |
| `Unsupported file type: .fasta` | `rf3` wants `.json/.yaml/.cif/.pdb` | convert, or pass a PDB |
| `no recognised interface-confidence field` | Boltz output lacks the key | pass `extract_ipde(conf, key=...)` |
| `LigandMPNN mutated catalytic residues` | `--fixed_residues` missing | pass it; the assert is working as intended |
| `np.int` AttributeError in LigandMPNN | repo re-cloned, patch lost | reapply the NumPy alias patch (SETUP.md §3) |
| OpenMM `CUDA_ERROR_UNSUPPORTED_PTX_VERSION` | conda CUDA newer than driver | `conda install -n mm cuda-version=13.0` |
| Q6 binaries won't start | missing libgfortran | `source ~/opt/Q6/env.sh` |
| Two conda commands at once corrupt an env | known conda bug | **never** run two conda transactions concurrently |
