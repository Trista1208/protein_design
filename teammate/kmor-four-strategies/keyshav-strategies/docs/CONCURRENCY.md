# How much can run at once?

**Short answer: run the GPU stages sequentially. Overlap CPU work freely.**

There is **one** GPU — a single RTX 3090 with 24 GB. There is no 96 GB device on
this machine; the 92 GB figure is system RAM. Everything below follows from that.

```
GPU   1 x RTX 3090, 24576 MiB   driver 580.178.04 (caps at CUDA 13.0)
CPU   16 cores
RAM   92 GB
```

Roughly **1.2 GB of VRAM is already taken by the desktop session**, so plan
against **~23 GB usable**, not 24.

---

## Measured, not estimated

Each tool run alone, peak VRAM sampled at 0.4 s intervals, on a short test input
(33–50 aa + ligand). `make vram` reproduces this table.

| Tool | Peak Δ VRAM | Wall clock | Notes |
|---|---|---|---|
| **ESMFold** | **14.0 GB** | 27 s | by far the heaviest; fp16 ESM trunk |
| **Chai-1** | **8.2 GB** | 73 s | includes ESM embedding pass |
| **RoseTTAFold3** | ~3.2 GB | 17 s+ | measured on a partial run; re-measure with a full fold |
| **Boltz-2** | **2.4 GB** | 95 s | cheapest folder per GB; 24 GB is ample (handover 1.6) |
| **RFdiffusion3** | ~3.3 GB observed | first run dominated by Triton JIT | see caveat below |
| **LigandMPNN** | **0.55 GB** | 4 s | effectively free |
| **ProteinMPNN** | ~0.5 GB | <1 s | effectively free |
| Rosetta / Q6 / ORCA / Psi4 | **0 GB** | — | CPU only |

### Two caveats that matter more than the numbers

1. **These are lower bounds.** The test inputs were 33–50 residues. Real designs
   are **190–260 aa plus a ligand**. Attention memory grows super-linearly with
   length, so expect materially higher peaks. **Re-run `make vram` with a
   250-residue sequence before trusting any schedule.**
2. **First RFD3 run is not representative.** It JIT-compiles Triton kernels,
   which dominates wall clock and inflates apparent latency. Measure the
   *second* run.

---

## What fits together

Against ~23 GB usable:

| Combination | Σ VRAM | Verdict |
|---|---|---|
| LigandMPNN + anything | +0.55 GB | ✅ always safe — treat as free |
| Boltz-2 + Chai-1 | 10.6 GB | ✅ comfortable |
| Boltz-2 ×2 (two seeds) | 4.8 GB | ✅ comfortable |
| Boltz-2 ×4 | 9.6 GB | ✅ fits, but see the CPU/IO caveat |
| RFD3 + Boltz-2 | ~5.7 GB | ✅ fits at test size; unproven at 250 aa |
| ESMFold + Boltz-2 | 16.4 GB | ⚠️ fits with ~6 GB headroom |
| ESMFold + Chai-1 | 22.2 GB | ❌ **do not** — under 1 GB spare, will OOM on real lengths |
| ESMFold + anything at 250 aa | — | ❌ assume ESMFold needs the card alone |

**ESMFold is the constraint.** At 14 GB on a 50-mer it already occupies 60% of
the card. On a 250-mer it should be assumed to need the GPU to itself.

---

## Recommended schedule

### Sequential on the GPU, parallel on the CPU

```
      ┌─ GPU lane (one at a time) ─────────────────────────────┐
M1/M4 │ RFD3 backbones                                        │
M2    │ LigandMPNN  (0.55 GB - can piggyback on anything)      │
M2/M4 │ Boltz-2 prefilter x1                                   │
M4    │ Boltz-2 x5 on finalists only                           │
      └────────────────────────────────────────────────────────┘
      ┌─ CPU lanes (16 cores, run alongside freely) ───────────┐
M0    │ theozyme audit, geometry, gates            (seconds)   │
M3    │ scoring, novelty gates (Foldseek/MMseqs2)  (minutes)   │
M5    │ FieldTools electrostatics                  (minutes)   │
      │ Rosetta relax / rosetta_scripts            (14 s each) │
      │ Q6 EVB, ORCA, Psi4                         (CPU-bound) │
      └────────────────────────────────────────────────────────┘
```

**Why sequential and not 2–3 GPU jobs?** Even where VRAM fits, concurrent CUDA
contexts on one consumer card contend for SM time and PCIe bandwidth. Two jobs
that each take 90 s alone typically take ~200 s together — slower than running
them back to back, and with OOM risk added. You gain throughput only when a job
*cannot* saturate the GPU, which is not the case for any folder here.

### The one genuinely free parallelism

**LigandMPNN at 0.55 GB and ~4 s is noise.** Run sequence design for batch *N*
while Boltz folds batch *N−1*. That is real pipelining with no contention.

### CPU saturation is the other real win

Rosetta, Q6, ORCA and Psi4 touch no VRAM. With 16 cores you can run the entire
scoring and novelty-gating stack for finished designs while the GPU keeps
generating. Reserve ~2 cores for the GPU jobs' data loaders — RFD3 uses
multi-worker loaders and starving them stalls the GPU.

---

## Budget arithmetic for M4

Using the measured Boltz-2 figure (95 s on a short input; assume **2–4 min** at
250 aa with a ligand) and 5 replicates per finalist:

- 1 finalist ≈ **10–20 GPU-minutes**
- One 12 h overnight run ≈ **36–72 finalists** *if Boltz had the card to itself*
- It does not — RFD3 and LigandMPNN need the same GPU

**So:** Boltz ×1 as a cheap prefilter over a few hundred designs, then Boltz ×5
on **20–40 finalists only**. This matches handover §4.1, now with measured
numbers behind it.

Do not run 5-replicate folding over hundreds of designs. You will not finish,
and the handover's warning stands: measure on your own inputs and re-plan.

---

## If a second GPU appears

The handover mentions "optionally one more cluster GPU". The natural split:

| Device | Job |
|---|---|
| Bigger card (≥40 GB) | ESMFold and RFD3 at full length |
| RTX 3090 | Boltz-2 prefilter + LigandMPNN pipelining |

Nothing in `kemp/` assumes one device — the runners take explicit output
directories and are safe to invoke per-GPU with `CUDA_VISIBLE_DEVICES`. There is
no multi-GPU scheduler; use one process per device.

---

## Reproducing this table

```bash
make vram      # sequential, one tool at a time, 0.4 s sampling
```

Writes `outputs/vram/summary.txt` plus per-tool logs and raw memory traces.
Edit `scripts/measure_vram.sh` to swap in a 250-residue sequence — **do this
before committing to an overnight schedule.**
