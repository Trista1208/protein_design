# Setup

Everything lives under `$HOME`. **No `sudo`, no `apt`, no system-wide installs.**
System `/opt`, `/usr/local` and `/usr/bin/gcc` are untouched, as are the ten
pre-existing conda envs (`her2_*`, `tdmpc*`, `alim-study`, `base`).

Hardware this was built and measured on:

```
GPU   1 x NVIDIA GeForce RTX 3090, 24 GB   driver 580.178.04  -> caps at CUDA 13.0
CPU   16 cores
RAM   92 GB
Disk  1.8 TB (~1.0 TB free after install)
```

There is **one** GPU. If you read "96 GB VRAM" somewhere, that was system RAM.

---

## 0. The one rule that will save you hours

> **Ubuntu 26.04's toolchain is newer than this scientific software expects.**
> Every build failure here had that one cause. Pin something older.

| Broke | Cause | Fix |
|---|---|---|
| Q6 | gfortran 16.2 rejects function-as-subroutine | conda gfortran **9.5.0**, and do **not** override the makefile's `FFLAGS` (it sets `PREC=-DQDOUBLE`) |
| Rosetta | GCC 15 `-Wdangling-pointer` inside bundled Boost Spirit + blanket `-Werror` | comment out the blanket `-Werror` |
| LigandMPNN | NumPy 2.5 removed `np.int` | patch the aliases |
| OpenMM | conda resolved CUDA **13.5** > driver 13.0 | pin `cuda-version=13.0` |
| RFdiffusion3 / RoseTTAFold3 | PyPI `rc-foundry 0.2.0` needs atomworks symbols in **no** released version | install `rc-foundry` **1.0.0** from git |

`-fallow-argument-mismatch` does **not** exist before GCC 10.

---

## 1. Environments

Conda is at `~/anaconda3` and is **not on PATH**. Use `~/anaconda3/bin/conda`.
All envs were created with `--override-channels -c conda-forge` (Anaconda's
`defaults` channel deliberately avoided).

| Env | Python | Purpose |
|---|---|---|
| `mm` | 3.12 | **the pipeline's home env** — OpenMM, pytraj (FieldTools), MDAnalysis, mdtraj, pdbfixer, openff-toolkit, RDKit |
| `foundry` | 3.12 | RFdiffusion3 (`rfd3`), RoseTTAFold3 (`rf3`), MPNN — `rc-foundry 1.0.0` |
| `aizymes` | 3.11 | ESMFold, ProteinMPNN, LigandMPNN — torch 2.6.0+cu124, fair-esm |
| `chai1` | 3.12 | Chai-1 0.6.1 |
| `her2_boltz2` | 3.12 | Boltz-2 2.2.1 (pre-existing) |
| `qm` | 3.12 | Psi4, PySCF, xtb, ASE — **ORCA + OpenMPI auto-loaded on activate** |
| `viz` | 3.12 | PyMOL |
| `q6` | — | gfortran 9.5.0 toolchain (build **and runtime** for Q6) |
| `rosettabuild` | 3.12 | scons, for building Rosetta only |

Recreate any env from `../envs/<name>.yml`. These are full `conda env export`
snapshots (exact builds) — they are lock files, not loose specs:

```bash
~/anaconda3/bin/conda env create -n mm -f ../envs/mm.yml
```

The pip-only pieces are additionally pinned in `../envs/pip-*.txt`.

### Non-conda installs, and why

Four things could not come from conda or uv and live in `~/opt`:

| Path | What | Install route |
|---|---|---|
| `~/opt/rosetta` | Rosetta source build (`relax`, `rosetta_scripts`, `match`) | public GitHub clone, academic licence, scons |
| `~/opt/orca` | ORCA 6.1.1, 17 GB | registered download from the ORCA forum (per-user EULA) |
| `~/opt/openmpi-4.1.8` | OpenMPI 4.1.8 | source build — conda-forge stops at 4.1.6 and ORCA needs **exactly** 4.1.8 |
| `~/opt/Q6` | Q6 EVB engine (`Qprep6`, `Qdyn6`, `Qfep6`) | git clone + `make all COMP=gcc` under gfortran 9.5.0 |

Each has an `env.sh`:

```bash
source ~/opt/rosetta/env.sh     # ROSETTA, ROSETTA3_DB, bin on PATH
source ~/opt/Q6/env.sh          # REQUIRED - binaries link conda's libgfortran.so.5
conda activate qm               # ORCA + OpenMPI come in via activate.d/orca.sh
```

---

## 2. Model weights

| Model | Location | Size |
|---|---|---|
| RFdiffusion3 | `~/.foundry/checkpoints/rfd3_latest.ckpt` | 2.5 G |
| RoseTTAFold3 | `~/.foundry/checkpoints/rf3_*.ckpt` | 2.8 G |
| RFD4-Proteina | HF cache (`flow-ema` + `ae-ema`) | 2.8 G |
| ESMFold | HF cache | ~2.6 G |
| Boltz-2 | `~/.boltz` | 7.6 G |
| ProteinMPNN | in-repo | 180 M |
| LigandMPNN | `tools/LigandMPNN/model_params/` | 118 M |

Fetch the foundry ones with:

```bash
~/anaconda3/envs/foundry/bin/foundry install rfd3 rf3
~/anaconda3/envs/foundry/bin/foundry list-installed
```

> **Never** run `snapshot_download("nvidia/RFD4-Proteina")` unfiltered. That repo
> is **107 GB** — it bundles AlphaFold2/OpenFold3 Apptainer images. Use
> `allow_patterns=["flow_checkpoints/*", "ae_checkpoints/*"]` (~2.8 GB).

---

## 3. Repo-local patches — a re-clone undoes these

1. **LigandMPNN** (`tools/LigandMPNN/`) — vendored OpenFold uses `np.int`,
   `np.bool`, `np.object`, removed in NumPy 1.24. Patched in
   `openfold/np/residue_constants.py`, `openfold/np/relax/utils.py`,
   `openfold/data/templates.py`. Word-boundary match, so `np.int64` is untouched.
   Its pinned torch 2.2.1 / numpy 1.23.5 are intentionally **not** honoured.
2. **Rosetta** (`~/opt/rosetta/source/tools/build/basic.settings`) — blanket
   `"Werror"` commented out; pristine copy at `basic.settings.orig`. Note
   `cxx_ver=12` is **invalid**; Rosetta accepts only `4.7`…`6.3` or `*`.

---

## 4. Verify the install

```bash
cd /home/mandalorian/Hackathons/Enzymes-BSSE2026
make verify          # every tool runs real work, not just imports
```

Import checks are near-worthless for this stack: `import rfd3` succeeded while
the engine was completely broken. `make verify` executes real code paths.

---

## 5. Licensing

- **Rosetta** — Non-Commercial Licence; ETH qualifies. Commercial use needs UW
  CoMotion (`license@uw.edu`). PyRosetta's licence is separate.
- **ORCA** — academic use only, registered per user. The 2025 EULA also restricts
  using ORCA output for **machine learning outside academic purposes**.
- **Q6** — © Åqvist, Marelius, Kamerlin, Bauer (Uppsala).
- **AI.zymes** — MIT-NC. We use only `FieldTools.py` from it.
- **RFD4-Proteina** — BSD-3 code, non-commercial generator weights.
