# tam.py: run RFdiffusion3 or RFD4-Proteina on Tamarind

Sources: official RFdiffusion3 docs (RosettaCommons/foundry, `models/rfd3/docs/input.md`), official RFD4-Proteina docs
(rfd4-proteina.org, `docs/reference/yaml-json.md`), and Tamarind's tool schema (`GET /api/tools`) and API guide.

## 1. Quick start

1. Put your Tamarind key in `settings.json` (`api_key`). `base_url` stays as it is.
2. Check one of the run files: `run_rfd3.json` or `run_rfd4.json`.
3. Validate first (free, no job created), then run:

```
python3 hackaton_scripts/tam.py hackaton_scripts/run_rfd3.json validate
python3 hackaton_scripts/tam.py hackaton_scripts/run_rfd3.json run        # validates again, submits ONLY if valid (costs credits)
```

Same commands with `run_rfd4.json` for RFD4-Proteina. `run` submits the exact payload that was just validated.

## 2. The three files

| File | You put in it |
|---|---|
| `settings.json` | `base_url`, `api_key` |
| `run_*.json` | `type` (tool), `jobName`, `input` (path to the design input file), `settings` (hyperparameters, see sections 3 and 4) |
| the design input file (`input`) | the design itself: structure, ligand, motif, lengths... |

`tam.py` never invents values: everything sent is the design input file, plus the `settings` you wrote, plus the uploaded file names.
The structure file named inside the input file must sit in the same folder as the input file.

## 3. RFdiffusion3 (`"type": "rfdiffusion3"`)

**Official input:** a JSON/YAML file in RFdiffusion3's `InputSpecification` format (`input`, `contig`, `unindex`, `length`, `ligand`,
`select_fixed_atoms`, `ori_token`, `allow_ligand_on_existing_chain`, ...). `test_pipeline/1qji_json.json` is already in this format.

**On Tamarind:** task `json`. The file is uploaded unchanged (`jsonFile`) together with its structure (`inputFiles`), so every field in it is
interpreted by RFdiffusion3 itself, including `ori_token`.

`run_rfd3.json` settings:

| Setting | Meaning | Tamarind default |
|---|---|---|
| `task` | must be `json` for this route | |
| `jsonInputType` | `file` | |
| `numDesigns` | total designs generated | 8 |
| `diffusionBatchSize` | designs per GPU pass | 8 |
| `stepScale` | diffusion step size; lower = more diversity, higher = more designable | 1.5 |
| `gammaZero` | extra noise per step; lower = more designable, less diverse | 0.6 |
| `claude_inference_optimization` | `off`, `exact` or `fast` inference kernels | `off` |

The official RFdiffusion3 CLI has more options (e.g. `num_timesteps`). Tamarind's schema does not expose them.

## 4. RFD4-Proteina (`"type": "rfd4-proteina"`)

**Official input:** a YAML/JSON *condition spec*: one named design with `input` (a `.cif`), `contig`, and `conditions`. Conditions are
annotations such as `C_CRD` (fix coordinates), `C_SEQ` (fix sequence), `C_IDX` (unindexed motif), `C_CENTROID` (chain centroid,
the counterpart of RFdiffusion3's `ori_token`), `C_HOT`, `C_RAS`, `C_HB`. Full list: rfd4-proteina.org, "YAML/JSON condition-spec parameters".
Residue numbers in the spec refer to the `label_seq_id` column of the mmCIF, not `auth_seq_id`.

**On Tamarind there are two routes:**

| Route | Tamarind task | What you provide |
|---|---|---|
| A. Native (not implemented in `tam.py`) | `custom-conditions` | the condition spec file (`conditionSpec`) + the `.cif` |
| B. Form settings (what `tam.py` does now) | `enzyme-design` | `tam.py` translates an RFdiffusion3-style input file into Tamarind's form settings |

Route B translation: `input`→`pdbFile`, `ligand`→`ligandCode`, `length`→`minScaffoldLength`/`maxScaffoldLength`,
`unindex`→`motifResidues`, `select_fixed_atoms`→`catalyticAtomPicks`.
**Route B cannot express `ori_token`/`C_CENTROID`, so it is not sent.** Route B is not an official RFD4 input format; route A is.

`run_rfd4.json` settings (Tamarind schema; anything omitted gets Tamarind's default):

| Setting | Meaning | Default |
|---|---|---|
| `task` | `enzyme-design` (route B) | |
| `numDesigns` | number of designs | 8 |
| `samplingEffort` | `thinking_low`, `thinking_medium`, `thinking_high` | `thinking_medium` |
| `samplingSteps` | sampling steps per design, 60 to 120 | 60 |
| `seed` | repeat a run | none |
| `searchPreset` | `off`, `interface-contacts`, `low-clash`, `oracle-final` | `off` |
| `searchBeamWidth` / `searchBranches` | reward-search breadth, 1 to 4 | 1 / 2 |
| `deterministicSampling` | sample without added noise | false |
| `backboneTemperature` | 0.01 to 1.5, higher = more varied folds | 1 |
| `sequenceTemperature` | 0.01 to 5, higher = more varied sequences | 1 |
| `tsrBackboneK`, `tsrBackboneSigma`, `tsrLatentK`, `tsrLatentSigma` | temporal score rescaling | 16, 3, 4, 12 |
| `catalyticAtoms` | which atoms of catalytic residues stay fixed: `tip-atom`, `functional-tip`, `sidechain`, `backbone` | `functional-tip` |
| `unindexMotif` | model chooses where the motif sits | true |
| `fixMotifSequence` | keep the motif's amino acids | true |
| `ligandHandling` | `fixed` or `sampled` ligand pose | `fixed` |
| `excludeAminoAcids` | amino acids to avoid (a bias, not a ban) | none |
| `refoldDesigns` | re-predict each design with OpenFold3 and compare | true |
| `preOracleFilters` | drop designs that break a condition before scoring | true |
| `passMotifRmsd`, `passBackboneRmsd`, `passInterfacePae` | pass limits (Å) for refolded designs | authors' limits |
| `clusterDesigns`, `developabilityMetrics`, `saveTrajectory` | extra outputs | false |

## 5. What you see in the terminal

| Step | What Tamarind returns |
|---|---|
| `upload` | the upload response (a `308 Redirecting...` is followed automatically) |
| `validate` | `{"valid": ..., "normalized": {...}}`. `normalized` is your payload plus the defaults Tamarind fills in. No design is computed, so there are no design outputs here. |
| `run` (submit) | a short confirmation, e.g. `<jobName> submitted to queue.` The designs are produced later by the job itself. |

Validation and submission use the same validator, so a payload that validates will be accepted by `run`.

## 6. Not handled / not verified

- Results are not polled or downloaded. Check the job on the Tamarind website.
- Route A (`custom-conditions`) is not implemented, and how Tamarind maps the spec's `input` path to an uploaded `.cif` is not verified.
- Uploading a file whose name already exists in your account: whether it is replaced is not verified. Rename the file if its content changed.
- `inputFiles` as a list of names and `jsonFile` as a plain file name (RFdiffusion3) have not been tested against Tamarind.
