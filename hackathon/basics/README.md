# Basics (only things we can verify)

No Boltz, no invented cutoffs. Sequence gates use the organizer numbers. Geometry numbers are **measured on PDB 4BS0 chain A**, not DFT.

## Verified facts (from RCSB files)

| Item | Value |
|---|---|
| HG3.17 complex | `4BS0` |
| Catalytic base | **A ASP 127** |
| H-bond donor | **A GLN 50** |
| Scoring ligand in crystal | **6NT** = 6-nitrobenzotriazole, residue 1305 |
| Crystal contact (base) | ASP127 OD2 … 6NT N3 = **2.53 Å** |
| Crystal contact (donor) | GLN50 NE2 … 6NT N1 = **3.10 Å** |
| CA within 12 Å of ASP127 CA | **47** (excluding self) |
| True substrate in `7EPO` | ligand **H5J** = 5-nitro-1,2-benzoxazole (different pocket: Glu/His — not our default motif) |

Use `theozyme/4bs0_A_GLN50_ASP127_6NT.pdb` as the **teacher motif**. Do not submit 4BS0.

## Commands

```bash
# 1) sequence hard gates (length ≤260, identity ≤70%)
python3 hackathon/basics/scripts/validate_design.py --self-test
python3 hackathon/basics/scripts/validate_design.py your_designs.fasta

# 2) rebuild motif from the crystal (needs references/4BS0.pdb)
python3 hackathon/basics/scripts/extract_theozyme.py

# 3) distances on a complex (your residue numbers if not 50/127)
python3 hackathon/basics/scripts/measure_active_site.py hackathon/basics/references/4BS0.pdb
python3 hackathon/basics/scripts/measure_active_site.py design.pdb --base-res 12 --donor-res 40

# 4) burial count (slides: CA within 12 Å of the catalytic base)
python3 hackathon/basics/scripts/burial.py hackathon/basics/references/4BS0.pdb --base-res 127

# 5) TM > 0.7 gate — only if USalign/TMalign is installed (else exit 2, no fake score)
python3 hackathon/basics/scripts/tm_gate.py design.pdb
```

Redownload PDB/FASTA: `bash hackathon/basics/scripts/download_references.sh`

## What we still do not claim

- Official **numeric** pass/fail for geometry or burial (slides did not give one). We report vs 4BS0.
- iPDE / Boltz-2 consistency — need their scorer.
- 4BS0 is a **bound TSA**, not the DFT transition state. Geometry is a correct crystal reference, not the live TS.
