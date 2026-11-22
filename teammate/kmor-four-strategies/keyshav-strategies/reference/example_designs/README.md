# Kemp eliminase designs — shareable structures

De novo designs for Kemp elimination of 5-nitrobenzisoxazole, scored against
the transition-state analogue **6-nitrobenzotriazole (6NT)**.
DBSSE Protein Design Hackathon, Basel, November 2026.

Chain A = designed protein (150 aa). Chain B = 6NT ligand.
Catalytic pair: **Asp102** (general base), **Gln46** (oxyanion donor).

---

## The files

| file | what it is |
|---|---|
| `ins080_seed1_S0.743.pdb` | **the good one.** Official score **0.761** on this structure, all seven terms non-zero |
| `ins002_H103F.pdb` | best single mutant; excellent distances, dihedral in a dead zone |
| `d045_parent.pdb` | the parent both were derived from |

Lineage: `d045` → `H103F` → `ins080` (= d045 + **F42Y + H103F→Y**).

---

## Read this before quoting a number

**`ins080` scores 0.761 on one seed out of five.** The other four score
0.000. This is a real structure with a real score, not an average — but it is
not reproducible yet, and the honest summary is "0.761 achievable, 1/5
reproducibility", never "scores 0.76".

Per-seed, same sequence, Boltz-2, 5 independent seeds:

| seed | official | d_base Å | d_oxy Å | oxyanion dihedral | limiting term |
|---|---|---|---|---|---|
| **1** | **0.7613** | 2.644 | 2.886 | **48.2°** (in band) | orient_oxy |
| 2 | 0.0000 | 2.452 | 3.076 | 94.0° (dead zone) | orient_oxy |
| 3 | 0.0000 | 2.498 | 2.970 | 95.9° (dead zone) | orient_oxy |
| 4 | 0.0000 | 2.513 | 4.017 | 249.8° (dead zone) | oxy |
| 5 | 0.0000 | 2.488 | 4.107 | 245.9° (dead zone) | oxy |

The catalytic base contact is excellent and stable on every seed
(`d_base` 2.45–2.64 Å, sd 0.073 Å against HG3.17's 0.068). What varies is the
**Gln46 side-chain rotamer**: it samples three states, and only one of them
puts the oxyanion contact in a scoring geometry.

For scale: HG3.17 crystal structure scores 0.959 under this metric, and
HG3.17 folded with an MSA scores 0.834.

## Independent literature filters

Measured with filters from Risso 2020 (Chem. Sci. 11, 6134), Privett 2012
(PNAS) and Alexandrova 2008 (JACS), all for this reaction:

| | ins080 | 4BS0 crystal | reference |
|---|---|---|---|
| G1 donor–acceptor distance | **2.644 Å** | 2.532 Å | WT 2.64 Å |
| G2 D–H···A angle | 150.9° | 167.2° | WT 166.5°, soft ≥150 |
| H3 charged residues within 4 Å of base | **2 (Arg127)** | 0 | must be 0 |
| S2 apolar fraction, 5–6 Å shell | 0.500 | 0.727 | ≥0.60 |

`ins080`'s donor–acceptor distance is essentially the published wild-type
value. Its two known defects are **Arg127 salt-bridging the catalytic
aspartate** (the KE16 failure mode — removing the equivalent bridge in KE07
improved kcat/KM 2.5-fold) and an insufficiently hydrophobic pocket. Both are
being addressed.

## Caveats worth stating to anyone you send this to

- These are **predicted** structures (Boltz-2, single-sequence mode — these
  are de novo sequences with MSA depth 1, verified). No experimental
  validation.
- A score under this metric is a **geometric proxy**, not a rate constant.
  Nothing here establishes catalytic activity.
- 6NT is a transition-state *analogue*; the real substrate is
  5-nitrobenzisoxazole. Binding the analogue is not the same as catalysing
  the reaction.
