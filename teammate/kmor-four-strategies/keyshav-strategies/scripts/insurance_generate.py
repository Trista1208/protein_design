"""Generate the d045 insurance family: parent-biased, few-edit variants.

PURPOSE, stated plainly: we have ZERO banked submissions and one parent that
folds above the pLDDT gate. The organisers accept five designs regardless of
score (slide 7). This produces candidates that are compliant and foldable, and
which additionally have a shot at the one defect worth fixing.

THE METHOD, and where each piece comes from
-------------------------------------------
* Perturb a foldable parent, never write a whole sequence
  (LESSONS 7.6; HER2/Catalyst in LESSONS 4.1). The composition pathology of
  LESSONS 3.7 - 24-27% Ala, zero Trp, pLDDT 0.41-0.47 - is structurally
  impossible when 97% of the chain is frozen.
* Open only a FEW positions at a time, so the edit count is bounded BY
  CONSTRUCTION rather than filtered for afterwards. AI.zymes caps mutational
  load and biases toward the input sequence for the same reason.
* Target the positions that contact the Gln46 side-chain tip. LOOP_DESIGN 6.2:
  `oxy` and `orient_oxy` are the same tip (d_oxy = NE2...N1, orient_oxy =
  CD-NE2...N1-N2), together worth 0.28 of the score, and both are currently
  0.000. A second-shell change that repacks that rotamer moves both at once.
* Catalytic Asp102/Gln46 are NEVER designable, and the result is asserted
  afterwards (handover 4.3.1 - LigandMPNN will mutate them given the chance).

WHAT THIS DOES NOT DO
---------------------
It does not claim any variant is better. It produces candidates; only a refold
scores them (LESSONS 6: 0/5 seeds sampled any of the 360 clash-free in-band
rotamers found by a coordinate scan, so coordinate reasoning does not survive).
Stage 2 (`insurance_fold.py`) folds and scores them.
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C                                   # noqa: E402
from kemp.runners import RunnerError, ligandmpnn               # noqa: E402

#: ABSOLUTE: LigandMPNN shells out to prody, which resolves a relative path
#: against its own cwd and then reports it as "not a valid PDB identifier".
PARENT_PDB = (Path(__file__).resolve().parents[1]
              / "outputs/refold_d045/seed2/folded.pdb")
PARENT_SEQ = ("SSFLELMHKAVSAAKQGDFTSAISDFVEACNLLLPDDKKDDFEAMQKALLSTLLAHHDQDRAKDALTTLLTVTLILFE"
              "KHHDDVVELMQKHYDHFHDMVSDDHKIVTEYLKGELDYEDYLTSLVHNRGDALNLSEAEITEIYEEAKKMLS")
CHAIN = "A"
BASE_IDX, OXY_IDX = 102, 46

#: Positions within ~5.1 A of the Gln46 side-chain TIP, measured on the parent
#: (scripts/insurance_shell.py). These are the residues that can physically
#: repack the oxyanion rotamer. Ordered by proximity.
ROTAMER_SHELL = [26, 103, 42, 99, 49, 45, 23, 27]

#: A few positions further out, opened only to create sequence diversity for
#: the novelty gate without touching the pocket. These are the true fallback:
#: if nothing fixes the rotamer we still need five distinct, compliant designs.
DIVERSITY_SHELL = [66, 70, 50, 53, 84, 87]

OUT = Path("outputs/insurance")
MAX_EDITS = 2

#: MEASURED, and the reason this script works at all: with the native residue
#: allowed, LigandMPNN returned the parent at 16 of 20 open sets and the whole
#: run produced 4 variants. Its reported `seq_rec` is 1.0000 - d045's sequence
#: is already at this model class's optimum, because it came from one. Asking
#: "what belongs here?" returns what is already there.
#:
#: So the native residue is FORBIDDEN at each open position
#: (--omit_AA_per_residue). That turns a proposal into a substitution.
#: Cysteine is forbidden everywhere: a free thiol in a designed pocket is a
#: liability nobody asked for.
OMIT_EVERYWHERE = "C"
TEMPERATURE = 0.3          # above the 0.1 default; natives are omitted anyway
N_SEQS = 8
N_BATCHES = 2              # 16 proposals per open set


#: A steric ladder, small -> large, for the rotamer-shell positions.
#:
#: WHY THIS EXISTS ALONGSIDE THE MODEL PROPOSALS, and it is the more honest
#: half: LigandMPNN optimises sequence-backbone compatibility. It has no
#: knowledge of which Gln46 chi2/chi3 rotamer BOLTZ will choose, which is the
#: only thing that moves `oxy` and `orient_oxy`. Its first choice is the
#: native (seq_rec 1.0) and forbidding the native merely yields its SECOND
#: choice - a slightly less likely residue, not one selected for our defect.
#:
#: Our defect is mechanical and specific: d_oxy misses the 3.6 A window by
#: 0.054 A and the dihedral misses the band by 3.4 deg (LOOP_DESIGN 6.2). A
#: neighbour that is too large crowds the Gln tip; one that is too small fails
#: to pin it. This ladder spans both directions explicitly, which is a
#: proposal aimed at the objective rather than at sequence likelihood.
#:
#: No claim is made that these are good sequences. Only a refold decides
#: (LESSONS 6); the pLDDT gate catches the ones that break the fold.
STERIC_LADDER = ["A", "V", "L", "I", "F", "W"]


def ladder_variants() -> dict[str, dict]:
    """One-edit variants spanning a size ladder at each rotamer-shell position."""
    out: dict[str, dict] = {}
    for pos in ROTAMER_SHELL:
        native = PARENT_SEQ[pos - 1]
        for aa in STERIC_LADDER:
            if aa == native:
                continue
            s = PARENT_SEQ[:pos - 1] + aa + PARENT_SEQ[pos:]
            out[s] = {"tag": f"ladder{pos}", "open": [pos],
                      "edits": [f"{native}{pos}{aa}"], "n_edits": 1,
                      "source": "steric_ladder"}
    return out


def edits(seq: str) -> list[tuple[int, str, str]]:
    """1-indexed differences from the parent."""
    return [(i + 1, PARENT_SEQ[i], seq[i])
            for i in range(len(PARENT_SEQ)) if seq[i] != PARENT_SEQ[i]]


def redesign_spec(designable: list[int]) -> str:
    """LigandMPNN `--redesigned_residues`: only these are open."""
    return " ".join(f"{CHAIN}{r}" for r in designable)


def omit_spec(designable: list[int]) -> dict[str, str]:
    """Forbid the NATIVE residue at each open position, plus Cys everywhere.

    Without this the model returns the parent (see OMIT_EVERYWHERE above).
    """
    return {f"{CHAIN}{r}": PARENT_SEQ[r - 1] + OMIT_EVERYWHERE
            for r in designable}


def main() -> int:
    if not PARENT_PDB.exists():
        print(f"parent structure missing: {PARENT_PDB}")
        return 2
    n_res = len(PARENT_SEQ)
    assert PARENT_SEQ[BASE_IDX - 1] == "D" and PARENT_SEQ[OXY_IDX - 1] == "Q"
    OUT.mkdir(parents=True, exist_ok=True)

    # One run per open set. Singles first (guaranteed <=1 edit), then pairs
    # drawn from the four positions nearest the tip (<=2 edits).
    open_sets: list[tuple[str, list[int]]] = []
    for p in ROTAMER_SHELL:
        open_sets.append((f"r{p}", [p]))
    for a, b in itertools.combinations(ROTAMER_SHELL[:4], 2):
        open_sets.append((f"r{a}_{b}", [a, b]))
    for p in DIVERSITY_SHELL:
        open_sets.append((f"d{p}", [p]))

    print(f"parent d045: {n_res} aa, Asp{BASE_IDX}/Gln{OXY_IDX}")
    print(f"{len(open_sets)} open sets x {N_SEQS} seqs, T={TEMPERATURE}, "
          f"<= {MAX_EDITS} edits\n")

    seen: dict[str, dict] = {}
    for tag, designable in open_sets:
        wdir = OUT / "mpnn" / tag
        wdir.mkdir(parents=True, exist_ok=True)
        try:
            res = ligandmpnn(PARENT_PDB, wdir,
                             redesigned_residues=redesign_spec(designable),
                             omit_per_residue=omit_spec(designable),
                             n_seqs=N_SEQS, n_batches=N_BATCHES,
                             temperature=TEMPERATURE, seed=37)
        except RunnerError as exc:
            print(f"  {tag:10} FAILED: {str(exc)[:120]}")
            continue
        if res.returncode != 0:
            print(f"  {tag:10} FAILED rc={res.returncode}, see {res.log}")
            continue
        kept = 0
        for _hdr, s in res.sequences:
            if len(s) != n_res:
                continue
            # the gate that matters: never ship a mutated catalytic pair
            if s[BASE_IDX - 1] != "D" or s[OXY_IDX - 1] != "Q":
                continue
            e = edits(s)
            if not e or len(e) > MAX_EDITS:
                continue
            if s in seen:
                continue
            seen[s] = {"tag": tag, "open": designable,
                       "edits": [f"{w}{p}{m}" for p, w, m in e],
                       "n_edits": len(e), "source": "ligandmpnn"}
            kept += 1
        print(f"  {tag:10} open={designable}  -> {kept} new variant(s)")

    n_mpnn = len(seen)
    # Mechanistic proposals, added second so an identical sequence keeps its
    # LigandMPNN provenance (the model found it too - worth knowing).
    added = 0
    for s, meta in ladder_variants().items():
        if s not in seen:
            seen[s] = meta
            added += 1
    print(f"\n  steric ladder -> {added} further variant(s) "
          f"({len(STERIC_LADDER)} rungs x {len(ROTAMER_SHELL)} positions, "
          "minus natives and duplicates)")

    if not seen:
        print("\nNo variants produced. Widen the shell or raise the temperature.")
        return 3

    manifest = []
    for i, (s, meta) in enumerate(sorted(seen.items(), key=lambda kv: kv[1]["n_edits"])):
        did = f"ins{i:03d}_{meta['tag']}_{'_'.join(meta['edits'])}"
        manifest.append({"design_id": did, "sequence": s, "length": len(s),
                         "parent": "rfd4_design_045", **meta})
    (OUT / "variants.json").write_text(json.dumps(manifest, indent=2))
    with (OUT / "variants.fasta").open("w") as fh:
        for m in manifest:
            fh.write(f">{m['design_id']} edits={','.join(m['edits'])}\n"
                     f"{m['sequence']}\n")

    by_n: dict[int, int] = {}
    by_src: dict[str, int] = {}
    for m in manifest:
        by_n[m["n_edits"]] = by_n.get(m["n_edits"], 0) + 1
        by_src[m["source"]] = by_src.get(m["source"], 0) + 1
    print(f"\n{len(manifest)} unique variants  "
          + "  ".join(f"{k} edit(s): {v}" for k, v in sorted(by_n.items())))
    print("  by source: " + "  ".join(f"{k}: {v}" for k, v in sorted(by_src.items()))
          + f"   (LigandMPNN found {n_mpnn} before the ladder was added)")
    print(f"wrote {OUT / 'variants.json'} and variants.fasta")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
