"""Build the targeted design library around the Tyr103 mechanism.

    python scripts/build_library.py [--max 600]

NOT a blind sweep. Every block below tests something specific that the
measurements support.

THE MECHANISM (LOOP_DESIGN 6.5)
-------------------------------
`ins080` scores 0.743 on 1 of 5 seeds. The discriminator is a single
hydrogen bond:

    scoring seed      Tyr103.OH ... Gln46.NE2   2.36 A
    non-scoring seeds no contact under 5 A

The Tyr103 hydroxyl pins the Gln46 rotamer in the only geometry that scores.
`ins079` is the same design with Phe at 103 - one atom different, no hydroxyl
- and scores 0/5. Tyr103's only persistent contacts are Met45 and Leu49.

THE BLOCKS
----------
  A  DONOR SCAN at 103. Falsifiable prediction: residues that can donate or
     accept an H-bond and reach Gln46 (Y, W, T, S, N, Q, H, E, D, K, R)
     should score; those that cannot (F, L, I, V, A, M, G) should not. If
     Phe and Leu score as often as Tyr, the mechanism is wrong and this
     whole direction should be abandoned - which is why the non-donors are
     included rather than assumed away.

  B  ANCHOR SCAN at Met45 and Leu49 - the two residues always in contact
     with Tyr103. Changing their bulk biases the Tyr103 chi angles, which is
     the lever on whether the H-bond forms at all.

  C  SEQUESTER FIX at Arg127. Literature hard filter H3: both crystal
     controls have zero charged residues within 4 A of the catalytic base;
     every design of ours has Arg127 at 3.3-4.0 A. KE07 gained 2.5-fold
     kcat/KM from removing the equivalent bridge. R127A already produced our
     first LIT-PASS.

  D  DOUBLES: anchor x donor, and R127 x donor. Combinations are the only
     thing that has produced a gain here - F42Y and H103Y each failed alone
     and together gave the first design with six of seven terms non-zero.
     Bhattacharya 2022 reports the same effect (triple 3x better than
     additivity predicted).

  E  HYDROPHOBIC SHELL. Privett S2 wants >=0.60 apolar in the 5-6 A shell
     around the base; the controls sit at 0.727 and ins080 at 0.500.

Catalytic Asp102 and Gln46 are never touched.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/library"
BASE_IDX, OXY_IDX = 102, 46
PARENT_ID = "ins080_r103_42_F42Y_H103Y"

#: A - can it donate/accept an H-bond to the Gln46 amide, and reach it?
DONORS = list("YWTSNQHEDKR")
NON_DONORS = list("FLIVAMG")          # the control arm of the prediction

#: B - Tyr103's persistent anchors (the only residues in contact in all
#: five ins080 structures)
ANCHORS = {45: list("AVLIFWMTGSY"), 49: list("AVLIFWMTGSY")}

#: C - the salt-bridging arginine
SEQUESTER = {127: list("AVLIMFQNSTGY")}

#: E - shell positions, apolar substitutions to raise the S2 fraction
SHELL = {26: list("AVLIFWM"), 42: list("AVLIFWMY"), 84: list("AVLIFWM"),
         99: list("AVLIFWM"), 23: list("AVLIFWM"), 27: list("AVLIFWM"),
         77: list("AVLIFWMY"), 88: list("AVLIFWM")}


def parent_sequence() -> str:
    man = json.loads((ROOT / "outputs/insurance/variants.json").read_text())
    for v in man:
        if v["design_id"] == PARENT_ID:
            return v["sequence"]
    raise SystemExit(f"{PARENT_ID} not found")


def mutate(seq: str, subs) -> str:
    s = list(seq)
    for pos, aa in subs:
        s[pos - 1] = aa
    return "".join(s)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max", type=int, default=600)
    a = ap.parse_args()
    seq0 = parent_sequence()
    assert seq0[BASE_IDX - 1] == "D" and seq0[OXY_IDX - 1] == "Q"
    OUT.mkdir(parents=True, exist_ok=True)

    lib: dict[str, dict] = {}

    def add(subs, block, note=""):
        s = mutate(seq0, subs)
        if s == seq0 or s in lib:
            return
        if s[BASE_IDX - 1] != "D" or s[OXY_IDX - 1] != "Q":
            return
        tag = "_".join(f"{seq0[p-1]}{p}{aa}" for p, aa in subs)
        lib[s] = {"design_id": f"L_{block}_{tag}", "block": block,
                  "mutations": [f"{seq0[p-1]}{p}{aa}" for p, aa in subs],
                  "note": note}

    # A - donor scan at 103, with the non-donor control arm
    for aa in DONORS:
        add([(103, aa)], "A_donor", "predicted to score: can H-bond Gln46")
    for aa in NON_DONORS:
        add([(103, aa)], "A_control", "predicted NOT to score: no H-bond")

    # B - anchors
    for pos, aas in ANCHORS.items():
        for aa in aas:
            add([(pos, aa)], "B_anchor", "biases the Tyr103 chi angles")

    # C - sequestering Arg
    for pos, aas in SEQUESTER.items():
        for aa in aas:
            add([(pos, aa)], "C_sequester", "removes the Arg127-Asp102 bridge")

    # E - hydrophobic shell
    for pos, aas in SHELL.items():
        for aa in aas:
            add([(pos, aa)], "E_shell", "raises the apolar fraction (S2)")

    # D - doubles: anchor x donor, sequester x donor, anchor x sequester
    best_donors = list("YWTNQ")
    for pos, aas in ANCHORS.items():
        for aa, d in itertools.product("AVLIFWM", best_donors):
            add([(pos, aa), (103, d)], "D_anchor_donor",
                "lock the Tyr103 chi AND keep a donor at 103")
    for aa, d in itertools.product(SEQUESTER[127], best_donors):
        add([(127, aa), (103, d)], "D_seq_donor",
            "remove the salt bridge AND keep a donor")
    for aa, s in itertools.product("AVLIM", "AVLIFWM"):
        add([(127, aa), (45, s)], "D_seq_anchor",
            "remove the salt bridge AND bias the anchor")

    manifest = [{"sequence": s, "parent": PARENT_ID, **m}
                for s, m in lib.items()][:a.max]
    (OUT / "library.json").write_text(json.dumps(manifest, indent=2))
    with (OUT / "library.fasta").open("w") as fh:
        for m in manifest:
            fh.write(f">{m['design_id']} {'+'.join(m['mutations'])}\n"
                     f"{m['sequence']}\n")

    from collections import Counter
    c = Counter(m["block"] for m in manifest)
    print(f"{len(manifest)} designs (parent {PARENT_ID})")
    for k in sorted(c):
        print(f"   {k:18} {c[k]:4}")
    print(f"\nwrote {OUT/'library.json'}")
    print("\nblock A is a falsifiable test: if A_control (F/L/I/V/A/M/G at")
    print("103) scores as often as A_donor, the Tyr103 H-bond mechanism is")
    print("wrong and this direction should be abandoned.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
