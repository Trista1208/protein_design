"""Library 4: cross the proven singles with the shell positions never paired.

WHY THIS AND NOT MORE OF THE SAME
---------------------------------
SUBMISSION.md allows 20 designs, not 5. That changes what is worth running.
Ranks 16-20 of the current submission score 0 on three of five seeds
(reproducibility 2/5, median 0.0000). Any new design that reaches 3/5
displaces one of them, so the marginal slot is cheap to beat.

The 23 x 45 x 84 x 88 ladder is exhausted - pairs lifted the in-band rate
from 10% to 47%, triples added nothing (LOOP_DESIGN 6.11), so stacking more
of those positions is not where the remaining gain is.

What was never tried: positions 26, 27, 49, 77 and 99 were generated as
SINGLES in library 1 and never crossed with a winning single. Library 2 only
paired winners with winners. This block fills that gap.

Catalytic Asp102 and Gln46 are never touched, and neither are the two pinner
positions 42/103 - the Tyr42/Tyr103 hydroxyls are the mechanism that makes
the shell score at all (LOOP_DESIGN 6.5), so changing them here would confound
the test.

Cysteine stays excluded: a lone free Cys cannot form a disulfide and a second
one could scramble in a reducing cytoplasm (LOOP_DESIGN 6.12).
"""
from __future__ import annotations
import argparse, itertools, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/library4"
PARENT_ID = "ins080_r103_42_F42Y_H103Y"

#: singles that reached >=2/5 reproducibility on their own or inside a
#: scoring pair - the proven half of each new combination
WINNERS = [(23, "L"), (23, "V"), (84, "I"), (84, "L"), (84, "M"),
           (45, "Y"), (88, "F")]

#: shell positions generated as singles in library 1 but NEVER paired
UNPAIRED = {26: list("AVLIFWM"), 27: list("AVLIFWM"),
            49: list("AVLIFWMT"), 77: list("AVLIFWMY"),
            99: list("AVLIFWM")}

FROZEN = {102, 46, 42, 103}
OMIT = "C"


def parent_sequence() -> str:
    v = json.loads((ROOT / "outputs/finalists" / PARENT_ID / "verdict.json").read_text())
    return v["sequence"]


def mutate(seq: str, subs) -> str:
    s = list(seq)
    for pos, aa in subs:
        s[pos - 1] = aa
    return "".join(s)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max", type=int, default=260)
    a = ap.parse_args()
    parent = parent_sequence()
    seen, out = {parent}, []

    for (wp, wa), (up, ual) in itertools.product(WINNERS, UNPAIRED.items()):
        if wp in FROZEN or up in FROZEN or wp == up:
            continue
        for ua in ual:
            if ua in OMIT or parent[up - 1] == ua:
                continue
            subs = [(wp, wa), (up, ua)]
            seq = mutate(parent, subs)
            if seq in seen:
                continue
            seen.add(seq)
            wt_w, wt_u = parent[wp - 1], parent[up - 1]
            out.append({
                "design_id": f"L4_{wt_w}{wp}{wa}_{wt_u}{up}{ua}",
                "sequence": seq, "parent": PARENT_ID,
                "mutations": [f"{wt_w}{wp}{wa}", f"{wt_u}{up}{ua}"],
                "block": "winner x unpaired-shell"})
            if len(out) >= a.max:
                break
        if len(out) >= a.max:
            break

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "library.json").write_text(json.dumps(out, indent=2))
    pos = sorted({int(m[1:-1]) for d in out for m in d["mutations"]})
    print(f"{len(out)} designs, all 2 mutations from {PARENT_ID}")
    print(f"positions touched: {pos}")
    print(f"-> {OUT/'library.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
