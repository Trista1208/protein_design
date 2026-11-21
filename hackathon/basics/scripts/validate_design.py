#!/usr/bin/env python3
"""Hard-gate checker for the Kemp hackathon. Stdlib only.

Fails (exit 1) only on rules the organizers stated:
  length, alphabet, identity > 70% vs published reference FASTA.

Warns (does not fail) on missing D/Q, odd Cys, empty sequence.
TM-score is NOT computed here — that needs a structure vs PDB.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

AA = set("ACDEFGHIKLMNPQRSTVWY")
MAX_LEN = 260
MAX_ID = 0.70

HERE = Path(__file__).resolve().parent
DEFAULT_REFS = HERE.parent / "references" / "published_kemp.fasta"


def parse_fasta(text: str) -> list[tuple[str, str]]:
    recs: list[tuple[str, str]] = []
    name, chunks = None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(chunks).upper()))
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append("".join(c for c in line.upper() if c.isalpha()))
    if name is not None:
        recs.append((name, "".join(chunks).upper()))
    return recs


def load_fasta(path: Path) -> list[tuple[str, str]]:
    return parse_fasta(path.read_text())


def needleman_identity(a: str, b: str) -> float:
    """Global alignment identity = matches / alignment length (including gaps)."""
    if not a or not b:
        return 0.0
    return _identity_with_pointers(a, b, gap=-1, match=1, mismatch=0)


def _identity_with_pointers(a: str, b: str, gap: int, match: int, mismatch: int) -> float:
    n, m = len(a), len(b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    ptr = [[0] * (m + 1) for _ in range(n + 1)]  # 0 diag, 1 up, 2 left
    for i in range(1, n + 1):
        score[i][0] = i * gap
        ptr[i][0] = 1
    for j in range(1, m + 1):
        score[0][j] = j * gap
        ptr[0][j] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d = score[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch)
            u = score[i - 1][j] + gap
            l = score[i][j - 1] + gap
            if d >= u and d >= l:
                score[i][j], ptr[i][j] = d, 0
            elif u >= l:
                score[i][j], ptr[i][j] = u, 1
            else:
                score[i][j], ptr[i][j] = l, 2
    i, j = n, m
    matches = aln = 0
    while i > 0 or j > 0:
        p = ptr[i][j] if i and j else (1 if i else 2)
        if i > 0 and j > 0 and p == 0:
            if a[i - 1] == b[j - 1]:
                matches += 1
            aln += 1
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or p == 1):
            aln += 1
            i -= 1
        else:
            aln += 1
            j -= 1
    return matches / aln if aln else 0.0


def max_identity(seq: str, refs: list[tuple[str, str]]) -> tuple[float, str]:
    best, who = -1.0, ""
    for name, ref in refs:
        ident = needleman_identity(seq, ref)
        if ident > best:
            best, who = ident, name
    return best, who


def check_one(name: str, seq: str, refs: list[tuple[str, str]]) -> dict:
    errors: list[str] = []
    warns: list[str] = []
    bad = [c for c in seq if c not in AA]
    if not seq:
        errors.append("empty sequence")
    if bad:
        errors.append(f"non-canonical letters: {''.join(sorted(set(bad)))}")
    if len(seq) > MAX_LEN:
        errors.append(f"length {len(seq)} > {MAX_LEN}")
    if len(seq) < 40:
        warns.append(f"very short ({len(seq)} aa) — unusual for a pocket + scaffold")
    if "D" not in seq:
        warns.append("no Asp (D); slides use Asp as the catalytic base")
    if "Q" not in seq:
        warns.append("no Gln (Q); slides use Gln as the H-bond donor")
    n_c = seq.count("C")
    if n_c % 2 == 1:
        warns.append(f"odd number of Cys ({n_c}) — unpaired cysteine risk")
    ident, who = max_identity(seq, refs) if refs else (0.0, "")
    if ident > MAX_ID:
        errors.append(f"identity {ident:.1%} vs {who} exceeds {MAX_ID:.0%}")
    return {
        "id": name,
        "length": len(seq),
        "max_identity": round(ident, 4),
        "closest_published": who,
        "n_asp": seq.count("D"),
        "n_gln": seq.count("Q"),
        "n_cys": n_c,
        "errors": errors,
        "warnings": warns,
        "pass": not errors,
    }


def main() -> int:
    p = argparse.ArgumentParser(description="Kemp hackathon hard gates (length + identity).")
    p.add_argument("fasta", nargs="?", help="design FASTA to check")
    p.add_argument("--refs", type=Path, default=DEFAULT_REFS)
    p.add_argument("--seq", help="single sequence instead of a FASTA file")
    p.add_argument("--name", default="design")
    p.add_argument("--json", action="store_true")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()

    refs = load_fasta(args.refs) if args.refs.exists() else []
    if not refs:
        print(f"warning: no references at {args.refs}", file=sys.stderr)

    if args.self_test:
        return _self_test(refs)

    if args.seq:
        designs = [(args.name, args.seq.upper())]
    elif args.fasta:
        designs = load_fasta(Path(args.fasta))
    else:
        p.error("provide a FASTA path, --seq, or --self-test")

    reports = [check_one(n, s, refs) for n, s in designs]
    if args.json:
        json.dump(reports, sys.stdout, indent=2)
        print()
    else:
        for r in reports:
            flag = "PASS" if r["pass"] else "FAIL"
            print(f"{flag}  {r['id']}  len={r['length']}  ident={r['max_identity']:.1%} vs {r['closest_published']}")
            for e in r["errors"]:
                print(f"  ERROR  {e}")
            for w in r["warnings"]:
                print(f"  WARN   {w}")
    return 0 if all(r["pass"] for r in reports) else 1


def _self_test(refs: list[tuple[str, str]]) -> int:
    assert refs, "need published_kemp.fasta"
    clone = refs[0][1]
    r_fail = check_one("clone_of_published", clone, refs)
    r_ok = check_one(
        "toy_new",
        "MEEKLKEKLEKALKEKLKEALKEKLKEALKDKLKEKLQEKALKEKLKEALKEK",
        refs,
    )
    print("self-test clone (must FAIL identity):", r_fail["pass"], r_fail["errors"])
    print("self-test toy  (must PASS identity):", r_ok["pass"], f"ident={r_ok['max_identity']:.1%}")
    ok = (not r_fail["pass"]) and r_ok["pass"] and r_ok["length"] <= MAX_LEN
    print("self-test", "OK" if ok else "BROKEN")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
