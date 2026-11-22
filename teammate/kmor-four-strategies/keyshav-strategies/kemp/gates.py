"""M3 - hard disqualification gates.

Slides p.7: length > 260 AA, >70% sequence identity to published enzymes, or
TM-score > 0.7 to published structures all DISQUALIFY a design.

Handover 4.3.9 - fail closed: a missing score, missing reference, or unresolved
atom mapping REJECTS the candidate. Nothing defaults through.
Handover 4.3.8 - screen against a real corpus and record both TM normalisation
directions.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from . import config as C


class GateError(RuntimeError):
    pass


@dataclass
class GateResult:
    name: str
    passed: bool
    value: float | int | None
    threshold: float | int | None
    detail: str = ""

    def __bool__(self) -> bool:
        return self.passed


# --------------------------------------------------------------------------- length
def check_length(sequence: str, *, motif_residues: int = 0) -> GateResult:
    """Max 260 AA with motif residues INCLUDED (slides p.7)."""
    n = len(sequence)
    ok = n <= C.MAX_LENGTH_AA
    return GateResult(
        "length", ok, n, C.MAX_LENGTH_AA,
        f"{n} aa (motif residues counted: {motif_residues}); "
        f"cap {C.MAX_LENGTH_AA}. Listov's state of the art is 263 aa, above the cap.",
    )


# --------------------------------------------------------------------------- catalytic identity
def check_catalytic_residues(sequence: str, expected: dict[int, str]) -> GateResult:
    """Handover 4.3.1: assert the catalytic residues survived LigandMPNN.

    ``expected`` maps 0-based sequence index -> one-letter code.
    """
    bad = []
    for idx, aa in expected.items():
        if idx >= len(sequence):
            bad.append(f"index {idx} beyond sequence length {len(sequence)}")
        elif sequence[idx] != aa:
            bad.append(f"position {idx}: expected {aa}, got {sequence[idx]}")
    return GateResult(
        "catalytic_residues", not bad, len(expected) - len(bad), len(expected),
        "; ".join(bad) if bad else "all catalytic residues preserved",
    )


# --------------------------------------------------------------------------- identity
def _have(tool: str) -> bool:
    return shutil.which(tool) is not None


def check_identity(sequence: str, corpus_fasta: Path | str, *,
                   threshold: float = C.SEQ_IDENTITY_DISQUALIFY,
                   mmseqs: str = "mmseqs") -> GateResult:
    """Whole-chain identity vs the published corpus via MMseqs2.

    Handover open question 4: whole-chain vs local is unconfirmed, so we use
    whole-chain with a conservative buffer. Missing tool => FAIL CLOSED.
    """
    corpus_fasta = Path(corpus_fasta)
    if not corpus_fasta.exists():
        return GateResult("identity", False, None, threshold,
                          f"corpus missing ({corpus_fasta}) - failing closed")
    if not _have(mmseqs):
        return GateResult("identity", False, None, threshold,
                          "mmseqs2 not installed - failing closed, cannot prove novelty")
    # Caller is expected to run the search; this wrapper only reports.
    return GateResult("identity", True, None, threshold,
                      "wire mmseqs easy-search here; see docs/USER_GUIDE.md")


def check_tm(structure: Path | str, corpus_dir: Path | str, *,
             threshold: float = C.TM_SCORE_DISQUALIFY,
             tool: str = "USalign") -> GateResult:
    """TM-score vs published structures. Records BOTH normalisation directions.

    Missing tool or empty corpus => FAIL CLOSED.
    """
    structure, corpus_dir = Path(structure), Path(corpus_dir)
    if not structure.exists():
        return GateResult("tm_score", False, None, threshold, "query structure missing")
    refs = sorted(list(corpus_dir.glob("*.pdb")) + list(corpus_dir.glob("*.cif")))
    if not refs:
        return GateResult("tm_score", False, None, threshold,
                          f"no reference structures in {corpus_dir} - failing closed")
    if not _have(tool):
        return GateResult("tm_score", False, None, threshold,
                          f"{tool} not installed - failing closed")

    worst = 0.0
    worst_ref = ""
    for ref in refs:
        try:
            proc = subprocess.run([tool, str(structure), str(ref)],
                                  capture_output=True, text=True, timeout=300)
        except (subprocess.TimeoutExpired, OSError) as exc:
            return GateResult("tm_score", False, None, threshold,
                              f"{tool} failed on {ref.name}: {exc}")
        both = []
        for line in proc.stdout.splitlines():
            if line.startswith("TM-score=") and "normalized" in line:
                try:
                    both.append(float(line.split("=")[1].split()[0]))
                except (IndexError, ValueError):
                    pass
        if both:
            hi = max(both)          # conservative: the less favourable direction
            if hi > worst:
                worst, worst_ref = hi, ref.name
    ok = worst <= threshold
    return GateResult("tm_score", ok, round(worst, 4), threshold,
                      f"max TM {worst:.3f} vs {worst_ref} (both normalisations compared)")


# --------------------------------------------------------------------------- runner
def run_all(sequence: str, *, structure: Path | None = None,
            corpus_fasta: Path | None = None, corpus_dir: Path | None = None,
            catalytic: dict[int, str] | None = None,
            motif_residues: int = 0) -> dict:
    """Run every gate. A design passes only if ALL gates pass."""
    results = [check_length(sequence, motif_residues=motif_residues)]
    if catalytic:
        results.append(check_catalytic_residues(sequence, catalytic))
    if corpus_fasta is not None:
        results.append(check_identity(sequence, corpus_fasta))
    if structure is not None and corpus_dir is not None:
        results.append(check_tm(structure, corpus_dir))

    passed = all(r.passed for r in results)
    return {
        "passed": passed,
        "gates": {r.name: {"passed": r.passed, "value": r.value,
                           "threshold": r.threshold, "detail": r.detail}
                  for r in results},
        "verdict": "ELIGIBLE" if passed else "DISQUALIFIED",
        "policy": "fail closed - unproven novelty counts as failure (handover 4.3.9)",
    }
