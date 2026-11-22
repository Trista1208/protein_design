"""Novelty gates wired to real tools (mmseqs2 + foldseek), in the `novelty` env.

Handover 4.3.8: screen against a real corpus and record BOTH TM normalisation
directions. 4.3.9: fail closed.

foldseek gives us qtmscore and ttmscore in one search, which is exactly the
"both directions" the handover asks for. We take the MAX as the conservative
answer, because either direction exceeding 0.7 is a disqualification risk.
"""

from __future__ import annotations

import csv
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from . import config as C

ENVS = Path("/home/mandalorian/anaconda3/envs")
MMSEQS = ENVS / "novelty/bin/mmseqs"
FOLDSEEK = ENVS / "novelty/bin/foldseek"
CORPUS = C.ROOT / "inputs" / "corpus"
SEQDB = CORPUS / "seqdb"
STRUCTDB = CORPUS / "structdb"
PUBLISHED_FASTA = CORPUS / "published.fasta"


class NoveltyError(RuntimeError):
    pass


@dataclass
class NoveltyHit:
    target: str
    value: float
    kind: str          # "identity" | "tm"
    detail: str = ""


def tools_available() -> dict[str, bool]:
    return {
        "mmseqs": MMSEQS.exists(),
        "foldseek": FOLDSEEK.exists(),
        "seqdb": SEQDB.exists() or (CORPUS / "seqdb.dbtype").exists(),
        "structdb": STRUCTDB.exists() or (CORPUS / "structdb.dbtype").exists(),
        "published_fasta": PUBLISHED_FASTA.exists(),
    }


# --------------------------------------------------------------------------- identity
def max_identity(sequence: str, *, tmp: Path | None = None) -> tuple[float, str]:
    """Highest whole-chain identity against the published corpus.

    Returns (identity_fraction, target). Raises if the toolchain is missing -
    the caller must treat that as a FAILED gate, not a pass.
    """
    avail = tools_available()
    if not avail["mmseqs"]:
        raise NoveltyError("mmseqs not installed; cannot prove sequence novelty")
    if not avail["published_fasta"]:
        raise NoveltyError(f"no corpus at {PUBLISHED_FASTA}")

    work = Path(tmp or tempfile.mkdtemp(prefix="kemp_ident_"))
    work.mkdir(parents=True, exist_ok=True)
    qf = work / "query.fasta"
    qf.write_text(f">query\n{sequence}\n")
    out = work / "hits.tsv"

    cmd = [str(MMSEQS), "easy-search", str(qf), str(PUBLISHED_FASTA), str(out),
           str(work / "tmp"), "--format-output", "query,target,fident,alnlen,evalue",
           "-s", "7.5", "--max-seqs", "1000", "-e", "1000"]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if proc.returncode != 0:
        raise NoveltyError(f"mmseqs easy-search failed: {proc.stderr[-400:]}")

    best, target = 0.0, ""
    if out.exists():
        for row in csv.reader(out.read_text().splitlines(), delimiter="\t"):
            if len(row) < 3:
                continue
            try:
                fid = float(row[2])
            except ValueError:
                continue
            if fid > best:
                best, target = fid, row[1]
    return best, target


# --------------------------------------------------------------------------- TM
def max_tm(structure: Path | str, *, tmp: Path | None = None) -> tuple[float, str, dict]:
    """Highest TM-score against the published corpus, worst of both normalisations.

    foldseek reports qtmscore (normalised by query) and ttmscore (by target);
    the handover asks for both to be recorded.
    """
    avail = tools_available()
    if not avail["foldseek"]:
        raise NoveltyError("foldseek not installed; cannot prove structural novelty")
    if not avail["structdb"]:
        raise NoveltyError(f"no structure corpus at {STRUCTDB}")

    structure = Path(structure)
    if not structure.exists():
        raise NoveltyError(f"query structure missing: {structure}")

    work = Path(tmp or tempfile.mkdtemp(prefix="kemp_tm_"))
    work.mkdir(parents=True, exist_ok=True)
    out = work / "tm.tsv"
    cmd = [str(FOLDSEEK), "easy-search", str(structure), str(STRUCTDB), str(out),
           str(work / "tmp"), "--format-output",
           "query,target,alntmscore,qtmscore,ttmscore,fident,evalue",
           "--exhaustive-search", "1", "-e", "inf"]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if proc.returncode != 0:
        raise NoveltyError(f"foldseek easy-search failed: {proc.stderr[-400:]}")

    best, target = 0.0, ""
    rows: list[dict] = []
    if out.exists():
        for row in csv.reader(out.read_text().splitlines(), delimiter="\t"):
            if len(row) < 5:
                continue
            try:
                q, t = float(row[3]), float(row[4])
            except ValueError:
                continue
            rows.append({"target": row[1], "qtmscore": q, "ttmscore": t})
            worst = max(q, t)          # conservative
            if worst > best:
                best, target = worst, row[1]
    return best, target, {"n_hits": len(rows), "hits": rows[:10],
                          "note": "value is max(qtmscore, ttmscore) - the "
                                  "conservative direction (handover 4.3.8)"}


# --------------------------------------------------------------------------- combined
def screen(sequence: str, structure: Path | str | None = None) -> dict:
    """Run both novelty gates. Any failure or missing tool => not eligible."""
    result: dict = {"eligible": True, "reasons": []}

    try:
        ident, itarget = max_identity(sequence)
        result["identity"] = {"max": round(ident, 4), "target": itarget,
                              "threshold": C.SEQ_IDENTITY_DISQUALIFY}
        if ident > C.SEQ_IDENTITY_DISQUALIFY:
            result["eligible"] = False
            result["reasons"].append(
                f"sequence identity {ident:.1%} to {itarget} exceeds "
                f"{C.SEQ_IDENTITY_DISQUALIFY:.0%}")
    except NoveltyError as exc:
        result["eligible"] = False
        result["identity"] = {"error": str(exc)}
        result["reasons"].append(f"identity gate unprovable: {exc}")

    if structure is not None:
        try:
            tm, ttarget, detail = max_tm(structure)
            result["tm"] = {"max": round(tm, 4), "target": ttarget,
                            "threshold": C.TM_SCORE_DISQUALIFY, **detail}
            if tm > C.TM_SCORE_DISQUALIFY:
                result["eligible"] = False
                result["reasons"].append(
                    f"TM-score {tm:.3f} to {ttarget} exceeds {C.TM_SCORE_DISQUALIFY}")
        except NoveltyError as exc:
            result["eligible"] = False
            result["tm"] = {"error": str(exc)}
            result["reasons"].append(f"TM gate unprovable: {exc}")

    result["verdict"] = "NOVEL" if result["eligible"] else "DISQUALIFIED"
    return result
