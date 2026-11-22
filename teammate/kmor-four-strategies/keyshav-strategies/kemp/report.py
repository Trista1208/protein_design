"""Live reporting: metrics table, per-metric plots, staged PDB files.

Reads outputs/ledger.jsonl (which the sweep appends to after every design) and
regenerates the published artefacts. Safe to call repeatedly while a sweep runs.

Chart design notes, so the choices are auditable rather than taste:

* ONE PANEL PER METRIC, single series each. The seven official sub-scores live
  on different scales and mean different things, so putting them on shared axes
  would be a dual-axis chart - the single worst chart mistake. Small multiples
  instead. A single series also needs no legend (the panel title names it) and
  no categorical palette, which is why no categorical validator run is needed.
* Colour does one job per panel: magnitude. Single sequential hue (blue 400).
* Status colour is reserved and used only for eligibility, always with a text
  label beside it, never colour alone.
* A reference line marks the crystallographic control (HG3.17 = 0.959) on the
  official-score panel: the number to beat, and the only honest yardstick.
* Thin marks, recessive grid, no value printed on every bar.

Palette values are the validated reference instance from the dataviz skill.
"""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path

from . import config as C
from . import objective
from . import official

# ---------------------------------------------------------------- design tokens
LIGHT = {
    "surface": "#fcfcfb",
    "text_primary": "#0b0b0b",
    "text_secondary": "#52514e",
    "text_muted": "#898781",
    "grid": "#e1e0d9",
    "series": "#3987e5",       # sequential blue 400
    "series_soft": "#86b6ef",  # blue 250
    "good": "#0ca30c",
    "critical": "#d03b3b",
    "warning": "#fab219",
}
DARK = {
    "surface": "#1a1a19",
    "text_primary": "#f0efec",
    "text_secondary": "#c3c2b7",
    "text_muted": "#898781",
    "grid": "#2c2c2a",
    "series": "#3987e5",
    "series_soft": "#1c5cab",
    "good": "#0ca30c",
    "critical": "#d03b3b",
    "warning": "#fab219",
}

#: the seven official terms plus the scalar, in weight order
PANELS = [
    ("official_score", "Official score  S = prod s_i^w_i", None),
    ("base", "base  (w 0.25)  -  d_base", 0.25),
    ("ipde", "ipde  (w 0.20)  -  complex_ipde", 0.20),
    ("oxy", "oxy  (w 0.20)  -  d_oxy", 0.20),
    ("orient_base", "orient_base  (w 0.12)", 0.12),
    ("orient_oxy", "orient_oxy  (w 0.08)", 0.08),
    ("bury", "bury  (w 0.08)  -  ca_neighbors_mean", 0.08),
    ("cons", "cons  (w 0.07)  -  sd_base", 0.07),
]

#: crystallographic controls, measured by APNovo on the real structures
CONTROL_POSITIVE_SCORE = 0.9591   # 4BS0, HG3.17 with a real Gln50
CONTROL_NEGATIVE_SCORE = 0.0      # 6TU6, Q50M knockout - no oxyanion donor

PUBLISH = C.ROOT / "published"


# ---------------------------------------------------------------- data
def load_rows(ledger: Path | None = None) -> list[dict]:
    ledger = ledger or (C.OUTPUTS / "ledger.jsonl")
    if not ledger.exists():
        return []
    rows = []
    for line in ledger.read_text().splitlines():
        if not line.strip():
            continue
        try:
            rows.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    # de-duplicate by design_id, keeping the richest (latest) record
    best: dict[str, dict] = {}
    for r in rows:
        did = r.get("design_id")
        if not did:
            continue
        prev = best.get(did)
        if prev is None or len(json.dumps(r)) >= len(json.dumps(prev)):
            best[did] = r
    return [best[k] for k in sorted(best)]


def _sub(row: dict, key: str):
    if key == "official_score":
        return (row.get("official") or {}).get("official_score") or row.get("score")
    return ((row.get("official") or {}).get("subscores") or {}).get(key)


# ---------------------------------------------------------------- table
#: The dihedral columns are NOT optional. Six of eight designs died on a term
#: that no table showed (LESSONS 1.3, 7.3); `dihedral_oxy_deg` and its band gap
#: are that term. `msa_mode` is here because a score without its fold protocol
#: is not comparable to anything (LESSONS 3.2).
TABLE_COLS = ["design_id", "variant", "length", "stage", "eligible",
              "official_score", "dense_rank", "base", "oxy", "ipde",
              "orient_base", "orient_oxy", "bury", "cons", "limiting_term",
              "worst_term", "d_base_A", "d_oxy_A", "dihedral_base_deg",
              "dihedral_oxy_deg", "gap_dihedral_base_deg", "gap_dihedral_oxy_deg",
              "ca_neighbors_mean", "mean_plddt", "msa_mode", "wall_s"]


def build_table(rows: list[dict]) -> list[dict]:
    out = []
    for r in rows:
        off = r.get("official") or {}
        det = off.get("detail") or {}
        dense = r.get("dense") or {}
        out.append({
            "design_id": r.get("design_id"),
            "variant": r.get("variant"),
            "length": r.get("length"),
            "stage": r.get("stage"),
            "eligible": r.get("eligible"),
            "official_score": _sub(r, "official_score"),
            **{k: _sub(r, k) for k in ("base", "oxy", "ipde", "orient_base",
                                       "orient_oxy", "bury", "cons")},
            "limiting_term": off.get("limiting_term"),
            "d_base_A": (det.get("base") or {}).get("d_base_A"),
            "d_oxy_A": (det.get("oxy") or {}).get("d_oxy_A"),
            **{k: dense.get(k) for k in
               ("dihedral_base_deg", "dihedral_oxy_deg",
                "gap_dihedral_base_deg", "gap_dihedral_oxy_deg")},
            "ca_neighbors_mean": (det.get("bury") or {}).get("ca_neighbors_mean"),
            "mean_plddt": (r.get("foldability") or {}).get("mean_plddt"),
            "msa_mode": (r.get("protocol") or {}).get("msa_mode"),
            "wall_s": r.get("wall_s"),
        })
    add_dense_rank(rows, out)
    return out


def add_dense_rank(rows: list[dict], table: list[dict]) -> None:
    """Attach the internal percentile ranking, in place.

    Ranks only the rows that actually carry dense defects, so unscored and
    rejected designs do not dilute the batch percentiles. Designs without a
    measurement get ``None`` rather than a middling rank they did not earn.
    """
    idx = [i for i, r in enumerate(rows) if (r.get("dense") or {}).get("badness")]
    for row in table:
        row.setdefault("dense_rank", None)
        row.setdefault("worst_term", None)
    if not idx:
        return
    ranked = objective.rank_batch([rows[i]["dense"]["badness"] for i in idx])
    for i, rk in zip(idx, ranked):
        table[i]["dense_rank"] = rk["dense_rank"]
        table[i]["worst_term"] = rk["worst_term"]


def write_csv(table: list[dict], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TABLE_COLS)
        w.writeheader()
        for row in table:
            w.writerow({k: row.get(k) for k in TABLE_COLS})
    return path


def _fmt(v, nd=4):
    if v is None:
        return "-"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.{nd}f}"
    return str(v)


def write_markdown(table: list[dict], path: Path, *, extra: dict | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    scored = [r for r in table if r.get("official_score") is not None]
    nonzero = [r for r in scored if (r.get("official_score") or 0) > 0]
    scored.sort(key=lambda r: -(r.get("official_score") or 0))

    L = ["# Kemp eliminase designs - live results", ""]
    if extra:
        L += [f"_Updated {extra.get('updated')} - "
              f"{extra.get('n_total', 0)} designs seen, "
              f"{len(nonzero)} with a non-zero score._", ""]
    L += ["## Yardstick", "",
          "| structure | official score | note |",
          "|---|---|---|",
          f"| HG3.17 (4BS0, crystal) | **{CONTROL_POSITIVE_SCORE:.4f}** | "
          "real refined Gln50 - the number to beat |",
          f"| HG3.17-Q50M (6TU6, crystal) | **{CONTROL_NEGATIVE_SCORE:.4f}** | "
          "oxyanion knockout - no donor, so S = 0 |",
          "",
          "The score is a weighted **geometric** mean, so a single zero "
          "sub-score makes the whole design 0.000. `limiting_term` names the "
          "component holding each design back.", ""]

    L += ["## Designs", "",
          "| " + " | ".join(TABLE_COLS) + " |",
          "|" + "|".join(["---"] * len(TABLE_COLS)) + "|"]
    ordered = scored + [r for r in table if r.get("official_score") is None]
    for r in ordered:
        L.append("| " + " | ".join(_fmt(r.get(k)) for k in TABLE_COLS) + " |")
    L += ["", "## Plots", ""]
    for key, title, _w in PANELS:
        L.append(f"### {title}")
        L.append(f"![{key}](plots/{key}.png)")
        L.append("")
    path.write_text("\n".join(L) + "\n")
    return path


# ---------------------------------------------------------------- plots
def _style(ax, tok, *, title: str, ylabel: str):
    ax.set_facecolor(tok["surface"])
    ax.set_title(title, color=tok["text_primary"], fontsize=10, loc="left", pad=8)
    ax.set_ylabel(ylabel, color=tok["text_secondary"], fontsize=8)
    ax.grid(axis="y", color=tok["grid"], linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(tok["grid"])
    ax.tick_params(colors=tok["text_muted"], labelsize=7, length=0)


def make_plots(table: list[dict], out_dir: Path, *, dark: bool = False) -> list[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tok = DARK if dark else LIGHT
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    rows = [r for r in table if r.get("official_score") is not None] or table
    rows = sorted(rows, key=lambda r: -(r.get("official_score") or 0))[:40]
    labels = [str(r.get("design_id", ""))[:18] for r in rows]

    for key, title, _w in PANELS:
        vals = [r.get(key) for r in rows]
        have = [(l, v) for l, v in zip(labels, vals) if v is not None]
        fig, ax = plt.subplots(figsize=(max(6, 0.38 * max(len(have), 6)), 3.2),
                               facecolor=tok["surface"])
        if not have:
            ax.text(0.5, 0.5, "no measurements yet", ha="center", va="center",
                    color=tok["text_muted"], fontsize=9, transform=ax.transAxes)
            _style(ax, tok, title=title, ylabel="")
            ax.set_xticks([])
        else:
            xs = range(len(have))
            ys = [v for _, v in have]
            # single series, single sequential hue: colour encodes magnitude only
            ax.bar(xs, ys, width=0.62, color=tok["series"], zorder=3,
                   edgecolor=tok["surface"], linewidth=1.0)
            _style(ax, tok, title=title,
                   ylabel="score (0-1)" if key != "official_score" else "S")
            ax.set_xticks(list(xs))
            ax.set_xticklabels([l for l, _ in have], rotation=60,
                               ha="right", fontsize=6)
            if key == "official_score":
                ax.axhline(CONTROL_POSITIVE_SCORE, color=tok["good"],
                           linewidth=2, linestyle="--", zorder=4)
                ax.text(len(have) - 0.4, CONTROL_POSITIVE_SCORE,
                        f" HG3.17 crystal {CONTROL_POSITIVE_SCORE:.3f}",
                        color=tok["good"], fontsize=7, va="bottom", ha="right")
                ax.set_ylim(0, max(1.0, max(ys) * 1.15, CONTROL_POSITIVE_SCORE * 1.15))
            else:
                ax.set_ylim(0, 1.05)
        fig.tight_layout()
        p = out_dir / (f"{key}_dark.png" if dark else f"{key}.png")
        fig.savefig(p, dpi=150, facecolor=tok["surface"])
        plt.close(fig)
        written.append(p)
    return written


# ---------------------------------------------------------------- pdb staging
def stage_structures(rows: list[dict], out_dir: Path, *, limit: int = 60) -> list[Path]:
    """Copy each design's folded structure (and backbone) into published/."""
    out_dir.mkdir(parents=True, exist_ok=True)
    staged: list[Path] = []
    for r in rows[:limit]:
        did = r.get("design_id")
        if not did:
            continue
        bb = r.get("backbone")
        if bb and Path(bb).exists():
            dst = out_dir / f"{did}_backbone.pdb"
            if not dst.exists() or dst.stat().st_mtime < Path(bb).stat().st_mtime:
                shutil.copy2(bb, dst)
            staged.append(dst)
        # the Boltz-folded model, converted to PDB by the pipeline
        wdir = C.RUNS / "m4" / did
        for cand in list(wdir.rglob("folded.pdb"))[:1]:
            dst = out_dir / f"{did}_folded.pdb"
            if not dst.exists() or dst.stat().st_mtime < cand.stat().st_mtime:
                shutil.copy2(cand, dst)
            staged.append(dst)
        # and the sequence
        if r.get("sequence"):
            (out_dir / f"{did}.fasta").write_text(
                f">{did} variant={r.get('variant')} "
                f"score={_sub(r,'official_score')}\n{r['sequence']}\n")
    return staged


# ---------------------------------------------------------------- driver
def publish(*, ledger: Path | None = None, dark_too: bool = True) -> dict:
    import time
    rows = load_rows(ledger)
    table = build_table(rows)
    PUBLISH.mkdir(parents=True, exist_ok=True)

    meta = {"updated": time.strftime("%F %T"), "n_total": len(table)}
    csv_p = write_csv(table, PUBLISH / "designs.csv")
    md_p = write_markdown(table, PUBLISH / "README.md", extra=meta)
    plots = make_plots(table, PUBLISH / "plots")
    if dark_too:
        plots += make_plots(table, PUBLISH / "plots", dark=True)
    pdbs = stage_structures(rows, PUBLISH / "structures")

    # Degeneracy alarm over the SCORED rows (LESSONS 4.4: "Would have fired on
    # day 1 here"). It is published rather than raised, because publish() runs
    # live during a campaign and must not kill a running batch - the launch
    # gate in objective.assert_ok() is what refuses to spend GPU time.
    scored = [r for r in table if r.get("official_score") is not None
              or r.get("dense_rank") is not None]
    degeneracy = objective.audit(scored, columns=[
        c for c in TABLE_COLS if c not in
        ("design_id", "variant", "stage", "eligible", "limiting_term",
         "worst_term", "msa_mode", "wall_s", "length")])
    (PUBLISH / "degeneracy.json").write_text(json.dumps(degeneracy, indent=2))

    # Protocol comparability: a batch folded under two different MSA settings
    # cannot be ranked against itself (LESSONS 3.1-3.2).
    msa_modes = sorted({r.get("msa_mode") for r in scored if r.get("msa_mode")})

    (PUBLISH / "status.json").write_text(json.dumps({
        **meta,
        "n_scored": sum(1 for r in table if r.get("official_score") is not None),
        "n_nonzero": sum(1 for r in table
                         if (r.get("official_score") or 0) > 0),
        "n_eligible": sum(1 for r in table if r.get("eligible")),
        "best": max((r.get("official_score") or 0) for r in table) if table else None,
        "best_dense_rank": max((r.get("dense_rank") or 0) for r in table) if table else None,
        "control_positive": CONTROL_POSITIVE_SCORE,
        "control_negative": CONTROL_NEGATIVE_SCORE,
        # `degeneracy_passed` on its own would read as a clean bill of health
        # on a batch too small to judge. An unmeasured check is not a passed
        # check (LESSONS 1.5), so the skip is published beside it and any gate
        # must require judged AND passed.
        "degeneracy_passed": degeneracy["passed"],
        "degeneracy_judged": not degeneracy["skipped"],
        "degeneracy_note": degeneracy.get("reason", ""),
        "degenerate_columns": degeneracy["failures"],
        "msa_modes_present": msa_modes,
        "protocol_mixed": len(msa_modes) > 1,
    }, indent=2))
    return {"csv": str(csv_p), "markdown": str(md_p),
            "plots": len(plots), "structures": len(pdbs), **meta}
