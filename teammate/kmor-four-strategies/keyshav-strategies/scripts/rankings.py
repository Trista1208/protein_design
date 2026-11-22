"""Rank every scoring design by every metric, not just the official score.

    python scripts/rankings.py

The official score decides the submission. But a design can be first on score
and middling on stability, and the other way round, and that is worth seeing
rather than hiding behind one number. This ranks all scoring designs
separately under each metric and shows where each one sits.

Metrics whose direction differs are labelled, so "rank 1" always means best:
  official median / best / pLDDT / core burial      higher is better
  iPDE / relax-per-residue / ThermoMPNN ddG /
  MD RMSD / MD drift / instability index / TM-score  lower is better

A design missing a metric is listed as "-" and excluded from THAT ranking
only - never imputed, never dropped from the others.
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: (column, label, higher_is_better, decimals)
METRICS = [
    ("official_median", "Official score (median)", True, 4),
    ("official_best", "Official score (best seed)", True, 4),
    ("plddt", "Boltz-2 pLDDT", True, 3),
    ("ipde", "Interface error iPDE", False, 2),
    ("relax_per_res", "Rosetta relax REU/residue", False, 3),
    ("thermompnn_ddg_sum", "ThermoMPNN ddG kcal/mol", False, 2),
    ("md_rmsd_A", "MD RMSD at 350 K (A)", False, 2),
    ("md_drift_A_per_ns", "MD drift (A/ns)", False, 2),
    ("core_frac", "Hydrophobic core burial", True, 3),
    ("instability_index", "Instability index", False, 1),
    ("tm_max", "TM-score to published (novelty)", False, 3),
]


def main() -> int:
    f = ROOT / "published/biophysics.json"
    if not f.exists():
        print("no biophysics.json; run biophysics_table.py first")
        return 1
    rows = json.loads(f.read_text())
    designs = [r for r in rows if not r["design_id"].startswith("REF_")]
    ref = next((r for r in rows if r["design_id"].startswith("REF_")), None)

    L = ["# Rankings — every design under every metric", "",
         f"_{len(designs)} designs. **Rank 1 is always best**, whichever "
         f"direction the metric runs._", "",
         "The official score decides the submission. The rest are shown "
         "because a design can lead on score and sit mid-table on stability, "
         "or the reverse, and that is worth seeing rather than averaging "
         "away. **None of the stability metrics was used to rank the "
         "submission.**", ""]
    if ref:
        L += [f"The HG3.17 crystal (`{ref['design_id']}`) is shown as the "
              f"control arm wherever it has a value. It is not ranked - it is "
              f"the yardstick, not an entry.", ""]

    # per-metric leaderboards
    place: dict[str, dict[str, int]] = {}
    for col, label, hi, dp in METRICS:
        have = [r for r in designs if r.get(col) is not None]
        if not have:
            continue
        have.sort(key=lambda r: r[col], reverse=hi)
        for i, r in enumerate(have, 1):
            place.setdefault(r["design_id"], {})[col] = i
        arrow = "higher is better" if hi else "lower is better"
        L += [f"## {label}", "",
              f"_{arrow}; {len(have)} of {len(designs)} designs have this "
              f"metric._", "",
              "| rank | design | value | repro | official |",
              "|---|---|---|---|---|"]
        for i, r in enumerate(have, 1):
            med = r.get("official_median")
            med_s = "–" if med is None else f"{med:.4f}"
            L.append(f"| {i} | `{r['design_id']}` | **{r[col]:.{dp}f}** | "
                     f"{r.get('repro','–')} | {med_s} |")
        if ref and ref.get(col) is not None:
            L.append(f"| – | **`{ref['design_id']}`** (control) | "
                     f"*{ref[col]:.{dp}f}* | – | – |")
        L.append("")

    # the cross-metric summary
    #: short, UNAMBIGUOUS column headers - "Official score (median)" and
    #: "(best seed)" both abbreviate to "Score" otherwise
    SHORT = {
        "official_median": "score&nbsp;med",
        "official_best": "score&nbsp;best",
        "plddt": "pLDDT",
        "ipde": "iPDE",
        "relax_per_res": "relax",
        "thermompnn_ddg_sum": "ddG",
        "md_rmsd_A": "MD&nbsp;RMSD",
        "md_drift_A_per_ns": "MD&nbsp;drift",
        "core_frac": "core",
        "instability_index": "instab",
        "tm_max": "TM",
    }
    cols = [(c, l) for c, l, _, _ in METRICS if any(c in v for v in place.values())]
    L += ["## Where each design places, across every metric", "",
          "Each cell is that design's **rank** under that metric (1 = best). "
          "A dash means the metric has not been measured for it. Sorted by "
          "official median score.", "",
          "| design | " + " | ".join(SHORT.get(c, c) for c, _ in cols) + " |",
          "|---" * (len(cols) + 1) + "|"]
    order = sorted(designs,
                   key=lambda r: (-(r.get("official_median") or 0),
                                  -(r.get("official_best") or 0)))
    for r in order:
        p = place.get(r["design_id"], {})
        cells = " | ".join(str(p.get(c, "–")) for c, _ in cols)
        L.append(f"| `{r['design_id']}` | {cells} |")
    L.append("")

    (ROOT / "published/RANKINGS.md").write_text("\n".join(L))
    (ROOT / "published/rankings.json").write_text(json.dumps(place, indent=2))
    print(f"wrote published/RANKINGS.md ({len(designs)} designs, "
          f"{len(cols)} metrics)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
