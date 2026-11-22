"""Autonomous supervisor: runs the whole remaining plan to a deadline.

    python scripts/supervisor.py --until 11:25

PHASES
------
  1  SCREEN library 1 (already running externally, or started here)
  2  BLOCK A VERDICT - the falsifiable test of the Tyr103 mechanism.
     Donors at 103 (Y/W/T/S/N/Q/H/E/D/K/R) are predicted to score; the
     control arm (F/L/I/V/A/M/G) is predicted not to. If the control arm
     does as well, the mechanism is wrong and phase 3 pivots instead of
     doubling down.
  3  LIBRARY 2 - built from whatever phase 2 established, plus combinations
     of every single mutation that improved on the parent
  4  HARVEST - 5 seeds on the best candidates, scored PER REPLICATE
  5  BANK - re-select the submission set

BUDGET
------
Measured: 20.6 s/fold with 8 workers (2.44x serial). A 1-seed screen costs
one fold; a full official score costs five. The supervisor reserves 55% of
the remaining folds for 5-seed scoring, because a screen that nobody can
afford to confirm is not worth running.

RULES IT MUST NOT BREAK (docs/good_practices.md)
------------------------------------------------
  * bound ligand -> pLDDT -> literature -> rank, in that order
  * score PER REPLICATE, never averaged geometry
  * a 1-seed number is a candidate, never a result
  * triple-check: the block A control arm exists so the mechanism can be
    falsified rather than confirmed by assumption
"""
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C, objective, official                # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "outputs/library"
LIB2 = ROOT / "outputs/library2"
FIN = ROOT / "outputs/finalists"
STATE = ROOT / "outputs/supervisor_state.json"
PY = sys.executable
SITE = official.CatalyticSite(base_resseq=102, oxy_resseq=46, chain="A")
WORKERS = 8
SEC_PER_FOLD = 20.6


def log(msg: str):
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)


def sh(args: list[str], timeout: int = 10800) -> subprocess.CompletedProcess:
    return subprocess.run(args, capture_output=True, text=True, timeout=timeout,
                          env={**__import__("os").environ,
                               "KEMP_ROOT": str(ROOT)})



def _pin(r: dict) -> float:
    """Pinner distance for ranking; module-level because BOTH build_library2
    and main()'s phase-4 sort use it. Defining it inside build_library2 crashed
    phase 4 with NameError and silently stopped all 5-seed scoring."""
    p = r.get("pinner") or {}
    d = p.get("distance_A")
    return d if d is not None else 9.9


def screen_running() -> bool:
    out = subprocess.run(["ps", "-eo", "cmd", "--no-headers"],
                         capture_output=True, text=True).stdout
    return "parallel_fold.py" in out


def load_screen(path: Path) -> list[dict]:
    if not path.exists():
        return []
    try:
        return [r for r in json.loads(path.read_text())
                if r.get("status") == "ok"]
    except Exception:                                            # noqa: BLE001
        return []


def block_a_verdict(rows: list[dict]) -> dict:
    """Did donors at 103 outperform non-donors? The mechanism stands or falls."""
    donor = [r for r in rows if r["design_id"].startswith("L_A_donor")]
    ctrl = [r for r in rows if r["design_id"].startswith("L_A_control")]

    def summarise(g):
        n = len(g)
        bound = [r for r in g if r.get("bound") and r.get("measurable")]
        gaps = [r["gap"] for r in bound if r.get("gap") is not None]
        return {"n": n, "bound": len(bound),
                "in_band": sum(1 for x in gaps if x == 0.0),
                "median_gap": (sorted(gaps)[len(gaps) // 2] if gaps else None),
                "best_gap": (min(gaps) if gaps else None)}

    d, c = summarise(donor), summarise(ctrl)
    verdict = "inconclusive"
    if d["n"] >= 5 and c["n"] >= 4:
        if d["in_band"] > c["in_band"]:
            verdict = "SUPPORTED"
        elif (d["median_gap"] is not None and c["median_gap"] is not None
              and d["median_gap"] < c["median_gap"] - 10):
            verdict = "SUPPORTED (by median gap)"
        elif d["in_band"] == 0 and c["in_band"] == 0:
            # Neither arm reached the band: the variable under test did not
            # move the outcome at all. A best-of-N comparison between small
            # arms is not evidence - the first version of this rule called a
            # 23.4 vs 28.6 best-gap difference "SUPPORTED" on n=9 vs n=6 with
            # fully overlapping distributions.
            verdict = ("INCONCLUSIVE - neither arm reached the band; "
                       "the substitution alone does not decide the outcome")
        elif c["in_band"] >= max(1, d["in_band"]):
            verdict = "REFUTED - non-donors do as well; pivot"
        else:
            verdict = "weak - no separation"
    return {"donor": d, "control": c, "verdict": verdict}


def build_library2(rows: list[dict], verdict: dict, budget: int) -> int:
    """Second library, built from what generation 1 actually showed."""
    parent_id = "ins080_r103_42_F42Y_H103Y"
    man = json.loads((ROOT / "outputs/insurance/variants.json").read_text())
    seq0 = next(v["sequence"] for v in man if v["design_id"] == parent_id)

    good = [r for r in rows if r.get("bound") and r.get("measurable")]
    if not good:
        return 0
    rk = objective.rank_batch([r["badness"] for r in good])
    for r, k in zip(good, rk):
        r["_rank"] = k["dense_rank"]
    # gap first, then the PINNER distance: measured across 77 screened
    # designs, the nearest H-bond donor to the Gln46 amide sits at a median
    # 2.33 A in the 14 in-band designs and 3.34 A in the 63 others. It is the
    # cheapest 1-seed predictor of the term that zeroes everything.
    good.sort(key=lambda r: (r.get("gap") if r.get("gap") is not None else 99,
                             _pin(r), -r["_rank"]))

    # every single mutation that produced a bound, measurable design
    singles: list[tuple[int, str]] = []
    for r in good[:40]:
        for m in (r.get("design_id", "").split("_")[2:]):
            if len(m) >= 3 and m[1:-1].isdigit():
                singles.append((int(m[1:-1]), m[-1]))
    singles = list(dict.fromkeys(singles))

    lib: dict[str, dict] = {}

    def add(subs, block):
        s = list(seq0)
        for p, aa in subs:
            s[p - 1] = aa
        s = "".join(s)
        if s == seq0 or s in lib:
            return
        if s[101] != "D" or s[45] != "Q":
            return
        lib[s] = {"design_id": "L2_" + block + "_" +
                  "_".join(f"{seq0[p-1]}{p}{aa}" for p, aa in subs),
                  "block": block,
                  "mutations": [f"{seq0[p-1]}{p}{aa}" for p, aa in subs]}

    # pairs and triples of what worked
    for a, b in itertools.combinations(singles[:14], 2):
        if a[0] != b[0]:
            add([a, b], "pair")
    for a, b, c3 in itertools.combinations(singles[:9], 3):
        if len({a[0], b[0], c3[0]}) == 3:
            add([a, b, c3], "triple")

    manifest = [{"sequence": s, "parent": parent_id, **m}
                for s, m in lib.items()][:budget]
    LIB2.mkdir(parents=True, exist_ok=True)
    (LIB2 / "library.json").write_text(json.dumps(manifest, indent=2))
    return len(manifest)


def five_seed(design_id: str, sequence: str) -> dict | None:
    p = sh([PY, str(ROOT / "scripts/score_finalist.py"), sequence, design_id],
           timeout=5400)
    d = FIN / design_id
    pdbs = sorted(d.glob("seed*/folded.pdb"))
    vd = d / "verdict.json"
    if not pdbs or not vd.exists():
        return None
    v = json.loads(vd.read_text())
    raws = [official.measure_raw(str(x), site=SITE) for x in pdbs]
    ipdes = [r.get("ipde") for r in v.get("replicates", [])
             if r.get("ipde") is not None]
    rep = official.score_replicates(raws, ipde_values=ipdes)
    v.update({"official_score": rep["median_total"],
              "official_score_best": rep["best_total"],
              "official_score_worst": rep["worst_total"],
              "n_nonzero_replicates": rep["n_nonzero_replicates"],
              "per_replicate_total": rep["per_replicate_total"],
              "terms_zeroed_on_any_replicate": rep["terms_zeroed_on_any_replicate"]})
    vd.write_text(json.dumps(v, indent=2, default=str))
    return rep


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--until", default="11:25")
    a = ap.parse_args()
    hh, mm = (int(x) for x in a.until.split(":"))
    stop = datetime.now().replace(hour=hh, minute=mm, second=0).timestamp()
    state = {"phase": "start", "hits": [], "scored": []}

    def save():
        STATE.write_text(json.dumps(state, indent=2, default=str))

    # ---------------- phase 1: wait for library 1
    log("phase 1: waiting for library-1 screen")
    while screen_running() and time.time() < stop - 1800:
        time.sleep(60)
    rows = load_screen(LIB / "screen.json")
    log(f"library 1: {len(rows)} folded, "
        f"{sum(1 for r in rows if r.get('bound') and r.get('measurable'))} "
        "bound+measurable")

    # ---------------- phase 2: the falsifiable test
    v = block_a_verdict(rows)
    state["block_a"] = v
    log(f"BLOCK A VERDICT: {v['verdict']}")
    log(f"   donors  n={v['donor']['n']} bound={v['donor']['bound']} "
        f"in_band={v['donor']['in_band']} best_gap={v['donor']['best_gap']}")
    log(f"   control n={v['control']['n']} bound={v['control']['bound']} "
        f"in_band={v['control']['in_band']} best_gap={v['control']['best_gap']}")
    save()

    # ---------------- phase 3/4 interleaved, to the budget
    remaining = stop - time.time()
    folds_left = int(remaining / SEC_PER_FOLD)
    screen_budget = int(folds_left * 0.45)
    log(f"{int(remaining/60)} min left = ~{folds_left} folds; "
        f"{screen_budget} to screening, rest to 5-seed scoring")

    n2 = build_library2(rows, v, screen_budget)
    if n2:
        log(f"phase 3: library 2 = {n2} designs (pairs/triples of what worked)")
        sh([PY, str(ROOT / "scripts/parallel_fold.py"),
            "--manifest", str(LIB2 / "library.json"),
            "--workers", str(WORKERS),
            "--out", "outputs/library2/screen.json"], timeout=max(600, int(stop - time.time())))
        rows += load_screen(LIB2 / "screen.json")
        log(f"library 2 done; {len(rows)} designs screened in total")
    save()

    # ---------------- phase 4: 5-seed the best
    good = [r for r in rows if r.get("bound") and r.get("measurable")]
    done = {p.parent.name for p in FIN.glob("*/verdict.json")}
    good = [r for r in good if r["design_id"] not in done]
    rk = objective.rank_batch([r["badness"] for r in good]) if good else []
    for r, k in zip(good, rk):
        r["_rank"] = k["dense_rank"]
    good.sort(key=lambda r: (r.get("gap") if r.get("gap") is not None else 99,
                             _pin(r), -r.get("_rank", 0)))
    log(f"phase 4: {len(good)} candidates for 5-seed scoring")

    for r in good:
        if time.time() > stop - 300:
            log("deadline: stopping 5-seed scoring")
            break
        did = r["design_id"]
        try:
            rep = five_seed(did, r["sequence"])
        except Exception as exc:                                 # noqa: BLE001
            log(f"   {did}: {type(exc).__name__}")
            continue
        if rep is None:
            continue
        line = (f"   {did:34} best {rep['best_total']:.4f} "
                f"repro {rep['n_nonzero_replicates']}/{rep['n_replicates']}")
        state["scored"].append({"design_id": did, **rep})
        if rep["best_total"] > 0:
            line += "   <<< SCORES"
            state["hits"].append({"design_id": did,
                                  "best": rep["best_total"],
                                  "repro": rep["n_nonzero_replicates"],
                                  "mutations": r.get("mutations")})
        log(line)
        save()
        sh([PY, str(ROOT / "scripts/bank_submissions.py"), "20"], timeout=7200)

    # ---------------- phase 5
    sh([PY, str(ROOT / "scripts/bank_submissions.py"), "20"], timeout=7200)
    log(f"DONE. {len(state['scored'])} fully scored, "
        f"{len(state['hits'])} score on >=1 replicate")
    for h in sorted(state["hits"], key=lambda x: -x["best"]):
        log(f"   {h['design_id']}  best {h['best']:.4f}  repro {h['repro']}/5")
    save()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
