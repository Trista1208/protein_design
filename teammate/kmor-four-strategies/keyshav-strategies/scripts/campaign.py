"""Generational design campaign, time-budgeted, runs unattended.

    python scripts/campaign.py --until 11:00 [--parent <design_id>]

WHAT IT DOES EACH GENERATION
----------------------------
  propose -> fold (1 seed) -> filter -> rank -> take the winners as parents

PROPOSAL SOURCES, in priority order, each with its reason:

  1. R127 SCAN. Arg127 sits 3.3-4.0 A from the catalytic Asp102 in every
     design we have. Both crystal controls (4BS0, 7K4Z) have ZERO charged
     residues within 4 A of their base. This is Privett/Alexandrova hard
     filter H3, and it has a worked precedent: in KE07 a Lys222-base salt
     bridge was removed and kcat/KM improved 2.5-fold; in KE16 the
     Lys201-Asp48 salt bridge competed with the substrate in the apo form.
     Rosetta's untargeted enzdes run independently chose R127I/L/M.
     Two unrelated methods pointing at one residue is the strongest signal
     we have.

  2. HYDROPHOBIC SHELL. Privett soft filter S2 wants >= 0.60 apolar residues
     in the 5-6 A shell around the base; a nonpolar pocket raises the base
     pKa (~6 vs ~4 when exposed). The controls sit at 0.727. ins080 is at
     0.500 and fails.

  3. ROTAMER SHELL. The positions that move the oxyanion dihedral, measured
     on an 84-variant screen: 23, 26, 27, 42, 45, 49, 84, 99, 103.

  4. COMBINATIONS of whatever won the previous generation. This is the one
     that has actually worked: F42Y and H103Y each failed alone and together
     produced the first design with six of seven terms non-zero. Bhattacharya
     2022 reports the same - a triple mutant 3x more active than additivity
     predicted.

FILTERS, in this order and no other (docs/good_practices.md §2):
  bound ligand -> pLDDT >= C.PLDDT_MIN -> literature filters -> dense rank.
A dihedral measured on a drifted ligand is noise, and that filter has already
overturned three apparent breakthroughs.

NOT USED: Rosetta score as an activity ranker. Risso et al. 2020 tested it on
this exact reaction and report "No correlation with activity. Good for
stability ranking, not for catalysis."
"""
from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kemp import config as C, litfilters, objective, official, pipeline  # noqa: E402
from kemp.runners import RunnerError, boltz_predict, ligandmpnn          # noqa: E402
from kemp.theozyme import load                                           # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/campaign"
BASE_IDX, OXY_IDX, CHAIN = 102, 46, "A"
IN_WINDOW = (2.0, 3.6)

SEQUESTER_POS = 127                     # the Arg that salt-bridges the base
SEQUESTER_SUBS = list("AVLIMFQNSTGY")   # remove the charge, vary the bulk
ROTAMER_SHELL = [26, 103, 42, 99, 49, 45, 23, 27, 84]
LADDER = ["A", "V", "L", "I", "F", "W", "M", "T", "Y"]
OMIT = "C"


def bound(d_base, d_oxy) -> bool:
    return (d_base is not None and d_oxy is not None
            and IN_WINDOW[0] <= d_base <= IN_WINDOW[1]
            and IN_WINDOW[0] <= d_oxy <= IN_WINDOW[1])


def mutate(seq: str, subs: list[tuple[int, str]]) -> str:
    s = list(seq)
    for pos, aa in subs:
        s[pos - 1] = aa
    return "".join(s)


def diff(a: str, b: str) -> list[str]:
    return [f"{a[i]}{i+1}{b[i]}" for i in range(min(len(a), len(b)))
            if a[i] != b[i]]


def propose(parent_seq: str, parent_pdb: Path | None, gen: int,
            winners: list[dict], budget: int) -> dict[str, dict]:
    """Build this generation's variant set."""
    out: dict[str, dict] = {}

    def add(seq, tag, source):
        if seq == parent_seq or seq in out:
            return
        if seq[BASE_IDX - 1] != "D" or seq[OXY_IDX - 1] != "Q":
            return
        out[seq] = {"tag": tag, "source": source}

    # 1. the sequestering Arg - highest-priority, literature-backed
    if parent_seq[SEQUESTER_POS - 1] in "RKH":
        for aa in SEQUESTER_SUBS:
            add(mutate(parent_seq, [(SEQUESTER_POS, aa)]),
                f"seq{SEQUESTER_POS}{aa}", "sequester_fix")

    # 2/3. rotamer + hydrophobic shell ladder
    for pos in ROTAMER_SHELL:
        for aa in LADDER:
            if parent_seq[pos - 1] != aa:
                add(mutate(parent_seq, [(pos, aa)]),
                    f"l{pos}{aa}", "ladder")

    # 4. combinations of the previous generation's winning single edits
    singles = []
    for w in winners:
        for m in w.get("mutations", []):
            if len(m) >= 3:
                singles.append((int(m[1:-1]), m[-1]))
    singles = list(dict.fromkeys(singles))
    for a, b in itertools.combinations(singles[:10], 2):
        if a[0] != b[0]:
            add(mutate(parent_seq, [a, b]),
                f"c{a[0]}{a[1]}_{b[0]}{b[1]}", "combination")

    # 5. model proposals, native forbidden (it returns the parent otherwise)
    if parent_pdb and parent_pdb.exists() and len(out) < budget:
        for pos in ROTAMER_SHELL[:5] + [SEQUESTER_POS]:
            if len(out) >= budget:
                break
            w = OUT / f"gen{gen}" / "mpnn" / f"p{pos}"
            w.mkdir(parents=True, exist_ok=True)
            try:
                res = ligandmpnn(parent_pdb, w,
                                 redesigned_residues=f"{CHAIN}{pos}",
                                 omit_per_residue={
                                     f"{CHAIN}{pos}": parent_seq[pos - 1] + OMIT},
                                 n_seqs=8, n_batches=2, temperature=0.3, seed=37)
                if res.returncode == 0:
                    for _h, s in res.sequences:
                        if len(s) == len(parent_seq) and len(diff(parent_seq, s)) == 1:
                            add(s, f"m{pos}", "ligandmpnn")
            except RunnerError:
                pass

    # Priority, not a lottery. A plain shuffle dropped 4 of the 12 R127
    # variants from the first generation - the one proposal class with both a
    # literature hard filter and an independent Rosetta vote behind it.
    # Combinations rank next because they are the only thing that has produced
    # a real gain so far (F42Y+H103Y).
    rank = {"sequester_fix": 0, "combination": 1, "ligandmpnn": 2, "ladder": 3}
    items = sorted(out.items(), key=lambda kv: rank.get(kv[1]["source"], 9))
    head = [kv for kv in items if rank.get(kv[1]["source"], 9) <= 1]
    tail = [kv for kv in items if rank.get(kv[1]["source"], 9) > 1]
    random.Random(gen).shuffle(tail)        # diversify only within the filler
    return dict((head + tail)[:budget])


def fold_and_score(seq: str, workdir: Path, site) -> dict | None:
    r = boltz_predict(seq, workdir, seed=C.BOLTZ_SEEDS[0],
                      ligand_ccd=C.LIGAND_CCD)
    f = next((x for x in r.structures), None)
    if f is None:
        return None
    p = f if f.suffix == ".pdb" else pipeline.cif_to_pdb(f, workdir / "folded.pdb")
    raw = official.measure_raw(str(p), site=site)
    conf = r.confidence or {}
    raw.ipde = conf.get("complex_ipde")
    plddt = conf.get("complex_plddt") or 0.0
    rec = {"pdb": str(p), "d_base": raw.d_base, "d_oxy": raw.d_oxy,
           "dih_oxy": raw.dihedral12_oxy,
           "gap": objective.dihedral_band_gap(raw.dihedral12_oxy),
           "plddt": plddt, "ipde": raw.ipde,
           "bound": bound(raw.d_base, raw.d_oxy),
           "measurable": plddt >= C.PLDDT_MIN,
           "subscores": official.score(raw, strict=False)["subscores"],
           "badness": objective.badness(raw)}
    if rec["bound"] and rec["measurable"]:
        try:
            lf = litfilters.evaluate(load(str(p), require_ligand=False),
                                     base_resseq=BASE_IDX, oxy_resseq=OXY_IDX,
                                     chain=CHAIN)
            rec["lit"] = {"values": lf.values, "hard_fail": lf.hard_fail,
                          "soft_pass": lf.soft_pass}
        except Exception as exc:                                # noqa: BLE001
            rec["lit"] = {"error": str(exc)[:120]}
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--until", default="11:00", help="stop time, HH:MM local")
    ap.add_argument("--parent-seq-id", default="ins080_r103_42_F42Y_H103Y")
    ap.add_argument("--parent-pdb",
                    default="outputs/finalists/ins080_F42Y_H103Y/seed1/folded.pdb")
    ap.add_argument("--per-gen", type=int, default=60)
    a = ap.parse_args()

    hh, mm = (int(x) for x in a.until.split(":"))
    now = datetime.now()
    deadline = now.replace(hour=hh, minute=mm, second=0, microsecond=0).timestamp()
    if deadline <= time.time():
        print("deadline already passed")
        return 2

    man = json.loads((ROOT / "outputs/insurance/variants.json").read_text())
    parent = next((v for v in man if v["design_id"] == a.parent_seq_id), None)
    if parent is None:
        print(f"parent {a.parent_seq_id} not in variants.json")
        return 2
    parent_seq = parent["sequence"]
    parent_pdb = ROOT / a.parent_pdb
    site = official.CatalyticSite(base_resseq=BASE_IDX, oxy_resseq=OXY_IDX,
                                  chain=CHAIN)
    OUT.mkdir(parents=True, exist_ok=True)

    print(f"campaign start {now:%H:%M}, stop by {a.until}")
    print("GOAL: ins080 reaches official 0.743 on 1 seed of 5. Its catalytic")
    print("base contact is stable (d_base sd 0.073, HG3.17 ref 0.068); what")
    print("varies is the Gln46 rotamer, which samples ~48 / ~95 / ~247 deg and")
    print("only scores at ~48. Make that state DOMINANT, not occasional.")
    print(f"parent {a.parent_seq_id}  R{SEQUESTER_POS}="
          f"{parent_seq[SEQUESTER_POS-1]} (salt-bridges the base: fix first)\n")

    all_rows: list[dict] = []
    best_overall: list[dict] = []
    winners: list[dict] = []
    gen = 0
    while time.time() < deadline - 300:
        gen += 1
        remaining = int((deadline - time.time()) / 45)
        budget = max(10, min(a.per_gen, remaining - 20))
        variants = propose(parent_seq, parent_pdb, gen, winners, budget)
        if not variants:
            print("no further proposals")
            break
        print(f"--- generation {gen}: {len(variants)} variants "
              f"({int((deadline-time.time())/60)} min left) ---", flush=True)

        rows = []
        for i, (seq, meta) in enumerate(variants.items()):
            if time.time() > deadline - 120:
                print("  deadline reached mid-generation")
                break
            muts = diff(parent_seq, seq)
            did = f"g{gen}_{meta['tag']}"
            w = OUT / f"gen{gen}" / did
            w.mkdir(parents=True, exist_ok=True)
            try:
                rec = fold_and_score(seq, w, site)
            except Exception as exc:                            # noqa: BLE001
                print(f"  [{i+1}] {did}: {type(exc).__name__}")
                continue
            if rec is None:
                continue
            rec.update({"design_id": did, "gen": gen, "sequence": seq,
                        "mutations": muts, "parent": a.parent_seq_id, **meta})
            rows.append(rec)
            all_rows.append(rec)
            flag = ""
            if rec["bound"] and rec["measurable"]:
                lit = rec.get("lit", {})
                if not lit.get("hard_fail"):
                    flag += "  LIT-PASS"
                if rec["gap"] == 0.0:
                    flag += "  <<< IN BAND"
            print(f"  [{i+1}/{len(variants)}] {did:22} {'+'.join(muts)[:26]:26} "
                  f"db {rec['d_base'] or 0:5.2f} do {rec['d_oxy'] or 0:5.2f} "
                  f"dih {rec['dih_oxy'] or 0:6.1f} gap {rec['gap'] or 0:5.1f} "
                  f"pL {rec['plddt']:.3f}{flag}", flush=True)
            (OUT / "campaign.json").write_text(
                json.dumps(all_rows, indent=2, default=str))

        good = [r for r in rows if r["bound"] and r["measurable"]]
        if not good:
            print(f"  generation {gen}: nothing bound+measurable; "
                  "keeping the same parent")
            winners = []
            continue
        rk = objective.rank_batch([r["badness"] for r in good])
        for r, k in zip(good, rk):
            r["dense_rank"] = k["dense_rank"]
        good.sort(key=lambda r: -r["dense_rank"])
        winners = good[:8]
        best_overall.extend(good)
        print(f"  generation {gen}: {len(good)} bound+measurable; "
              f"best rank {good[0]['dense_rank']:.3f} "
              f"({'+'.join(good[0]['mutations'])})")

        # promote: the new parent is the best design that also clears the
        # literature hard filters, else simply the best ranked
        lit_ok = [r for r in good if not r.get("lit", {}).get("hard_fail")]
        promote = (lit_ok or good)[0]
        if promote["dense_rank"] > 0.5:
            parent_seq = promote["sequence"]
            parent_pdb = Path(promote["pdb"])
            print(f"  -> new parent: {promote['design_id']} "
                  f"({'+'.join(promote['mutations'])})")

    best_overall.sort(key=lambda r: -r.get("dense_rank", 0))
    (OUT / "best.json").write_text(json.dumps(best_overall[:40], indent=2,
                                              default=str))
    print(f"\n{len(all_rows)} folded over {gen} generation(s); "
          f"{len(best_overall)} bound+measurable")
    print(f"top candidates -> {OUT/'best.json'}")
    for r in best_overall[:12]:
        lit = r.get("lit", {})
        print(f"  {r['design_id']:22} {'+'.join(r['mutations'])[:30]:30} "
              f"gap {r['gap']:5.1f}  rank {r.get('dense_rank', 0):.3f}  "
              f"lit={'PASS' if not lit.get('hard_fail') else 'fail'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
