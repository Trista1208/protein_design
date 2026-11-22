"""Thermal analysis: does the fold hold at elevated temperature?

    python scripts/thermal_md.py [--designs 5] [--ns 2.0] [--temp 350]

WHAT THIS IS
------------
Implicit-solvent (GBn2) molecular dynamics at 350 K, protein only. We measure
how far the backbone drifts from its starting fold and whether the chain stays
compact. A design that unfolds at 350 K shows RMSD climbing without plateau
and radius of gyration rising; one that holds shows both flat after
equilibration.

CONTROL ARM: the HG3.17 (4BS0) crystal structure is run under identical
settings. Its melting temperature is published (~70 C / 343 K), so it tells us
what "holds at 350 K" looks like for a protein that demonstrably folds and
functions. Without it these RMSD numbers mean nothing (good_practices 5b).

WHAT IT IS NOT
--------------
This is NOT a melting temperature and NOT a free energy. Two nanoseconds of
implicit-solvent MD cannot unfold a stable protein, so a flat trace is weak
evidence of stability and a rising one is strong evidence of a problem. Read
it as a screen for gross instability, not as a Tm.

Implicit solvent also has no viscosity from explicit water, so absolute
timescales are not physical; only the comparison between runs is meaningful,
which is why every run uses identical settings.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/stability"
REF = ROOT / "outputs/runs/controls_full/4BS0_A.pdb"


def protein_only(src: Path, dst: Path) -> None:
    with dst.open("w") as f:
        for ln in src.read_text().splitlines():
            if ln.startswith(("ATOM", "TER")):
                f.write(ln + "\n")
        f.write("END\n")


def simulate(pdb_in: Path, tag: str, ns: float, temp: float) -> dict | None:
    from openmm import app, unit, LangevinMiddleIntegrator, Platform
    from pdbfixer import PDBFixer

    work = OUT / "md" / tag
    work.mkdir(parents=True, exist_ok=True)
    stripped = work / "prot.pdb"
    protein_only(pdb_in, stripped)

    fixer = PDBFixer(filename=str(stripped))
    fixer.findMissingResidues()
    fixer.missingResidues = {}           # do not build loops that were not solved
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    fixed = work / "fixed.pdb"
    with fixed.open("w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)

    pdb = app.PDBFile(str(fixed))
    ff = app.ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")
    system = ff.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds,
                             hydrogenMass=1.5 * unit.amu)
    integ = LangevinMiddleIntegrator(temp * unit.kelvin,
                                     1.0 / unit.picosecond,
                                     0.004 * unit.picoseconds)
    try:
        plat = Platform.getPlatformByName("CUDA")
        sim = app.Simulation(pdb.topology, system, integ, plat)
    except Exception:                                            # noqa: BLE001
        sim = app.Simulation(pdb.topology, system, integ)
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(temp * unit.kelvin)
    sim.step(25000)                      # 100 ps equilibration, discarded

    # CA indices for RMSD / Rg
    ca = [a.index for a in pdb.topology.atoms() if a.name == "CA"]
    st = sim.context.getState(getPositions=True)
    ref = np.array(st.getPositions().value_in_unit(unit.nanometer))[ca]
    ref -= ref.mean(axis=0)

    n_steps = int(ns * 1000 / 0.004)
    chunk = 2500                          # 10 ps
    rmsd, rg, t = [], [], []
    for i in range(n_steps // chunk):
        sim.step(chunk)
        s = sim.context.getState(getPositions=True)
        x = np.array(s.getPositions().value_in_unit(unit.nanometer))[ca]
        x = x - x.mean(axis=0)
        # Kabsch superposition onto the post-equilibration reference
        u, _, vt = np.linalg.svd(x.T @ ref)
        d = np.sign(np.linalg.det(u @ vt))
        r = u @ np.diag([1, 1, d]) @ vt
        xa = x @ r
        rmsd.append(float(np.sqrt(((xa - ref) ** 2).sum(axis=1).mean())) * 10)
        rg.append(float(np.sqrt((x ** 2).sum(axis=1).mean())) * 10)
        t.append((i + 1) * chunk * 0.004)

    half = len(rmsd) // 2
    return {"tag": tag, "n_ca": len(ca), "temp_K": temp, "ns": ns,
            "rmsd_mean_last_half": float(np.mean(rmsd[half:])),
            "rmsd_max": float(np.max(rmsd)),
            "rmsd_final": float(rmsd[-1]),
            "rmsd_drift_per_ns": float(
                np.polyfit(t[half:], rmsd[half:], 1)[0] * 1000),
            "rg_mean_last_half": float(np.mean(rg[half:])),
            "rg_start": float(rg[0]),
            "rmsd_series": rmsd, "rg_series": rg, "time_ps": t}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--designs", type=int, default=5)
    ap.add_argument("--ns", type=float, default=2.0)
    ap.add_argument("--temp", type=float, default=350.0)
    a = ap.parse_args()

    top = json.loads((ROOT / "outputs/top20.json").read_text())[:a.designs]
    targets = []
    for r in top:
        d = ROOT / "outputs/finalists" / r["id"]
        per = json.loads((d / "verdict.json").read_text()).get("per_replicate_total") or []
        pdbs = sorted(d.glob("seed*/folded.pdb"))
        if pdbs:
            i = max(range(len(per)), key=lambda k: per[k]) if per else 0
            targets.append((pdbs[min(i, len(pdbs) - 1)], r["id"]))
    if REF.exists():
        targets.append((REF, "REF_HG3.17_4BS0"))

    print(f"{len(targets)} runs, {a.ns} ns each at {a.temp} K, implicit solvent",
          flush=True)
    rows = []
    for p, tag in targets:
        t0 = time.time()
        try:
            r = simulate(p, tag, a.ns, a.temp)
        except Exception as exc:                                 # noqa: BLE001
            print(f"  {tag}: {type(exc).__name__}: {exc}", flush=True)
            continue
        if r:
            rows.append(r)
            print(f"  {tag:30} RMSD {r['rmsd_mean_last_half']:5.2f} A  "
                  f"max {r['rmsd_max']:5.2f}  drift {r['rmsd_drift_per_ns']:+5.2f} A/ns  "
                  f"Rg {r['rg_start']:5.2f}->{r['rg_mean_last_half']:5.2f} A  "
                  f"({time.time()-t0:.0f}s)", flush=True)
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / "thermal_md.json").write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {OUT/'thermal_md.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
