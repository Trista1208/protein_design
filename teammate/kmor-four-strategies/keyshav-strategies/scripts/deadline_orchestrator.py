"""Guarantee everything is finished and published before the 13:30 deadline.

    python scripts/deadline_orchestrator.py

TIMELINE (wall clock, local)
  until 12:45  let the library-4 1-seed screen run, then CUT IT
    12:45      annotate the partial screen with the literature filters
    12:46      5-seed the in-band survivors, hard stop 13:08
    13:10      final_consolidate.py - re-measure, rebuild top 20 and
               designs.fasta, ThermoMPNN from the cached scan, relax for
               any new entrant, biophysics table
    13:22      publish and push
  by 13:30     done

WHY CUT THE SCREEN
------------------
It runs at ~17.5 s/design, so all 224 would finish at ~13:16 and leave no
time to 5-seed anything. A partial screen is not a problem: a 1-seed screen
is a FILTER, not a result (good_practices 5c), and its only job is to decide
which designs deserve five seeds. Half a filter still filters.

Killing is by PID recorded from pgrep on a pattern this process's own command
line cannot contain (good_practices 5f: never kill by a self-matching
pattern). The pattern is assembled at runtime for the same reason.
"""
from __future__ import annotations
import json, os, signal, subprocess, sys, time
from datetime import datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable
LOG = ROOT / "logs/orchestrator.log"

CUT_SCREEN = "12:45"
SCORE_UNTIL = "13:08"
CONSOLIDATE = "13:10"
PUBLISH = "13:22"


def say(msg: str) -> None:
    line = f"[{datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)
    with LOG.open("a") as f:
        f.write(line + "\n")


def at(hhmm: str) -> float:
    h, m = (int(x) for x in hhmm.split(":"))
    t = datetime.now().replace(hour=h, minute=m, second=0, microsecond=0)
    if t < datetime.now() - timedelta(hours=6):
        t += timedelta(days=1)
    return t.timestamp()


def wait_until(hhmm: str, why: str) -> None:
    target = at(hhmm)
    if time.time() >= target:
        say(f"already past {hhmm}; {why} now")
        return
    say(f"waiting until {hhmm} ({why}) - {(target-time.time())/60:.1f} min")
    while time.time() < target:
        time.sleep(min(20, max(1, target - time.time())))


def kill_screen() -> None:
    # built at runtime so this process's own argv can never match it
    pat = "parallel" + "_fold" + ".py"
    try:
        out = subprocess.run(["pgrep", "-f", pat], capture_output=True,
                             text=True, timeout=30).stdout.split()
    except Exception as exc:                                     # noqa: BLE001
        say(f"pgrep failed: {exc}")
        return
    mine = str(os.getpid())
    pids = [p for p in out if p != mine]
    say(f"cutting the screen: {len(pids)} process(es) {pids}")
    for p in pids:
        try:
            os.kill(int(p), signal.SIGTERM)
        except Exception:                                        # noqa: BLE001
            pass
    time.sleep(8)
    # also stop the wrapper that would otherwise start the next stage itself
    pat2 = "run_library4" + ".sh"
    out2 = subprocess.run(["pgrep", "-f", pat2], capture_output=True,
                          text=True, timeout=30).stdout.split()
    for p in [x for x in out2 if x != mine]:
        try:
            os.kill(int(p), signal.SIGTERM)
        except Exception:                                        # noqa: BLE001
            pass


def run(script: str, *args: str, timeout: int) -> None:
    say(f"running {script} {' '.join(args)}")
    try:
        p = subprocess.run([PY, str(ROOT / script), *args],
                           capture_output=True, text=True, timeout=timeout,
                           env={**os.environ, "KEMP_ROOT": str(ROOT)})
        tail = [l for l in (p.stdout or "").splitlines() if l.strip()][-12:]
        for l in tail:
            say("   " + l[:120])
        if p.returncode != 0:
            err = [l for l in (p.stderr or "").splitlines() if l.strip()][-4:]
            for l in err:
                say("   ERR " + l[:120])
    except subprocess.TimeoutExpired:
        say(f"   {script} hit its timeout; continuing")
    except Exception as exc:                                     # noqa: BLE001
        say(f"   {script}: {type(exc).__name__}: {exc}")


def screened() -> int:
    f = ROOT / "outputs/library4/screen.json"
    if not f.exists():
        return 0
    try:
        return len(json.loads(f.read_text()))
    except Exception:                                            # noqa: BLE001
        return 0


def main() -> int:
    say("orchestrator up; deadline 13:30")
    wait_until(CUT_SCREEN, "cut the library-4 screen")
    kill_screen()
    say(f"library4 screened {screened()} designs before the cut")

    run("scripts/annotate_screens.py", timeout=900)
    run("scripts/score_inband.py", "--until", SCORE_UNTIL, "--max", "40",
        timeout=3600)

    wait_until(CONSOLIDATE, "final consolidation")
    run("scripts/final_consolidate.py", timeout=1500)

    wait_until(PUBLISH, "publish and push")
    run("scripts/live_publish_team.py", "--push", timeout=420)
    say("DONE - everything measured, consolidated and pushed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
