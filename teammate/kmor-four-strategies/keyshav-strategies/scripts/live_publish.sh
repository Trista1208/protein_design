#!/bin/bash
# Live publisher: regenerate table/plots/structures and push to the branch.
#
# Runs until the sweep finishes AND one final publish has been done. Designed to
# survive logout: started with setsid+nohup, and `loginctl enable-linger` keeps
# user processes alive after the session ends.
#
#   bash scripts/live_publish.sh [interval_seconds]
set -u
D=/home/mandalorian/Hackathons/Enzymes-BSSE2026
T=/home/mandalorian/Hackathons/team-sharing
DST=$T/keyshav-strategies
BRANCH=kmor/four-strategies
REMOTE=git@github.com:nikafiltse/enzyme-design-hackathon-2026.git
PY=/home/mandalorian/anaconda3/envs/mm/bin/python
INTERVAL="${1:-120}"
cd "$D" || exit 1
export KEMP_ROOT="$D"

say(){ echo "[$(date '+%F %T')] $*"; }

publish_once() {
  "$PY" - <<'PYEOF' 2>&1 | tail -3
import sys; sys.path.insert(0, "/home/mandalorian/Hackathons/Enzymes-BSSE2026")
from kemp import report
import json
print(json.dumps(report.publish()))
PYEOF
}

sync_and_push() {
  mkdir -p "$DST"
  # code + docs (cheap, keeps the branch coherent with the results)
  for item in kemp docs tests scripts envs Makefile README.md .gitignore; do
    cp -r "$D/$item" "$DST/" 2>/dev/null
  done
  rm -rf "$DST/results"
  cp -r "$D/published" "$DST/results" 2>/dev/null
  find "$DST" -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null
  find "$DST" -name '.pytest_cache' -type d -exec rm -rf {} + 2>/dev/null

  cd "$T" || return 1
  git add -A keyshav-strategies >/dev/null 2>&1
  if git diff --cached --quiet 2>/dev/null; then
    cd "$D"; return 0        # nothing changed
  fi
  local n best
  n=$("$PY" -c "
import json,pathlib
p=pathlib.Path('$D/published/status.json')
d=json.loads(p.read_text()) if p.exists() else {}
print(d.get('n_total',0))" 2>/dev/null || echo 0)
  best=$("$PY" -c "
import json,pathlib
p=pathlib.Path('$D/published/status.json')
d=json.loads(p.read_text()) if p.exists() else {}
b=d.get('best')
print('n/a' if b is None else round(float(b),4))" 2>/dev/null || echo n/a)
  git -c user.email=keyshav.mor@lemna.bio commit -q -m \
"live: ${n} designs, best S=${best} (auto-published $(date '+%F %T'))

Table: keyshav-strategies/results/README.md
Plots: keyshav-strategies/results/plots/  (one panel per official metric)
PDBs:  keyshav-strategies/results/structures/
Yardstick: HG3.17 crystal 0.9591, Q50M knockout 0.0000

Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>" >/dev/null 2>&1
  if git push "$REMOTE" "$BRANCH" >/dev/null 2>&1; then
    say "pushed: ${n} designs, best S=${best}"
  else
    say "PUSH FAILED (will retry next cycle)"
  fi
  cd "$D"
}

say "live publisher started (interval ${INTERVAL}s, branch ${BRANCH})"
idle=0
while true; do
  publish_once >/dev/null 2>&1
  sync_and_push
  # stop once the sweep is gone and we have published twice more
  if ! pgrep -f 'kemp m[4]' >/dev/null 2>&1; then
    idle=$((idle+1))
    if [ "$idle" -ge 2 ]; then
      say "sweep finished and final publish done - exiting"
      break
    fi
  else
    idle=0
  fi
  sleep "$INTERVAL"
done
say "DONE"
