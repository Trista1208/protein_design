#!/bin/bash
# Wait for the calibration gate, then run the sweep. Fully detached (setsid), so
# it survives the session ending - plain nohup did not.
#
#   bash scripts/launch_campaign.sh 50          # n designs
#   bash scripts/launch_campaign.sh 500 --force # run uncalibrated (not advised)
set -u
D=/home/mandalorian/Hackathons/Enzymes-BSSE2026
PY=/home/mandalorian/anaconda3/envs/foundry/bin/python
N="${1:-50}"
FORCE="${2:-}"
cd "$D" || exit 1
export KEMP_ROOT="$D"

say(){ echo "[$(date '+%F %T')] $*"; }

say "waiting for the controls job to finish"
w=0
while pgrep -f 'kemp control[s]' >/dev/null 2>&1; do
  sleep 20; w=$((w+20))
  [ $((w % 300)) -eq 0 ] && say "  still folding controls (${w}s)"
  [ $w -ge 5400 ] && { say "  TIMEOUT waiting on controls"; break; }
done

say "controls verdict:"
$PY - <<'PYEOF'
import json, pathlib
p = pathlib.Path("outputs/state.json")
st = json.loads(p.read_text()) if p.exists() else {}
c = st.get("controls") or {}
for k in ("positive_median", "negative_median", "separation", "noise_floor",
          "ranks_correctly", "separation_exceeds_noise", "lower_is_better"):
    if k in c:
        print(f"    {k:26} {c[k]}")
print(f"    verdict: {str(c.get('verdict'))[:150]}")
PYEOF

say "starting the sweep: n=$N ${FORCE}"
$PY -m kemp m4 --n "$N" --timesteps 200 --seqs 8 --finalise --top 30 $FORCE \
   2>&1 | tee -a "$D/logs/campaign_n${N}.log"
RC=${PIPESTATUS[0]}
say "sweep finished rc=$RC"

if [ "$RC" -eq 0 ]; then
  say "selecting final designs"
  $PY -m kemp m5 2>&1 | tee -a "$D/logs/campaign_n${N}.log"
fi
say "DONE"
