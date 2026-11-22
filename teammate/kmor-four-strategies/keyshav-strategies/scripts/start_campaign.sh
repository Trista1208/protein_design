#!/bin/bash
# Start the design sweep AND the live publisher, both fully detached.
#
#   bash scripts/start_campaign.sh 100          # n designs
#   bash scripts/start_campaign.sh 500 60       # n designs, publish every 60 s
#
# Survives logout three ways:
#   * setsid        - new session, no controlling terminal
#   * nohup         - ignores SIGHUP
#   * enable-linger - systemd keeps user processes after the session ends
# (plain `nohup ... &` was NOT enough earlier: the job died at session teardown.)
set -u
D=/home/mandalorian/Hackathons/Enzymes-BSSE2026
PY=/home/mandalorian/anaconda3/envs/foundry/bin/python
N="${1:-100}"
INTERVAL="${2:-120}"
cd "$D" || exit 1
export KEMP_ROOT="$D"
mkdir -p logs

say(){ echo "[$(date '+%F %T')] $*"; }

if pgrep -f 'kemp m[4]' >/dev/null 2>&1; then
  say "a sweep is already running - refusing to start a second"
  exit 1
fi

loginctl enable-linger "$USER" 2>/dev/null
say "linger: $(loginctl show-user "$USER" --property=Linger 2>/dev/null)"

say "starting sweep: n=$N"
setsid nohup "$PY" -m kemp m4 --n "$N" --timesteps 200 --seqs 8 \
    --finalise --top 30 \
    > "$D/logs/campaign_n${N}.log" 2>&1 < /dev/null &
disown 2>/dev/null || true
sleep 3
say "sweep pid group: $(pgrep -f 'kemp m[4]' | tr '\n' ' ')"

say "starting live publisher (every ${INTERVAL}s)"
setsid nohup bash "$D/scripts/live_publish.sh" "$INTERVAL" \
    > "$D/logs/live_publish.log" 2>&1 < /dev/null &
disown 2>/dev/null || true
sleep 2
say "publisher pid: $(pgrep -f 'live_publis[h]' | tr '\n' ' ')"

cat <<EOF

  Both detached. Safe to log out.

  sweep log      tail -f $D/logs/campaign_n${N}.log
  publisher log  tail -f $D/logs/live_publish.log
  local results  $D/published/
  branch         https://github.com/nikafiltse/enzyme-design-hackathon-2026/tree/kmor/four-strategies
                 -> keyshav-strategies/results/

  stop:  pkill -f 'kemp m[4]' ; pkill -f 'live_publis[h]'
EOF
