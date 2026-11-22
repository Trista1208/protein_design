#!/bin/bash
# Emit one compact line per NEWLY scored design, plus terminal states.
# Used by the Monitor tool: every stdout line becomes a notification.
set -u
D=/home/mandalorian/Hackathons/Enzymes-BSSE2026
PY=/home/mandalorian/anaconda3/envs/mm/bin/python
SEEN=/tmp/kemp_seen_designs.txt
: > "$SEEN"
while true; do
  if [ -f "$D/outputs/ledger.jsonl" ]; then
    "$PY" - "$SEEN" <<'PYEOF'
import json, sys, pathlib
seen_p = pathlib.Path(sys.argv[1])
seen = set(seen_p.read_text().split()) if seen_p.exists() else set()
led = pathlib.Path("/home/mandalorian/Hackathons/Enzymes-BSSE2026/outputs/ledger.jsonl")
rows, new = {}, []
for line in led.read_text().splitlines():
    if not line.strip(): continue
    try: r = json.loads(line)
    except Exception: continue
    if r.get("design_id"): rows[r["design_id"]] = r
for did in sorted(rows):
    if did in seen: continue
    r = rows[did]
    off = r.get("official") or {}
    sub = off.get("subscores") or {}
    g = lambda k: ("-" if sub.get(k) is None else f"{sub[k]:.3f}")
    s = off.get("official_score")
    s = off.get("prefilter_score") if s is None else s
    print(f"DESIGN|{did}|{r.get('variant')}|{r.get('length')}|"
          f"{'-' if s is None else f'{s:.4f}'}|{g('base')}|{g('oxy')}|{g('ipde')}|"
          f"{g('orient_base')}|{g('orient_oxy')}|{g('bury')}|{g('cons')}|"
          f"{off.get('limiting_term') or '-'}|{r.get('eligible')}|{r.get('stage','')[:40]}",
          flush=True)
    new.append(did)
if new:
    with seen_p.open("a") as fh: fh.write("\n".join(new) + "\n")
PYEOF
  fi
  if ! pgrep -f 'kemp m[4]' >/dev/null 2>&1; then
    n=$(wc -l < "$D/outputs/ledger.jsonl" 2>/dev/null || echo 0)
    echo "SWEEP_ENDED|total_rows=$n" ; exit 0
  fi
  sleep 150
done
