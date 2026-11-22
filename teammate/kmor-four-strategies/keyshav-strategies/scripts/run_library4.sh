#!/bin/bash
# Screen library 4 at 1 seed, annotate, then 5-seed the in-band survivors.
set -u
PY=/home/mandalorian/anaconda3/envs/mm/bin/python
cd /home/mandalorian/Hackathons/Enzymes-BSSE2026
echo "[$(date +%H:%M)] 1-seed screen of library4"
$PY scripts/parallel_fold.py --manifest outputs/library4/library.json \
    --workers 8 --out outputs/library4/screen.json
echo "[$(date +%H:%M)] annotate"
$PY scripts/annotate_screens.py
echo "[$(date +%H:%M)] 5-seed the in-band set"
$PY scripts/score_inband.py --until 23:30 --max 40
echo "[$(date +%H:%M)] done"
