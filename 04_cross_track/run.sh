#!/bin/bash
# === 04_cross_track/run.sh ===
# Cross-track analyses requiring outputs from all 3 tracks (plasmid + MAG + UB).
set -e
BASE=$(cd "$(dirname "$0")" && pwd)

# 01 Mobile ARG (ARG×MGE co-localization + plasmid↔MAG↔UB cross-source HGT)
bash $BASE/01_mobile_arg/run.sh

echo "[$(date '+%F %T')] 04_cross_track DONE"
