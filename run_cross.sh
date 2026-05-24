#!/bin/bash
# === run_cross.sh — 04_cross_track ===
set -euo pipefail
REPO="$(cd "$(dirname "$0")" && pwd)"
source "$REPO/config.sh"
PROJECT=${1:?ERROR: usage: bash run_cross.sh /path/to/project}; export PROJECT
T="$REPO/04_cross_track"

echo "[$(date '+%F %T')] === 04_cross_track START ==="
bash "$T/01_mobile_arg/run.sh"
echo "[$(date '+%F %T')] === 04_cross_track DONE ==="
