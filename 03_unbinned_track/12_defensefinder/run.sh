#!/bin/bash
# === DefenseFinder: bacterial defense systems on master.faa ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/unbinned/03_master_orf/master.faa
OUT=$PROJECT/unbinned/12_defensefinder
mkdir -p $OUT

[ -s $OUT/defense_finder_systems.tsv ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_DEFENSEFINDER"

defense-finder run \
  --models-dir $DEFENSEFINDER_MODELS \
  --workers $THREADS \
  --out-dir $OUT \
  $FAA

echo "[$(date '+%F %T')] Defense systems: $(($(wc -l < $OUT/defense_finder_systems.tsv 2>/dev/null || echo 1) - 1))"
