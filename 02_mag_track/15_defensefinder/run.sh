#!/bin/bash
# === 15 DefenseFinder (anti-phage systems) ===
# Input:  $PROJECT/mag/03_master_orf/all/master.faa
# Output: $PROJECT/mag/15_defensefinder/defense_finder_systems.tsv (+ genes.tsv, hmmer.tsv)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/mag/03_master_orf/all/master.faa
OUT=$PROJECT/mag/15_defensefinder
mkdir -p "$OUT"

[ -s "$FAA" ] || { echo "ERROR: $FAA empty" >&2; exit 1; }

activate_env "$ENV_DEFENSEFINDER"

echo "[$(date '+%F %T')] defense-finder run"
defense-finder run \
    --models-dir "$DEFENSEFINDER_MODELS" \
    --out-dir "$OUT" \
    --workers "$THREADS" \
    "$FAA"

n=$(tail -n +2 "$OUT"/*defense_finder_systems.tsv 2>/dev/null | wc -l || echo 0)
echo "[$(date '+%F %T')] DONE — $n defense systems"
