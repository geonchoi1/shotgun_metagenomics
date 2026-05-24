#!/bin/bash
# === 08 AMRFinderPlus --plus ===
# Input:  $PROJECT/mag/03_master_orf/all/master.faa
# Output: $PROJECT/mag/08_amrfinder/amrfinder.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/mag/03_master_orf/all/master.faa
OUT=$PROJECT/mag/08_amrfinder
mkdir -p "$OUT"

[ -s "$FAA" ] || { echo "ERROR: $FAA empty" >&2; exit 1; }

activate_env "$ENV_AMRFINDER"

echo "[$(date '+%F %T')] amrfinder --plus -p (protein)"
amrfinder \
    -p "$FAA" \
    --plus \
    --threads "$THREADS" \
    -d "$AMRFINDER_DB" \
    -o "$OUT/amrfinder.tsv"

n=$(tail -n +2 "$OUT/amrfinder.tsv" 2>/dev/null | wc -l)
echo "[$(date '+%F %T')] DONE — $n AMR/virulence/stress hits"
