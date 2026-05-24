#!/bin/bash
# === AMRFinderPlus on master.faa (protein-only) ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/plasmid/04_master_orf/all/master.faa
OUT=$PROJECT/plasmid/09_amrfinder
mkdir -p $OUT

[ -s $OUT/amrfinder.tsv ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_AMRFINDER"

amrfinder --plus \
  -p $FAA \
  -d $AMRFINDER_DB \
  --threads $THREADS \
  -o $OUT/amrfinder.tsv \
  > $OUT/amrfinder.log 2>&1

echo "[$(date '+%F %T')] AMR/SCC hits: $(($(wc -l < $OUT/amrfinder.tsv) - 1))"
