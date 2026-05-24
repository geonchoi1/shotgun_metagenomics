#!/bin/bash
# === DBSCAN-SWA: prophage detection on plasmids ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FNA=$PROJECT/plasmid/02_drep/dereplicated.fna
OUT=$PROJECT/plasmid/19_dbscan_swa
mkdir -p $OUT

[ -s $OUT/bac_DBSCAN-SWA_prophage_summary.txt ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_DBSCANSWA"

# DBSCAN-SWA is provided via repo install; try both invocation styles
if command -v DBSCAN-SWA >/dev/null 2>&1; then
  DBSCAN-SWA --input $FNA --output $OUT --thread_num $THREADS --prefix plasmid
elif command -v dbscan-swa.py >/dev/null 2>&1; then
  dbscan-swa.py --input $FNA --output $OUT --thread_num $THREADS --prefix plasmid
else
  echo "ERROR: DBSCAN-SWA executable not found in env $ENV_DBSCANSWA" >&2
  exit 1
fi

echo "[$(date '+%F %T')] DBSCAN-SWA done"
