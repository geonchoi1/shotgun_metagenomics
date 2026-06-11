#!/bin/bash
# === 18 DBSCAN-SWA (prophage detection) on UB contigs ===
# Input:  $PROJECT/03_unbinned_track/01_raw_fasta/all/unbinned.fna
# Output: $PROJECT/03_unbinned_track/18_dbscan_swa/all/

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/03_unbinned_track/01_raw_fasta/all/unbinned.fna
OUT=$PROJECT/03_unbinned_track/18_dbscan_swa/all
mkdir -p "$OUT"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }

if [ -s "$OUT/bacteria_DBSCAN-SWA_prophage.txt" ] || [ -s "$OUT/_prophage_summary.txt" ]; then
    echo "[$(date '+%F %T')] UB DBSCAN-SWA already done — skip"; exit 0
fi

activate_env "$ENV_DBSCANSWA"

echo "[$(date '+%F %T')] dbscan-swa on UB contigs"
dbscan-swa \
    --input "$IN" \
    --output "$OUT" \
    --prefix unbinned \
    --thread_num "$THREADS" > "$OUT/dbscanswa.log" 2>&1 || \
    { echo "  dbscan-swa returned nonzero — check $OUT/dbscanswa.log"; }

echo "[$(date '+%F %T')] DONE — see $OUT/"
