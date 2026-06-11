#!/bin/bash
# === 13 dbCAN CAZyme annotation (protein mode) ===
# Input:  $PROJECT/02_mag_track/03_master_orf/all/master.faa
# Output: $PROJECT/02_mag_track/13_dbcan/overview.tsv (and dbCAN files)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/02_mag_track/03_master_orf/all/master.faa
OUT=$PROJECT/02_mag_track/13_dbcan
mkdir -p "$OUT"

[ -s "$FAA" ] || { echo "ERROR: $FAA empty" >&2; exit 1; }

activate_env "$ENV_DBCAN"

echo "[$(date '+%F %T')] run_dbcan CAZyme_annotation --mode protein"
run_dbcan CAZyme_annotation \
    --mode protein \
    --input_raw_data "$FAA" \
    --db_dir "$DBCAN_DB" \
    --threads "$THREADS" \
    --output_dir "$OUT"

echo "[$(date '+%F %T')] DONE"
