#!/bin/bash
# === 10 VFDB virulence factor on MAG ===
# DIAMOND blastp vs VFDB setB (extended) — pathogen / virulence gene detection
#
# Input:  $PROJECT/02_mag_track/03_master_orf/all/master.faa
# Output: $PROJECT/02_mag_track/10_vfdb/vfdb.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

MASTER=$PROJECT/02_mag_track/03_master_orf/all/master.faa
OUT=$PROJECT/02_mag_track/10_vfdb
mkdir -p "$OUT"

[ -s "$MASTER" ] || { echo "ERROR: missing $MASTER" >&2; exit 1; }
[ -s "$OUT/vfdb.tsv" ] && { echo "[$(date '+%F %T')] VFDB already done — skip"; exit 0; }

activate_env "$ENV_DIAMOND"

echo "[$(date '+%F %T')] DIAMOND blastp vs VFDB setB"
diamond blastp \
    -q "$MASTER" \
    -d "$VFDB_DMND" \
    --id 70 --query-cover 50 --subject-cover 50 -e 1e-5 \
    --outfmt 6 \
    -p "$THREADS" \
    -o "$OUT/vfdb.tsv" 2>&1 | tail -5

echo "[$(date '+%F %T')] DONE — hits: $(wc -l < $OUT/vfdb.tsv)"
