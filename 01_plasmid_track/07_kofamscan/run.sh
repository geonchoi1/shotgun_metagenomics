#!/bin/bash
# === KofamScan: KEGG Orthology HMM assignment ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/plasmid/04_master_orf/all/master.faa
OUT=$PROJECT/plasmid/07_kofamscan
mkdir -p $OUT $OUT/tmp

[ -s $OUT/kofam.mapper.tsv ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_KOFAMSCAN"

exec_annotation -p $KOFAM_PROFILES -k $KOFAM_KO_LIST \
  -f mapper \
  --cpu $THREADS \
  --tmp-dir $OUT/tmp \
  -o $OUT/kofam.mapper.tsv \
  $FAA

# Also detail format for downstream scoring
exec_annotation -p $KOFAM_PROFILES -k $KOFAM_KO_LIST \
  -f detail-tsv --no-report-unannotated \
  --cpu $THREADS \
  --tmp-dir $OUT/tmp \
  -o $OUT/kofam.detail.tsv \
  $FAA

echo "[$(date '+%F %T')] KO assigned: $(awk -F'\t' '$2!=""' $OUT/kofam.mapper.tsv | wc -l)"
