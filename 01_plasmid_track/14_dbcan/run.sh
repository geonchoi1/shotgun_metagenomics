#!/bin/bash
# === dbCAN: CAZyme annotation (diamond + HMM) ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/01_plasmid_track/04_master_orf/all/master.faa
OUT=$PROJECT/01_plasmid_track/14_dbcan
mkdir -p $OUT

[ -d $OUT/overview.tsv ] && { echo "skip (exists)"; exit 0; }
[ -s $OUT/overview.tsv ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_DBCAN"

run_dbcan CAZyme_annotation \
  --input_raw_data $FAA \
  --mode protein \
  --methods diamond hmm \
  --db_dir $DBCAN_DB \
  --output_dir $OUT \
  --threads $THREADS

echo "[$(date '+%F %T')] dbCAN CAZyme calls: $(wc -l < $OUT/overview.tsv 2>/dev/null || echo NA)"
