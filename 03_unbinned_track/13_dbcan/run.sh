#!/bin/bash
# === 13 dbCAN (CAZymes) on UB master ORFs ===
# Input:  $PROJECT/unbinned/03_master_orf/all/master.faa
# Output: $PROJECT/unbinned/13_dbcan/all/overview.txt + hmmer/dimaond/dbsub results

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/unbinned/03_master_orf/all/master.faa
OUT=$PROJECT/unbinned/13_dbcan/all
mkdir -p "$OUT"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }
[ -d "$DBCAN_DB" ] || { echo "ERROR: missing $DBCAN_DB" >&2; exit 1; }

if [ -s "$OUT/overview.txt" ]; then
    echo "[$(date '+%F %T')] UB dbCAN already done — skip"; exit 0
fi

activate_env "$ENV_DBCAN"

echo "[$(date '+%F %T')] run_dbcan on UB ORFs"
run_dbcan "$IN" protein \
    --db_dir "$DBCAN_DB" \
    --out_dir "$OUT" \
    --dia_cpu "$THREADS" \
    --hmm_cpu "$THREADS" \
    --tf_cpu  "$THREADS" \
    --stp_cpu "$THREADS" \
    > "$OUT/dbcan.log" 2>&1

n=$(awk 'NR>1' "$OUT/overview.txt" | wc -l)
echo "[$(date '+%F %T')] DONE — dbCAN overview rows: $n"
