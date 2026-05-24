#!/bin/bash
# === 06 KofamScan (KO) on UB master ORFs ===
# Input:  $PROJECT/unbinned/03_master_orf/all/master.faa
# Output: $PROJECT/unbinned/06_kofamscan/all/kofam.tsv  (detail-tsv)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/unbinned/03_master_orf/all/master.faa
OUT=$PROJECT/unbinned/06_kofamscan/all
TMP=$OUT/tmp
mkdir -p "$OUT" "$TMP"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }
[ -s "$KOFAM_KO_LIST" ] || { echo "ERROR: missing $KOFAM_KO_LIST" >&2; exit 1; }
[ -d "$KOFAM_PROFILES" ] || { echo "ERROR: missing $KOFAM_PROFILES" >&2; exit 1; }

if [ -s "$OUT/kofam.tsv" ]; then
    echo "[$(date '+%F %T')] UB KofamScan already done — skip"; exit 0
fi

activate_env "$ENV_KOFAMSCAN"

echo "[$(date '+%F %T')] KofamScan on UB ORFs"
exec_annotation \
    -p "$KOFAM_PROFILES" \
    -k "$KOFAM_KO_LIST" \
    --cpu "$THREADS" \
    --tmp-dir "$TMP" \
    -f detail-tsv \
    -o "$OUT/kofam.tsv" \
    "$IN"

rm -rf "$TMP"
n=$(awk -F'\t' 'NR>1 && $1=="*"' "$OUT/kofam.tsv" | wc -l)
echo "[$(date '+%F %T')] DONE — KO assignments (* significant): $n"
