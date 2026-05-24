#!/bin/bash
# === 06 KofamScan (exec_annotation -f mapper) ===
# Input:  $PROJECT/mag/03_master_orf/all/master.faa
# Output: $PROJECT/mag/06_kofamscan/kofam_mapper.tsv  (locus_tag<TAB>KO ...)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/mag/03_master_orf/all/master.faa
OUT=$PROJECT/mag/06_kofamscan
mkdir -p "$OUT" "$OUT/tmp"

[ -s "$FAA" ] || { echo "ERROR: $FAA empty" >&2; exit 1; }

activate_env "$ENV_KOFAMSCAN"

echo "[$(date '+%F %T')] exec_annotation -f mapper"
exec_annotation \
    -p "$KOFAM_PROFILES" \
    -k "$KOFAM_KO_LIST" \
    -f mapper \
    --cpu "$THREADS" \
    --tmp-dir "$OUT/tmp" \
    -o "$OUT/kofam_mapper.tsv" \
    "$FAA"

n=$(awk -F'\t' 'NF>=2 && $2!=""' "$OUT/kofam_mapper.tsv" | wc -l)
echo "[$(date '+%F %T')] DONE — $n ORFs with KO"
rm -rf "$OUT/tmp"
