#!/bin/bash
# === 20 IntegronFinder on UB contigs (--local-max --func-annot) ===
# Input:  $PROJECT/unbinned/01_raw_fasta/all/unbinned.fna
# Output: $PROJECT/unbinned/20_integronfinder/all/Results_Integron_Finder_unbinned/

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/unbinned/01_raw_fasta/all/unbinned.fna
OUT=$PROJECT/unbinned/20_integronfinder/all
mkdir -p "$OUT"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }

if [ -d "$OUT/Results_Integron_Finder_unbinned" ]; then
    echo "[$(date '+%F %T')] UB IntegronFinder already done — skip"; exit 0
fi

activate_env "$ENV_INTEGRONFINDER"

echo "[$(date '+%F %T')] integron_finder on UB contigs"
integron_finder \
    --local-max \
    --func-annot \
    --cpu "$THREADS" \
    --outdir "$OUT" \
    "$IN" > "$OUT/integronfinder.log" 2>&1

summ=$(find "$OUT/Results_Integron_Finder_unbinned" -name "*.summary" | head -1)
[ -n "$summ" ] && n=$(awk 'NR>1' "$summ" | wc -l) || n=0
echo "[$(date '+%F %T')] DONE — integron summary rows: $n"
