#!/bin/bash
# === 19 ISEScan (insertion sequences) on UB contigs ===
# Input:  $PROJECT/03_unbinned_track/01_raw_fasta/all/unbinned.fna
# Output: $PROJECT/03_unbinned_track/19_isescan/all/unbinned.fna.{tsv,sum,gff,is.fna}

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/03_unbinned_track/01_raw_fasta/all/unbinned.fna
OUT=$PROJECT/03_unbinned_track/19_isescan/all
mkdir -p "$OUT"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }

if ls "$OUT"/*.tsv >/dev/null 2>&1; then
    echo "[$(date '+%F %T')] UB ISEScan already done — skip"; exit 0
fi

activate_env "$ENV_ISESCAN"

echo "[$(date '+%F %T')] isescan on UB contigs"
( cd "$OUT" && isescan.py \
    --seqfile "$IN" \
    --output "$OUT" \
    --nthread "$THREADS" \
    > "$OUT/isescan.log" 2>&1 )

n=$(awk 'NR>1' "$OUT"/*.tsv 2>/dev/null | wc -l || true)
echo "[$(date '+%F %T')] DONE — IS predictions: $n"
