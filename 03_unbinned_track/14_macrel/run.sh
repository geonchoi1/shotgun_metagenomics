#!/bin/bash
# === 14 Macrel (antimicrobial peptides) on UB contigs ===
# Input:  $PROJECT/unbinned/01_raw_fasta/all/unbinned.fna
# Output: $PROJECT/unbinned/14_macrel/all/macrel.prediction

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/unbinned/01_raw_fasta/all/unbinned.fna
OUT=$PROJECT/unbinned/14_macrel/all
mkdir -p "$OUT"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }

if [ -s "$OUT/macrel.prediction" ] || [ -s "$OUT/macrel.smorfs.faa.gz" ]; then
    echo "[$(date '+%F %T')] UB Macrel already done — skip"; exit 0
fi

activate_env "$ENV_MACREL"

echo "[$(date '+%F %T')] macrel contigs on UB"
macrel contigs \
    --fasta "$IN" \
    --output "$OUT" \
    --tag macrel \
    --threads "$THREADS" \
    --force > "$OUT/macrel.log" 2>&1

n=$(zcat -f "$OUT/macrel.prediction" 2>/dev/null | awk 'NR>1' | wc -l || true)
echo "[$(date '+%F %T')] DONE — Macrel AMP predictions: $n"
