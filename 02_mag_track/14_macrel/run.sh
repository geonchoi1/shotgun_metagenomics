#!/bin/bash
# === 14 macrel peptides (AMP prediction) on ORF FAA ===
# Input:  $PROJECT/mag/03_master_orf/all/master.faa
# Output: $PROJECT/mag/14_macrel/macrel.prediction (and others)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/mag/03_master_orf/all/master.faa
OUT=$PROJECT/mag/14_macrel
mkdir -p "$OUT"

[ -s "$FAA" ] || { echo "ERROR: $FAA empty" >&2; exit 1; }

activate_env "$ENV_MACREL"

echo "[$(date '+%F %T')] macrel peptides"
macrel peptides \
    --fasta "$FAA" \
    --output "$OUT" \
    --tag macrel \
    --threads "$THREADS" \
    --keep-negatives \
    --force

n=$(zcat -f "$OUT"/macrel.prediction* 2>/dev/null | awk -F'\t' 'NR>1 && $3=="AMP"' | wc -l || echo 0)
echo "[$(date '+%F %T')] DONE — $n predicted AMPs"
