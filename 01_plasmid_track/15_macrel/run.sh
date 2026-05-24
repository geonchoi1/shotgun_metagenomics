#!/bin/bash
# === Macrel: AMP (antimicrobial peptide) prediction on master.faa ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/plasmid/04_master_orf/all/master.faa
OUT=$PROJECT/plasmid/15_macrel
mkdir -p $OUT

[ -s $OUT/macrel.prediction.gz ] || [ -s $OUT/macrel.prediction ] && { echo "skip (exists)"; exit 0; } || true

activate_env "$ENV_MACREL"

macrel peptides \
  --fasta $FAA \
  --output $OUT \
  --tag macrel \
  --threads $THREADS \
  --force

echo "[$(date '+%F %T')] Macrel AMP candidates: $(zcat -f $OUT/macrel.prediction.gz 2>/dev/null | grep -vc '^#' || echo NA)"
