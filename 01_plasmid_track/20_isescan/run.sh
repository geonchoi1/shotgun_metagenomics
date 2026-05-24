#!/bin/bash
# === ISEScan: insertion sequence detection on plasmids ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FNA=$PROJECT/plasmid/02_drep/dereplicated.fna
OUT=$PROJECT/plasmid/20_isescan
mkdir -p $OUT

[ -s $OUT/results/$(basename $FNA).tsv ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_ISESCAN"

cd $OUT
isescan.py \
  --removeShortIS \
  --seqfile $FNA \
  --output $OUT/results \
  --nthread $THREADS

echo "[$(date '+%F %T')] ISEScan done"
ls $OUT/results/ | head -5
