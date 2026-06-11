#!/bin/bash
# === dbAPIS (Anti-Phage Immune System) HMM search on master.faa ===
# Tool: HMMER hmmsearch
# DB: dbAPIS (Yan et al. 2024, https://bcb.unl.edu/dbAPIS)
# Output: ORF-level anti-defense system hits
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/01_plasmid_track/04_master_orf/all/master.faa
OUT=$PROJECT/01_plasmid_track/13_dbapis
mkdir -p $OUT

[ -s $OUT/apis.tblout ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_HMMER"

hmmsearch \
  --cpu $THREADS \
  --tblout $OUT/apis.tblout \
  --domtblout $OUT/apis.domtblout \
  -E 1e-5 \
  $DBAPIS_HMM \
  $FAA > $OUT/apis.full.out

echo "[$(date '+%F %T')] dbAPIS hits: $(grep -v '^#' $OUT/apis.tblout | wc -l)"
