#!/bin/bash
# === dbAPIS (anti-phage immune systems) HMM search ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/plasmid/04_master_orf/all/master.faa
OUT=$PROJECT/plasmid/13_dbapis
mkdir -p $OUT

[ -s $OUT/dbapis.tblout ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_HMMER"

for ext in h3i h3f h3m h3p; do
  [ -f ${DBAPIS_HMM}.${ext} ] || { hmmpress -f $DBAPIS_HMM; break; }
done

hmmsearch -E 1e-5 --cpu $THREADS_HMMSEARCH \
  --tblout $OUT/dbapis.tblout \
  --domtblout $OUT/dbapis.domtblout \
  $DBAPIS_HMM $FAA > $OUT/dbapis.stdout

echo "[$(date '+%F %T')] dbAPIS hits: $(grep -vc '^#' $OUT/dbapis.tblout)"
