#!/bin/bash
# === Pfam-A HMM search on master.faa ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/01_plasmid_track/04_master_orf/all/master.faa
OUT=$PROJECT/01_plasmid_track/05_pfam
mkdir -p $OUT

[ -s $OUT/pfam.tblout ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_HMMER"

# hmmpress if needed
for ext in h3i h3f h3m h3p; do
  [ -f ${PFAM_HMM}.${ext} ] || { hmmpress -f $PFAM_HMM; break; }
done

hmmsearch --cut_ga --cpu $THREADS_HMMSEARCH \
  --tblout $OUT/pfam.tblout \
  --domtblout $OUT/pfam.domtblout \
  $PFAM_HMM $FAA > $OUT/pfam.stdout

echo "[$(date '+%F %T')] Pfam hits: $(grep -vc '^#' $OUT/pfam.tblout)"
