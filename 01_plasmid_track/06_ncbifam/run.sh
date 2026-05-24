#!/bin/bash
# === NCBIfam HMM search (uses Pfam-style trusted cutoffs) ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/plasmid/04_master_orf/all/master.faa
OUT=$PROJECT/plasmid/06_ncbifam
mkdir -p $OUT

[ -s $OUT/ncbifam.tblout ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_HMMER"

for ext in h3i h3f h3m h3p; do
  [ -f ${NCBIFAM_HMM}.${ext} ] || { hmmpress -f $NCBIFAM_HMM; break; }
done

hmmsearch --cut_ga --cpu $THREADS_HMMSEARCH \
  --tblout $OUT/ncbifam.tblout \
  --domtblout $OUT/ncbifam.domtblout \
  $NCBIFAM_HMM $FAA > $OUT/ncbifam.stdout

echo "[$(date '+%F %T')] NCBIfam hits: $(grep -vc '^#' $OUT/ncbifam.tblout)"
