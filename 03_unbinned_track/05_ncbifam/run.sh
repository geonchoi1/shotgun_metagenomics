#!/bin/bash
# === 05 NCBIfam (PGAP) hmmsearch on UB master ORFs ===
# Input:  $PROJECT/03_unbinned_track/03_master_orf/all/master.faa
# Output: $PROJECT/03_unbinned_track/05_ncbifam/all/ncbifam.{tblout,domtblout}

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/03_unbinned_track/03_master_orf/all/master.faa
OUT=$PROJECT/03_unbinned_track/05_ncbifam/all
mkdir -p "$OUT"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }
[ -s "$NCBIFAM_HMM" ] || { echo "ERROR: missing $NCBIFAM_HMM" >&2; exit 1; }

if [ -s "$OUT/ncbifam.tblout" ]; then
    echo "[$(date '+%F %T')] UB NCBIfam already done — skip"; exit 0
fi

activate_env "$ENV_HMMER"

echo "[$(date '+%F %T')] hmmsearch NCBIfam on UB ORFs"
hmmsearch --cpu "$THREADS_HMMSEARCH" \
          --cut_tc \
          --tblout    "$OUT/ncbifam.tblout" \
          --domtblout "$OUT/ncbifam.domtblout" \
          -o          "$OUT/ncbifam.log" \
          "$NCBIFAM_HMM" "$IN"

n=$(grep -vc '^#' "$OUT/ncbifam.tblout" || true)
echo "[$(date '+%F %T')] DONE — NCBIfam hits: $n"
