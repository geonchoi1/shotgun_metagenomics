#!/bin/bash
# === 05 NCBIfam (PGAP) HMM search (--cut_ga) ===
# Input:  $PROJECT/mag/03_master_orf/all/master.faa
# Output: $PROJECT/mag/05_ncbifam/ncbifam.tblout, ncbifam.domtblout

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/mag/03_master_orf/all/master.faa
OUT=$PROJECT/mag/05_ncbifam
mkdir -p "$OUT"

[ -s "$FAA" ] || { echo "ERROR: $FAA empty" >&2; exit 1; }

activate_env "$ENV_HMMER"

if [ ! -f "${NCBIFAM_HMM}.h3i" ]; then
    hmmpress "$NCBIFAM_HMM"
fi

echo "[$(date '+%F %T')] hmmsearch --cut_ga vs NCBIfam"
hmmsearch --cpu "$THREADS_HMMSEARCH" --cut_ga \
    --tblout "$OUT/ncbifam.tblout" \
    --domtblout "$OUT/ncbifam.domtblout" \
    "$NCBIFAM_HMM" "$FAA" > "$OUT/ncbifam.stdout"

n=$(grep -vc '^#' "$OUT/ncbifam.tblout" || true)
echo "[$(date '+%F %T')] DONE — $n hits"
