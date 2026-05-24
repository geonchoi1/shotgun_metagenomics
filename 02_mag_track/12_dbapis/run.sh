#!/bin/bash
# === 12 dbAPIS (anti-phage immune systems) HMM search ===
# Input:  $PROJECT/mag/03_master_orf/all/master.faa
# Output: $PROJECT/mag/12_dbapis/dbapis.tblout

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/mag/03_master_orf/all/master.faa
OUT=$PROJECT/mag/12_dbapis
mkdir -p "$OUT"

[ -s "$FAA" ] || { echo "ERROR: $FAA empty" >&2; exit 1; }

activate_env "$ENV_HMMER"

if [ ! -f "${DBAPIS_HMM}.h3i" ]; then
    hmmpress "$DBAPIS_HMM"
fi

echo "[$(date '+%F %T')] hmmsearch -E 1e-5 vs dbAPIS"
hmmsearch --cpu "$THREADS_HMMSEARCH" -E 1e-5 \
    --tblout "$OUT/dbapis.tblout" \
    --domtblout "$OUT/dbapis.domtblout" \
    "$DBAPIS_HMM" "$FAA" > "$OUT/dbapis.stdout"

n=$(grep -vc '^#' "$OUT/dbapis.tblout" || true)
echo "[$(date '+%F %T')] DONE — $n hits"
