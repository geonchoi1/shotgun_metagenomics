#!/bin/bash
# === 12 dbAPIS (anti-phage immune system) hmmsearch on UB ORFs ===
# Input:  $PROJECT/unbinned/03_master_orf/all/master.faa
# Output: $PROJECT/unbinned/12_dbapis/all/dbapis.{tblout,domtblout}

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/unbinned/03_master_orf/all/master.faa
OUT=$PROJECT/unbinned/12_dbapis/all
mkdir -p "$OUT"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }
[ -s "$DBAPIS_HMM" ] || { echo "ERROR: missing $DBAPIS_HMM" >&2; exit 1; }

if [ -s "$OUT/dbapis.tblout" ]; then
    echo "[$(date '+%F %T')] UB dbAPIS already done — skip"; exit 0
fi

activate_env "$ENV_HMMER"

echo "[$(date '+%F %T')] hmmsearch dbAPIS on UB ORFs"
hmmsearch --cpu "$THREADS_HMMSEARCH" \
          -E 1e-10 \
          --tblout    "$OUT/dbapis.tblout" \
          --domtblout "$OUT/dbapis.domtblout" \
          -o          "$OUT/dbapis.log" \
          "$DBAPIS_HMM" "$IN"

n=$(grep -vc '^#' "$OUT/dbapis.tblout" || true)
echo "[$(date '+%F %T')] DONE — dbAPIS hits: $n"
