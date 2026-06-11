#!/bin/bash
# === 04 Pfam-A hmmsearch on UB master ORFs ===
# Input:  $PROJECT/03_unbinned_track/03_master_orf/all/master.faa
# Output: $PROJECT/03_unbinned_track/04_pfam/all/pfam.tblout, pfam.domtblout

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/03_unbinned_track/03_master_orf/all/master.faa
OUT=$PROJECT/03_unbinned_track/04_pfam/all
mkdir -p "$OUT"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }
[ -s "$PFAM_HMM" ] || { echo "ERROR: missing $PFAM_HMM" >&2; exit 1; }

if [ -s "$OUT/pfam.tblout" ]; then
    echo "[$(date '+%F %T')] UB Pfam already done — skip"; exit 0
fi

activate_env "$ENV_HMMER"

echo "[$(date '+%F %T')] hmmsearch Pfam-A on UB ORFs"
hmmsearch --cpu "$THREADS_HMMSEARCH" \
          --cut_ga \
          --tblout    "$OUT/pfam.tblout" \
          --domtblout "$OUT/pfam.domtblout" \
          -o          "$OUT/pfam.log" \
          "$PFAM_HMM" "$IN"

n=$(grep -vc '^#' "$OUT/pfam.tblout" || true)
echo "[$(date '+%F %T')] DONE — Pfam hits: $n"
