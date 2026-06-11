#!/bin/bash
# === 04 Pfam HMM search (--cut_ga) ===
# Input:  $PROJECT/02_mag_track/03_master_orf/all/master.faa
# Output: $PROJECT/02_mag_track/04_pfam/pfam.tblout, pfam.domtblout

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/02_mag_track/03_master_orf/all/master.faa
OUT=$PROJECT/02_mag_track/04_pfam
mkdir -p "$OUT"

[ -s "$FAA" ] || { echo "ERROR: $FAA empty" >&2; exit 1; }

activate_env "$ENV_HMMER"

if [ ! -f "${PFAM_HMM}.h3i" ]; then
    hmmpress "$PFAM_HMM"
fi

echo "[$(date '+%F %T')] hmmsearch --cut_ga vs Pfam-A"
hmmsearch --cpu "$THREADS_HMMSEARCH" --cut_ga \
    --tblout "$OUT/pfam.tblout" \
    --domtblout "$OUT/pfam.domtblout" \
    "$PFAM_HMM" "$FAA" > "$OUT/pfam.stdout"

n=$(grep -vc '^#' "$OUT/pfam.tblout" || true)
echo "[$(date '+%F %T')] DONE — $n hits"
