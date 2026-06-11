#!/bin/bash
# === 07 eggNOG-mapper ===
# Input:  $PROJECT/02_mag_track/03_master_orf/all/master.faa
# Output: $PROJECT/02_mag_track/07_eggnog/eggnog.emapper.annotations

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/02_mag_track/03_master_orf/all/master.faa
OUT=$PROJECT/02_mag_track/07_eggnog
mkdir -p "$OUT"

[ -s "$FAA" ] || { echo "ERROR: $FAA empty" >&2; exit 1; }

activate_env "$ENV_EGGNOG"

echo "[$(date '+%F %T')] emapper.py"
emapper.py \
    -i "$FAA" \
    --itype proteins \
    --data_dir "$EGGNOG_DATA" \
    --output eggnog \
    --output_dir "$OUT" \
    --cpu "$THREADS" \
    --override

n=$(grep -vc '^#' "$OUT/eggnog.emapper.annotations" 2>/dev/null || echo 0)
echo "[$(date '+%F %T')] DONE — $n annotated ORFs"
