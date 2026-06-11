#!/bin/bash
# === 07 eggNOG-mapper on UB master ORFs ===
# Input:  $PROJECT/03_unbinned_track/03_master_orf/all/master.faa
# Output: $PROJECT/03_unbinned_track/07_eggnog/all/eggnog.emapper.{annotations,seed_orthologs,hits}

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/03_unbinned_track/03_master_orf/all/master.faa
OUT=$PROJECT/03_unbinned_track/07_eggnog/all
TMP=$OUT/tmp
mkdir -p "$OUT" "$TMP"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }
[ -d "$EGGNOG_DATA" ] || { echo "ERROR: missing $EGGNOG_DATA" >&2; exit 1; }

if [ -s "$OUT/eggnog.emapper.annotations" ]; then
    echo "[$(date '+%F %T')] UB eggNOG already done — skip"; exit 0
fi

activate_env "$ENV_EGGNOG"

echo "[$(date '+%F %T')] emapper.py on UB ORFs"
emapper.py \
    -i "$IN" \
    --itype proteins \
    --output eggnog \
    --output_dir "$OUT" \
    --temp_dir "$TMP" \
    --data_dir "$EGGNOG_DATA" \
    --cpu "$THREADS" \
    --override

rm -rf "$TMP"
n=$(grep -vc '^#' "$OUT/eggnog.emapper.annotations" || true)
echo "[$(date '+%F %T')] DONE — eggNOG annotations: $n"
