#!/bin/bash
# === 19 ISEScan on master FNA (per-MAG concat) ===
# Input:  $PROJECT/02_mag_track/03_master_orf/all/master.fna  (built here from MAG raw fastas if missing)
# Output: $PROJECT/02_mag_track/19_isescan/prediction/<...>.is.fna, *.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

RAW_BASE=$PROJECT/02_mag_track/01_raw_fasta
OUT=$PROJECT/02_mag_track/19_isescan
WORK="$OUT/work"
mkdir -p "$OUT" "$WORK"

CONCAT="$WORK/all_mag.fna"
if [ ! -s "$CONCAT" ]; then
    echo "[$(date '+%F %T')] building per-MAG concat FNA"
    : > "$CONCAT"
    for topo in circ frag; do
        for fa in "$RAW_BASE/$topo"/*.fa; do
            [ -f "$fa" ] || continue
            mag=$(basename "$fa" .fa)
            awk -v m="$mag" '/^>/{sub(/^>/,">"m"|"); print; next}{print}' "$fa" >> "$CONCAT"
        done
    done
fi

activate_env "$ENV_ISESCAN"

echo "[$(date '+%F %T')] isescan.py --removeShortIS"
isescan.py \
    --seqfile "$CONCAT" \
    --output "$OUT" \
    --nthread "$THREADS" \
    --removeShortIS

n=$(find "$OUT" -name '*.tsv' -exec tail -n +2 {} + 2>/dev/null | wc -l)
echo "[$(date '+%F %T')] DONE — $n IS hits"
