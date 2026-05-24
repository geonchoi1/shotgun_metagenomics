#!/bin/bash
# === 20 IntegronFinder2 --local-max on per-MAG concat FNA ===
# Input:  $PROJECT/mag/19_isescan/work/all_mag.fna  (reuse), else builds it
# Output: $PROJECT/mag/20_integronfinder/Results_Integron_Finder_all_mag/

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

RAW_BASE=$PROJECT/mag/01_raw_fasta
OUT=$PROJECT/mag/20_integronfinder
mkdir -p "$OUT"

CONCAT=$PROJECT/mag/19_isescan/work/all_mag.fna
if [ ! -s "$CONCAT" ]; then
    echo "[$(date '+%F %T')] building per-MAG concat FNA (19_isescan not yet run)"
    mkdir -p "$(dirname "$CONCAT")"
    : > "$CONCAT"
    for topo in circ frag; do
        for fa in "$RAW_BASE/$topo"/*.fa; do
            [ -f "$fa" ] || continue
            mag=$(basename "$fa" .fa)
            awk -v m="$mag" '/^>/{sub(/^>/,">"m"|"); print; next}{print}' "$fa" >> "$CONCAT"
        done
    done
fi

activate_env "$ENV_INTEGRONFINDER"

echo "[$(date '+%F %T')] integron_finder --local-max"
integron_finder \
    --local-max \
    --cpu "$THREADS" \
    --outdir "$OUT" \
    --func-annot \
    --gbk \
    "$CONCAT"

result_dir=$(find "$OUT" -maxdepth 1 -name 'Results_Integron_Finder_*' -type d | head -1)
n=$(awk 'NR>1' "$result_dir"/*.integrons 2>/dev/null | wc -l || echo 0)
echo "[$(date '+%F %T')] DONE — $n integron rows"
