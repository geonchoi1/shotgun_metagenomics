#!/bin/bash
# === 07.07 tRNAscan-SE per MAG (-B bac; -A arc optional) ===
# Input:  $PROJECT/07_mag/04_dastool/<SAMPLE>/_DASTool_bins/*.fa
# Output: $PROJECT/07_mag/07_trna/<MAG>.bac.tsv
#         $PROJECT/07_mag/07_trna/<MAG>.arc.tsv (if RUN_ARC=1)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

DAS_BASE=$PROJECT/07_mag/04_dastool
OUT_DIR=$PROJECT/07_mag/07_trna
mkdir -p "$OUT_DIR"

RUN_ARC=${RUN_ARC:-1}

activate_env "$ENV_TRNASCAN"

echo "[$(date '+%F %T')] tRNAscan-SE per MAG"
for d in "$DAS_BASE"/*/_DASTool_bins/; do
    [ -d "$d" ] || continue
    sample=$(basename "$(dirname "$d")")
    for f in "$d"/*.fa; do
        [ -f "$f" ] || continue
        mag="${sample}__$(basename "$f" .fa)"
        bac="$OUT_DIR/${mag}.bac.tsv"
        arc="$OUT_DIR/${mag}.arc.tsv"
        if [ -s "$bac" ] && { [ "$RUN_ARC" = "0" ] || [ -s "$arc" ]; }; then
            continue
        fi
        echo "  $mag"
        [ -s "$bac" ] || tRNAscan-SE -B -q --thread "$THREADS" -o "$bac" "$f" \
            > "$OUT_DIR/${mag}.bac.log" 2>&1 || true
        if [ "$RUN_ARC" = "1" ] && [ ! -s "$arc" ]; then
            tRNAscan-SE -A -q --thread "$THREADS" -o "$arc" "$f" \
                > "$OUT_DIR/${mag}.arc.log" 2>&1 || true
        fi
    done
done

echo "[$(date '+%F %T')] DONE"
