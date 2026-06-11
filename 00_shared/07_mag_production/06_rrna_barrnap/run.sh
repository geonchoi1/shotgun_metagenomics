#!/bin/bash
# === 07.06 barrnap rRNA detection per MAG (bac + arc kingdoms) ===
# Each MAG fasta is scanned with barrnap --kingdom bac and --kingdom arc.
# Output GFFs contain 5S/16S/23S calls.
#
# Input:  $PROJECT/00_shared/07_mag_production/04_dastool/<SAMPLE>/_DASTool_bins/*.fa
# Output: $PROJECT/00_shared/07_mag_production/06_rrna/<MAG>.bac.gff
#         $PROJECT/00_shared/07_mag_production/06_rrna/<MAG>.arc.gff
# (MAG name = <SAMPLE>__<bin>.fa)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

DAS_BASE=$PROJECT/00_shared/07_mag_production/04_dastool
OUT_DIR=$PROJECT/00_shared/07_mag_production/06_rrna
mkdir -p "$OUT_DIR"

activate_env "$ENV_BARRNAP"

echo "[$(date '+%F %T')] barrnap (bac + arc) per MAG"
for d in "$DAS_BASE"/*/_DASTool_bins/; do
    [ -d "$d" ] || continue
    sample=$(basename "$(dirname "$d")")
    for f in "$d"/*.fa; do
        [ -f "$f" ] || continue
        mag="${sample}__$(basename "$f" .fa)"
        bac="$OUT_DIR/${mag}.bac.gff"
        arc="$OUT_DIR/${mag}.arc.gff"
        if [ -s "$bac" ] && [ -s "$arc" ]; then
            continue
        fi
        echo "  $mag"
        [ -s "$bac" ] || barrnap --kingdom bac --threads "$THREADS" "$f" > "$bac" 2>/dev/null || true
        [ -s "$arc" ] || barrnap --kingdom arc --threads "$THREADS" "$f" > "$arc" 2>/dev/null || true
    done
done

echo "[$(date '+%F %T')] DONE"
