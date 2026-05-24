#!/bin/bash
# === 18 DBSCAN-SWA (prophage detection) per MAG ===
# Input:  $PROJECT/mag/01_raw_fasta/{circ,frag}/*.fa
# Output: $PROJECT/mag/18_dbscan_swa/<MAG>/

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN_BASE=$PROJECT/mag/01_raw_fasta
OUT_BASE=$PROJECT/mag/18_dbscan_swa
mkdir -p "$OUT_BASE"

activate_env "$ENV_DBSCANSWA"

# DBSCAN-SWA install path (override DBSCAN_SWA_DIR if needed)
DBSCAN_SWA_DIR=${DBSCAN_SWA_DIR:-$HOME/tools/DBSCAN-SWA}
SCRIPT="$DBSCAN_SWA_DIR/bin/dbscan-swa.py"
[ -f "$SCRIPT" ] || { echo "ERROR: dbscan-swa.py not found at $SCRIPT (set DBSCAN_SWA_DIR)" >&2; exit 1; }

echo "[$(date '+%F %T')] DBSCAN-SWA per MAG"
for topo in circ frag; do
    for fa in "$IN_BASE/$topo"/*.fa; do
        [ -f "$fa" ] || continue
        mag=$(basename "$fa" .fa)
        out="$OUT_BASE/$mag"
        if [ -s "$out/bacteria_DBSCAN-SWA_prophage_summary.txt" ]; then
            echo "  $mag already done"
            continue
        fi
        mkdir -p "$out"
        echo "  $topo/$mag"
        python "$SCRIPT" \
            --input "$fa" \
            --output "$out" \
            --prefix "$mag" \
            --thread_num "$THREADS" \
            > "$out/dbscanswa.log" 2>&1 || { echo "  FAIL $mag — see $out/dbscanswa.log"; continue; }
    done
done

n=$(find "$OUT_BASE" -maxdepth 2 -name '*prophage_summary.txt' | wc -l)
echo "[$(date '+%F %T')] DONE — $n MAGs scanned"
