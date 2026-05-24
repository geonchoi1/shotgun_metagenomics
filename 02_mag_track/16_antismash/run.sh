#!/bin/bash
# === 16 antiSMASH per MAG (FNA) ===
# Input:  $PROJECT/mag/01_raw_fasta/{circ,frag}/*.fa
# Output: $PROJECT/mag/16_antismash/<MAG>/index.html, *.gbk, regions.js

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN_BASE=$PROJECT/mag/01_raw_fasta
OUT_BASE=$PROJECT/mag/16_antismash
mkdir -p "$OUT_BASE"

activate_env "$ENV_ANTISMASH"

echo "[$(date '+%F %T')] antiSMASH per MAG"
for topo in circ frag; do
    for fa in "$IN_BASE/$topo"/*.fa; do
        [ -f "$fa" ] || continue
        mag=$(basename "$fa" .fa)
        out="$OUT_BASE/$mag"
        if [ -s "$out/index.html" ]; then
            echo "  $mag already done"
            continue
        fi
        mkdir -p "$out"
        extra=""
        [ "$topo" = "frag" ] && extra="--allow-long-headers"
        echo "  $topo/$mag"
        antismash \
            --cpus "$THREADS" \
            --output-dir "$out" \
            --output-basename "$mag" \
            --genefinding-tool prodigal-m \
            --cb-general --cb-knownclusters --cb-subclusters \
            --asf --pfam2go --rre --tigrfam \
            "$fa" > "$out/antismash.log" 2>&1 || { echo "  FAIL $mag — see $out/antismash.log"; continue; }
    done
done

n=$(find "$OUT_BASE" -maxdepth 2 -name 'index.html' | wc -l)
echo "[$(date '+%F %T')] DONE — $n MAGs annotated"
