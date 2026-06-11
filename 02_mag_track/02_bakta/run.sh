#!/bin/bash
# === 02 Bakta per-MAG annotation (circ uses --complete) ===
# Input:  $PROJECT/02_mag_track/01_raw_fasta/{circ,frag}/*.fa
# Output: $PROJECT/02_mag_track/02_bakta/{circ,frag}/<MAG>/<MAG>.{faa,ffn,fna,gff3,tsv}

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN_BASE=$PROJECT/02_mag_track/01_raw_fasta
OUT_BASE=$PROJECT/02_mag_track/02_bakta
mkdir -p "$OUT_BASE/circ" "$OUT_BASE/frag"

activate_env "$ENV_BAKTA"

run_bakta() {
    local topo=$1
    local extra=""
    [ "$topo" = "circ" ] && extra="--complete"
    for fa in "$IN_BASE/$topo"/*.fa; do
        [ -f "$fa" ] || continue
        mag=$(basename "$fa" .fa)
        out="$OUT_BASE/$topo/$mag"
        if [ -s "$out/${mag}.gff3" ]; then
            echo "  $topo/$mag already done"
            continue
        fi
        mkdir -p "$out"
        echo "  bakta $topo/$mag"
        bakta --db "$BAKTA_DB" \
              --output "$out" \
              --prefix "$mag" \
              --threads "$THREADS_BAKTA" \
              --locus-tag "$mag" \
              --force \
              $extra \
              "$fa" > "$out/bakta.log" 2>&1 || { echo "  FAIL $topo/$mag — see $out/bakta.log"; continue; }
    done
}

echo "[$(date '+%F %T')] Bakta on circ MAGs (--complete)"
run_bakta circ
echo "[$(date '+%F %T')] Bakta on frag MAGs"
run_bakta frag
echo "[$(date '+%F %T')] DONE"
