#!/bin/bash
# === 07.04 DAS_Tool consensus per sample ===
# Input:  $PROJECT/07_mag/03_binner/{metabinner,metadecoder,semibin,metabat2}/<SAMPLE>/...bins
#         $PROJECT/06_chromosomal/<SAMPLE>/chromosomal.fasta
# Output: $PROJECT/07_mag/04_dastool/<SAMPLE>/_DASTool_bins/

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

REF_BASE=$PROJECT/06_chromosomal
BIN_BASE=$PROJECT/07_mag/03_binner
OUT_BASE=$PROJECT/07_mag/04_dastool
mkdir -p "$OUT_BASE"

activate_env "$ENV_DASTOOL"

# Build a contig2bin TSV from a bins directory ( fasta-extensions vary ).
make_c2b() {
    local bins_dir=$1 binner=$2 out=$3
    : > "$out"
    for f in "$bins_dir"/*.{fa,fasta,fna}; do
        [ -f "$f" ] || continue
        local bin
        bin=$(basename "$f"); bin=${bin%.*}
        grep '^>' "$f" | sed 's/^>//' | awk -v b="$binner.$bin" '{print $1"\t"b}' >> "$out"
    done
}

echo "[$(date '+%F %T')] DAS_Tool consensus"
for d in "$REF_BASE"/*/; do
    [ -d "$d" ] || continue
    sample=$(basename "$d")
    ref="$d/chromosomal.fasta"
    [ -s "$ref" ] || { echo "  skip $sample: no ref"; continue; }

    out_dir="$OUT_BASE/$sample"
    if [ -d "$out_dir/_DASTool_bins" ] && [ "$(ls -A "$out_dir/_DASTool_bins" 2>/dev/null)" ]; then
        echo "  $sample already done"; continue
    fi
    mkdir -p "$out_dir"

    declare -a C2B=() LBL=()
    # Try each binner; only include if bins exist
    for bn in metabinner metadecoder semibin metabat2; do
        case "$bn" in
            metabinner)  bd="$BIN_BASE/metabinner/$sample/metabinner_res"  ;;
            metadecoder) bd="$BIN_BASE/metadecoder/$sample/bins"           ;;
            semibin)     bd="$BIN_BASE/semibin/$sample/output_bins"        ;;
            metabat2)    bd="$BIN_BASE/metabat2/$sample"                   ;;
        esac
        [ -d "$bd" ] || continue
        ls "$bd"/*.{fa,fasta,fna} 2>/dev/null | grep -q . || continue
        c2b="$out_dir/${bn}.c2b.tsv"
        make_c2b "$bd" "$bn" "$c2b"
        [ -s "$c2b" ] || continue
        C2B+=("$c2b"); LBL+=("$bn")
    done

    [ ${#C2B[@]} -ge 1 ] || { echo "  skip $sample: no binner outputs"; continue; }

    echo "  $sample (binners: ${LBL[*]})"
    DAS_Tool \
        -i "$(IFS=, ; echo "${C2B[*]}")" \
        -l "$(IFS=, ; echo "${LBL[*]}")" \
        -c "$ref" \
        -o "$out_dir/" \
        --search_engine diamond \
        --write_bins \
        -t "$THREADS" \
        > "$out_dir/run.log" 2>&1 || true
done

echo "[$(date '+%F %T')] DONE"
