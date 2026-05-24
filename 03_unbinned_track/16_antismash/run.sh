#!/bin/bash
# === 16 antiSMASH (BGCs) on UB contigs (with Bakta gff3) ===
# Input:  $PROJECT/unbinned/02_bakta/unbinned/unbinned.gbff (preferred)
#         (fallback: $PROJECT/unbinned/01_raw_fasta/all/unbinned.fna + Bakta gff3)
# Output: $PROJECT/unbinned/16_antismash/all/<antismash project>/

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

GBFF=$PROJECT/unbinned/02_bakta/unbinned/unbinned.gbff
FNA=$PROJECT/unbinned/01_raw_fasta/all/unbinned.fna
GFF=$PROJECT/unbinned/02_bakta/unbinned/unbinned.gff3
OUT=$PROJECT/unbinned/16_antismash/all
mkdir -p "$OUT"

if [ -s "$OUT/index.html" ]; then
    echo "[$(date '+%F %T')] UB antiSMASH already done — skip"; exit 0
fi

activate_env "$ENV_ANTISMASH"

if [ -s "$GBFF" ]; then
    INPUT="$GBFF"
    EXTRA=""
elif [ -s "$FNA" ] && [ -s "$GFF" ]; then
    INPUT="$FNA"
    EXTRA="--genefinding-gff3 $GFF"
else
    echo "ERROR: need Bakta gbff (preferred) or fna+gff3" >&2; exit 1
fi

echo "[$(date '+%F %T')] antismash on UB ($INPUT)"
antismash \
    --output-dir "$OUT" \
    --cpus "$THREADS" \
    --genefinding-tool none \
    --cb-general --cb-knownclusters --cb-subclusters \
    --asf --pfam2go --smcog-trees \
    $EXTRA \
    "$INPUT" > "$OUT/antismash.log" 2>&1

echo "[$(date '+%F %T')] DONE — see $OUT/index.html"
