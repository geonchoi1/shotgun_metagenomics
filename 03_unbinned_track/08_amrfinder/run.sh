#!/bin/bash
# === 08 AMRFinderPlus on UB master ORFs (protein mode + nucleotide context) ===
# Input:  $PROJECT/unbinned/03_master_orf/all/master.faa
#         $PROJECT/unbinned/03_master_orf/all/master.gff (Bakta gff3)
#         $PROJECT/unbinned/03_master_orf/all/master.fna
# Output: $PROJECT/unbinned/08_amrfinder/all/amrfinder.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

ORF_DIR=$PROJECT/unbinned/03_master_orf/all
OUT=$PROJECT/unbinned/08_amrfinder/all
mkdir -p "$OUT"

[ -s "$ORF_DIR/master.faa" ] || { echo "ERROR: missing master.faa" >&2; exit 1; }
[ -s "$ORF_DIR/master.gff" ] || { echo "ERROR: missing master.gff" >&2; exit 1; }
[ -s "$ORF_DIR/master.fna" ] || { echo "ERROR: missing master.fna" >&2; exit 1; }
[ -d "$AMRFINDER_DB" ] || { echo "ERROR: missing AMRFINDER_DB ($AMRFINDER_DB)" >&2; exit 1; }

if [ -s "$OUT/amrfinder.tsv" ]; then
    echo "[$(date '+%F %T')] UB AMRFinder already done — skip"; exit 0
fi

activate_env "$ENV_AMRFINDER"

echo "[$(date '+%F %T')] amrfinder --protein + --nucleotide + --gff"
amrfinder \
    --protein    "$ORF_DIR/master.faa" \
    --nucleotide "$ORF_DIR/master.fna" \
    --gff        "$ORF_DIR/master.gff" \
    --annotation_format bakta \
    --database   "$AMRFINDER_DB" \
    --threads    "$THREADS" \
    --plus \
    -o "$OUT/amrfinder.tsv" > "$OUT/amrfinder.log" 2>&1

n=$(awk 'NR>1' "$OUT/amrfinder.tsv" | wc -l)
echo "[$(date '+%F %T')] DONE — AMR/stress/virulence hits: $n"
