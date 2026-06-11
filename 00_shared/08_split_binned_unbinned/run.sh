#!/bin/bash
# === 08 Split chromosomal contigs into binned vs unbinned ===
# Binned   = chromosomal contig captured in ANY DAS_Tool MAG (any sample)
# Unbinned = the rest
#
# Input:  $PROJECT/00_shared/06_chromosomal_extract/<SAMPLE>/chromosomal.fasta
#         $PROJECT/00_shared/07_mag_production/04_dastool/<SAMPLE>/_DASTool_bins/*.fa
# Output: $PROJECT/00_shared/08_split_binned_unbinned/binned/all.fna
#         $PROJECT/00_shared/08_split_binned_unbinned/unbinned/all.fna

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

CHR_BASE=$PROJECT/00_shared/06_chromosomal_extract
DAS_BASE=$PROJECT/00_shared/07_mag_production/04_dastool
OUT_BASE=$PROJECT/00_shared/08_split_binned_unbinned
mkdir -p "$OUT_BASE/binned" "$OUT_BASE/unbinned"

BINNED_FA=$OUT_BASE/binned/all.fna
UNBINNED_FA=$OUT_BASE/unbinned/all.fna

if [ -s "$BINNED_FA" ] && [ -s "$UNBINNED_FA" ]; then
    echo "[$(date '+%F %T')] 08 split already done"; exit 0
fi

: > "$BINNED_FA"; : > "$UNBINNED_FA"

echo "[$(date '+%F %T')] split per sample"
for d in "$CHR_BASE"/*/; do
    [ -d "$d" ] || continue
    sample=$(basename "$d")
    chr="$d/chromosomal.fasta"
    [ -s "$chr" ] || { echo "  skip $sample"; continue; }

    bin_ids="$OUT_BASE/binned/${sample}.ids"
    : > "$bin_ids"
    das_dir="$DAS_BASE/$sample/_DASTool_bins"
    if [ -d "$das_dir" ]; then
        for f in "$das_dir"/*.fa; do
            [ -f "$f" ] || continue
            grep '^>' "$f" | sed -e 's/^>//' -e 's/ .*//' >> "$bin_ids"
        done
    fi
    sort -u "$bin_ids" -o "$bin_ids"

    # Append binned & unbinned to global FAs, prefixing IDs with sample
    awk -v s="$sample" -v ids="$bin_ids" -v bo="$BINNED_FA" -v uo="$UNBINNED_FA" '
        BEGIN{ while((getline l < ids)>0) keep[l]=1 }
        /^>/{
            id=$1; sub(/^>/,"",id)
            out = (id in keep) ? bo : uo
            print ">" s "|" id > out
            next
        }
        { print > out }
    ' "$chr"

    n_b=$(grep -c '^>' "$BINNED_FA")
    n_u=$(grep -c '^>' "$UNBINNED_FA")
    echo "    $sample done (running totals: binned=$n_b, unbinned=$n_u)"
done

echo "[$(date '+%F %T')] DONE"
echo "  binned   : $(grep -c '^>' $BINNED_FA) contigs"
echo "  unbinned : $(grep -c '^>' $UNBINNED_FA) contigs"
