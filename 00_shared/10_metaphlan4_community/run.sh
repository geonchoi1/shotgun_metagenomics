#!/bin/bash
# === MetaPhlAn4 community taxonomy (read-level, species + strain markers) ===
# Complement to Kraken2 (09): MetaPhlAn4 uses clade-specific marker genes
# Output: per-sample profile + merged table
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

READ_DIR=$PROJECT/shared/01_read_qc
OUT=$PROJECT/shared/10_metaphlan4_community
mkdir -p $OUT

[ -s $OUT/merged_metaphlan.tsv ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_METAPHLAN"

for sample in $SAMPLES; do
    fq=$READ_DIR/${sample}.clean.fastq.gz
    [ ! -s "$fq" ] && { echo "  WARN: no read for $sample"; continue; }
    
    mkdir -p $OUT/$sample
    [ -s $OUT/$sample/${sample}.profile.tsv ] && continue
    
    metaphlan $fq \
        --input_type fastq \
        --bowtie2db $METAPHLAN_DB \
        --bowtie2out $OUT/$sample/${sample}.mapout \
        --nproc $THREADS \
        --output_file $OUT/$sample/${sample}.profile.tsv
    
    echo "  $sample DONE — $(grep -v '^#' $OUT/$sample/${sample}.profile.tsv | wc -l) rows"
done

# Merge all samples
merge_metaphlan_tables.py $OUT/*/*.profile.tsv -o $OUT/merged_metaphlan.tsv
echo "[$(date '+%F %T')] MetaPhlAn4 merged: $(wc -l < $OUT/merged_metaphlan.tsv) rows"
