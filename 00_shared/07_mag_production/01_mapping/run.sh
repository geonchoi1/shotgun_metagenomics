#!/bin/bash
# === 07.01 Mapping clean reads -> chromosomal.fasta (per sample, self-mapping) ===
# HiFi      -> minimap2 -ax map-hifi
# Illumina  -> bwa-mem (paired-end)
#
# Input:  $PROJECT/00_shared/06_chromosomal_extract/<SAMPLE>/chromosomal.fasta
#         $PROJECT/00_shared/01_read_qc/dehuman/<SAMPLE>_clean.fastq.gz             (HiFi)
#         $PROJECT/00_shared/01_read_qc/dehuman/<SAMPLE>_clean_R{1,2}.fastq.gz      (Illumina)
# Output: $PROJECT/00_shared/07_mag_production/01_mapping/<SAMPLE>.bam (sorted, indexed)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}
: ${READ_TYPE:?ERROR: export READ_TYPE=hifi or illumina}

REF_BASE=$PROJECT/00_shared/06_chromosomal_extract
READ_BASE=$PROJECT/00_shared/01_read_qc/dehuman
OUT_DIR=$PROJECT/00_shared/07_mag_production/01_mapping
mkdir -p "$OUT_DIR"

if [ "$READ_TYPE" = "hifi" ]; then
    activate_env "$ENV_MINIMAP2"
    echo "[$(date '+%F %T')] minimap2 -ax map-hifi (self mapping)"
    for d in "$REF_BASE"/*/; do
        [ -d "$d" ] || continue
        sample=$(basename "$d")
        ref="$d/chromosomal.fasta"
        fq="$READ_BASE/${sample}_clean.fastq.gz"
        [ -s "$ref" ] || { echo "  skip $sample: no ref"; continue; }
        [ -s "$fq"  ] || { echo "  skip $sample: no reads"; continue; }
        bam="$OUT_DIR/${sample}.bam"
        if [ -s "$bam" ] && [ -s "$bam.bai" ]; then
            echo "  $sample already done"; continue
        fi
        echo "  $sample"
        minimap2 -t "$THREADS_MINIMAP2" -ax map-hifi "$ref" "$fq" \
            | samtools sort -@ "$THREADS" -o "$bam" -
        samtools index -@ "$THREADS" "$bam"
    done
elif [ "$READ_TYPE" = "illumina" ]; then
    activate_env "$ENV_BOWTIE2"   # bwa + samtools commonly co-installed here; adjust if needed
    echo "[$(date '+%F %T')] bwa-mem (Illumina self mapping)"
    for d in "$REF_BASE"/*/; do
        [ -d "$d" ] || continue
        sample=$(basename "$d")
        ref="$d/chromosomal.fasta"
        r1="$READ_BASE/${sample}_clean_R1.fastq.gz"
        r2="$READ_BASE/${sample}_clean_R2.fastq.gz"
        [ -s "$ref" ] || { echo "  skip $sample: no ref"; continue; }
        [ -s "$r1"  ] || { echo "  skip $sample: no R1"; continue; }
        [ -s "$r2"  ] || { echo "  skip $sample: no R2"; continue; }
        bam="$OUT_DIR/${sample}.bam"
        if [ -s "$bam" ] && [ -s "$bam.bai" ]; then
            echo "  $sample already done"; continue
        fi
        echo "  $sample"
        [ -f "$ref.bwt" ] || bwa index "$ref"
        bwa mem -t "$THREADS" "$ref" "$r1" "$r2" \
            | samtools sort -@ "$THREADS" -o "$bam" -
        samtools index -@ "$THREADS" "$bam"
    done
else
    echo "ERROR: READ_TYPE must be hifi or illumina" >&2; exit 2
fi

echo "[$(date '+%F %T')] DONE"
