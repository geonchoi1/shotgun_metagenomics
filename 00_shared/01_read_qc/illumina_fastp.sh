#!/bin/bash
# === Illumina read trimming with fastp ===
# Input:  $PROJECT_DIR/00_input/reads/<SAMPLE>_R1.fastq.gz + <SAMPLE>_R2.fastq.gz
# Output: $PROJECT_DIR/01_qc/01_fastp/<SAMPLE>_R{1,2}.fastp.fastq.gz
#         $PROJECT_DIR/01_qc/01_fastp/<SAMPLE>.{json,html}  (QC reports)

set -e
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source $REPO/config/db_paths.sh
source $REPO/config/envs.sh
source $REPO/config/threads.sh

: ${PROJECT_DIR:?ERROR: export PROJECT_DIR=/path/to/your/project}

IN_DIR=$PROJECT_DIR/00_input/reads
OUT_DIR=$PROJECT_DIR/01_qc/01_fastp
mkdir -p $OUT_DIR

activate_env $ENV_FASTP

echo "[$(date '+%F %T')] fastp — Illumina trimming"
# Auto-discover samples from *_R1.fastq.gz pattern
for r1 in $IN_DIR/*_R1.fastq.gz; do
    [ -f "$r1" ] || continue
    sample=$(basename $r1 _R1.fastq.gz)
    r2=$IN_DIR/${sample}_R2.fastq.gz
    [ -f "$r2" ] || { echo "  skip $sample: $r2 missing"; continue; }

    out1=$OUT_DIR/${sample}_R1.fastp.fastq.gz
    out2=$OUT_DIR/${sample}_R2.fastp.fastq.gz
    [ -f $out1 ] && [ -f $out2 ] && { echo "  $sample already done"; continue; }

    echo "  $sample"
    fastp \
      -i $r1 -I $r2 \
      -o $out1 -O $out2 \
      --json $OUT_DIR/${sample}.json \
      --html $OUT_DIR/${sample}.html \
      -w $THREADS \
      --detect_adapter_for_pe \
      2> $OUT_DIR/${sample}.log
done

echo "[$(date '+%F %T')] DONE"
