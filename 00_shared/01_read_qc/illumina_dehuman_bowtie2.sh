#!/bin/bash
# === Illumina human read removal (bowtie2 vs GRCh38) ===
# Input:  $PROJECT/00_shared/01_read_qc/fastp/<SAMPLE>_R{1,2}.fastp.fastq.gz
# Output: $PROJECT/00_shared/01_read_qc/dehuman/<SAMPLE>_clean_R{1,2}.fastq.gz

set -e
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"




: ${PROJECT:?ERROR: export PROJECT=/path/to/your/project}

IN_DIR=$PROJECT/00_shared/01_read_qc/fastp
OUT_DIR=$PROJECT/00_shared/01_read_qc/dehuman
mkdir -p $OUT_DIR

activate_env $ENV_BOWTIE2

# Build bowtie2 index if missing
if [ ! -f ${GRCH38_BT2_IDX}.1.bt2 ] && [ ! -f ${GRCH38_BT2_IDX}.1.bt2l ]; then
    echo "[$(date '+%F %T')] Building bowtie2 index for GRCh38 (one-time, slow)"
    mkdir -p $(dirname $GRCH38_BT2_IDX)
    REF=$GRCH38
    [ "${REF##*.}" = "gz" ] && { gunzip -k $REF; REF=${REF%.gz}; }
    bowtie2-build --threads $THREADS $REF $GRCH38_BT2_IDX
fi

echo "[$(date '+%F %T')] bowtie2 dehuman"
for r1 in $IN_DIR/*_R1.fastp.fastq.gz; do
    [ -f "$r1" ] || continue
    sample=$(basename $r1 _R1.fastp.fastq.gz)
    r2=$IN_DIR/${sample}_R2.fastp.fastq.gz

    out1=$OUT_DIR/${sample}_clean_R1.fastq.gz
    out2=$OUT_DIR/${sample}_clean_R2.fastq.gz
    [ -f $out1 ] && [ -f $out2 ] && { echo "  $sample already done"; continue; }

    echo "  $sample"
    # --un-conc-gz writes pairs where neither mate mapped (= non-human)
    bowtie2 \
      -x $GRCH38_BT2_IDX \
      -1 $r1 -2 $r2 \
      -p $THREADS \
      --un-conc-gz $OUT_DIR/${sample}_clean_R%.fastq.gz \
      --no-unal \
      -S /dev/null \
      2> $OUT_DIR/${sample}.log
done

echo "[$(date '+%F %T')] DONE"
