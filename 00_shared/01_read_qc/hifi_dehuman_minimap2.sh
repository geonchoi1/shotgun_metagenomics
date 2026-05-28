#!/bin/bash
# === HiFi (PacBio long read) human removal (minimap2 vs GRCh38) ===
# Input:  $PROJECT/00_input/reads/<SAMPLE>.fastq.gz   (or .hifi.fastq.gz)
# Output: $PROJECT/00_shared/01_read_qc/dehuman/<SAMPLE>_clean.fastq.gz
#
# Note: HiFi reads are already adapter-trimmed by the PacBio CCS pipeline (lima).
# fastp is NOT used for HiFi — dehuman is the only QC step here.

set -e
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"




: ${PROJECT:?ERROR: export PROJECT=/path/to/your/project}

IN_DIR=$PROJECT/00_input/reads
OUT_DIR=$PROJECT/00_shared/01_read_qc/dehuman
mkdir -p $OUT_DIR

activate_env $ENV_MINIMAP2

echo "[$(date '+%F %T')] minimap2 dehuman (HiFi vs GRCh38)"
# Auto-discover samples: anything ending .fastq.gz that doesn't match _R[12]
for fq in $IN_DIR/*.fastq.gz; do
    [ -f "$fq" ] || continue
    bn=$(basename $fq .fastq.gz)
    [[ $bn == *_R1 || $bn == *_R2 ]] && continue        # skip Illumina paired
    sample=${bn%.hifi}                                   # tolerate .hifi.fastq.gz

    out=$OUT_DIR/${sample}_clean.fastq.gz
    [ -f $out ] && { echo "  $sample already done"; continue; }

    echo "  $sample"
    minimap2 -ax map-hifi -t $THREADS $GRCH38 $fq 2> $OUT_DIR/${sample}.minimap2.log \
      | samtools view -b -f 4 -@ 8 \
      | samtools fastq -0 $out -
done

echo "[$(date '+%F %T')] DONE"
