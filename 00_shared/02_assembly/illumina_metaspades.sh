#!/bin/bash
# === metaSPAdes assembly (Illumina paired-end) ===
# Input:  $PROJECT/00_shared/01_read_qc/dehuman/<SAMPLE>_clean_R{1,2}.fastq.gz
# Output: $PROJECT/00_shared/02_assembly/metaflye/<SAMPLE>/scaffolds.fasta (canonical assembly)
#         $PROJECT/00_shared/02_assembly/metaflye/<SAMPLE>/contigs.fasta
#         (per-sample assembly subfolder)

set -e
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"




: ${PROJECT:?ERROR: export PROJECT=/path/to/your/project}

IN_DIR=$PROJECT/00_shared/01_read_qc/dehuman
OUT_BASE=$PROJECT/00_shared/02_assembly/metaflye
mkdir -p $OUT_BASE

activate_env $ENV_METASPADES

MEM_GB=${MEM_GB:-200}

echo "[$(date '+%F %T')] metaSPAdes assembly"
for r1 in $IN_DIR/*_clean_R1.fastq.gz; do
    [ -f "$r1" ] || continue
    sample=$(basename $r1 _clean_R1.fastq.gz)
    r2=$IN_DIR/${sample}_clean_R2.fastq.gz

    out_dir=$OUT_BASE/$sample
    [ -f $out_dir/scaffolds.fasta ] && { echo "  $sample already done"; continue; }

    echo "  $sample"
    metaspades.py \
      -1 $r1 -2 $r2 \
      -o $out_dir \
      -t $THREADS \
      -m $MEM_GB \
      2>&1 | tail -50 > $out_dir/run.log 2>&1 || true

    # Canonical output handle: assembly.fasta -> scaffolds.fasta
    [ -f $out_dir/scaffolds.fasta ] && ln -sf scaffolds.fasta $out_dir/assembly.fasta
done

echo "[$(date '+%F %T')] DONE"
echo "  assembly counts:"
for d in $OUT_BASE/*/; do
    [ -f $d/scaffolds.fasta ] && echo "    $(basename $d): $(grep -c '^>' $d/scaffolds.fasta) contigs"
done
