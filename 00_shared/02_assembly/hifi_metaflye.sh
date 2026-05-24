#!/bin/bash
# === metaFlye assembly (PacBio HiFi) ===
# Input:  $PROJECT_DIR/01_qc/02_dehuman/<SAMPLE>_clean.fastq.gz
# Output: $PROJECT_DIR/02_assembly/<SAMPLE>/assembly.fasta
#         $PROJECT_DIR/02_assembly/<SAMPLE>/assembly_info.txt  (circ flag in col 4)
#         $PROJECT_DIR/02_assembly/<SAMPLE>/assembly_graph.gfa

set -e
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source $REPO/config/db_paths.sh
source $REPO/config/envs.sh
source $REPO/config/threads.sh

: ${PROJECT_DIR:?ERROR: export PROJECT_DIR=/path/to/your/project}

IN_DIR=$PROJECT_DIR/01_qc/02_dehuman
OUT_BASE=$PROJECT_DIR/02_assembly
mkdir -p $OUT_BASE

activate_env $ENV_METAFLYE

echo "[$(date '+%F %T')] metaFlye --meta --pacbio-hifi"
for fq in $IN_DIR/*_clean.fastq.gz; do
    [ -f "$fq" ] || continue
    sample=$(basename $fq _clean.fastq.gz)

    out_dir=$OUT_BASE/$sample
    [ -f $out_dir/assembly.fasta ] && { echo "  $sample already done"; continue; }

    echo "  $sample"
    flye \
      --meta --pacbio-hifi \
      $fq \
      --out-dir $out_dir \
      --threads $THREADS \
      2>&1 | tail -50 > $out_dir/run.log 2>&1 || true
done

echo "[$(date '+%F %T')] DONE"
echo "  assembly + circular counts:"
for d in $OUT_BASE/*/; do
    [ -f $d/assembly.fasta ] || continue
    s=$(basename $d)
    n=$(grep -c '^>' $d/assembly.fasta)
    c=$(tail -n +2 $d/assembly_info.txt 2>/dev/null | awk '$4=="Y"' | wc -l)
    echo "    $s: $n contigs, $c circular"
done
