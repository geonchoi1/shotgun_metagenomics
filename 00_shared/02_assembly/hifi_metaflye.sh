#!/bin/bash
# === metaFlye assembly (PacBio HiFi) + ≥1kb size filter ===
# Input:  $PROJECT/00_shared/01_read_qc/dehuman/<SAMPLE>_clean.fastq.gz
# Output: $PROJECT/00_shared/02_assembly/metaflye/<SAMPLE>/assembly.fasta       (raw, all contigs)
#         $PROJECT/00_shared/02_assembly/metaflye/<SAMPLE>/assembly_1kb.fasta   (≥1kb filter; used by ALL downstream)
#         $PROJECT/00_shared/02_assembly/metaflye/<SAMPLE>/assembly_info.txt    (circ flag in col 4)
#         $PROJECT/00_shared/02_assembly/metaflye/<SAMPLE>/assembly_graph.gfa
#
# All downstream steps (geNomad plasmid/virus identification, MAG binning,
# unbinned chromosomal extraction) operate on assembly_1kb.fasta.
# Contigs <1kb carry too little information for reliable annotation / binning
# and are removed here with seqkit.

set -e
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"




: ${PROJECT:?ERROR: export PROJECT=/path/to/your/project}

IN_DIR=$PROJECT/00_shared/01_read_qc/dehuman
OUT_BASE=$PROJECT/00_shared/02_assembly/metaflye
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

echo "[$(date '+%F %T')] metaFlye DONE. Applying ≥1kb size filter (seqkit)..."

# === ≥1kb size filter ===
activate_env "$ENV_SEQKIT"
for d in $OUT_BASE/*/; do
    [ -f $d/assembly.fasta ] || continue
    [ -f $d/assembly_1kb.fasta ] && continue
    seqkit seq --min-len 1000 $d/assembly.fasta -o $d/assembly_1kb.fasta
done

echo "  assembly + ≥1kb + circular counts:"
for d in $OUT_BASE/*/; do
    [ -f $d/assembly.fasta ] || continue
    s=$(basename $d)
    n=$(grep -c '^>' $d/assembly.fasta)
    n1=$(grep -c '^>' $d/assembly_1kb.fasta 2>/dev/null || echo 0)
    c=$(tail -n +2 $d/assembly_info.txt 2>/dev/null | awk '$4=="Y"' | wc -l)
    echo "    $s: $n contigs ($n1 ≥1kb), $c circular"
done
