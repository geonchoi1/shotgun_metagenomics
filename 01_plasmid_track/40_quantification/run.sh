#!/bin/bash
# === minimap2 + CoverM TPM quantification per sample ===
# Reads clean (dehuman) reads from $PROJECT/00_shared/01_read_qc/dehuman/
#   HiFi:     <SAMPLE>_clean.fastq.gz
#   Illumina: <SAMPLE>_clean_R1.fastq.gz + <SAMPLE>_clean_R2.fastq.gz
# Maps each sample's reads to dereplicated.fna
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

REF=$PROJECT/01_plasmid_track/02_drep/dereplicated.fna
READS_DIR=${READS_DIR:-$PROJECT/00_shared/01_read_qc/dehuman}
PRESET=${MINIMAP2_PRESET:-map-hifi}  # use sr for short reads
OUT=$PROJECT/01_plasmid_track/40_quantification
mkdir -p $OUT $OUT/bam

[ -s $REF ] || { echo "ERROR: $REF missing"; exit 1; }

activate_env "$ENV_COVERM"

# Auto-detect samples from clean dehuman reads (flat: <SAMPLE>_clean.fastq.gz, or
# <SAMPLE>_clean_R1/_R2.fastq.gz for Illumina).
shopt -s nullglob
SAMPLES=$(ls "$READS_DIR"/*_clean.fastq.gz "$READS_DIR"/*_clean_R1.fastq.gz 2>/dev/null \
          | xargs -n1 basename 2>/dev/null | sed -E 's/_clean(_R1)?\.fastq\.gz$//' | sort -u)
for sample in $SAMPLES; do
  bam=$OUT/bam/${sample}.bam
  if [ -s $bam ]; then echo "skip $sample (bam exists)"; continue; fi
  R1=$READS_DIR/${sample}_clean_R1.fastq.gz
  R2=$READS_DIR/${sample}_clean_R2.fastq.gz
  HIFI=$READS_DIR/${sample}_clean.fastq.gz
  echo "[$(date '+%F %T')] $sample → minimap2"
  if [ -e "$R1" ] && [ -e "$R2" ]; then
    minimap2 -ax sr -t $THREADS_MINIMAP2 $REF "$R1" "$R2" \
      | samtools sort -@ $THREADS_MINIMAP2 -o $bam -
  elif [ -e "$HIFI" ]; then
    minimap2 -ax $PRESET -t $THREADS_MINIMAP2 $REF "$HIFI" \
      | samtools sort -@ $THREADS_MINIMAP2 -o $bam -
  else
    echo "no reads for $sample"; continue
  fi
  samtools index -@ $THREADS_MINIMAP2 $bam
done

# CoverM per sample
for bam in $OUT/bam/*.bam; do
  sample=$(basename $bam .bam)
  outfile=$OUT/coverm_${sample}.tsv
  [ -s $outfile ] && continue
  echo "[$(date '+%F %T')] CoverM $sample"
  coverm contig \
    --bam-files $bam \
    --methods mean tpm rpkm \
    --min-read-percent-identity 90 \
    --min-covered-fraction 50 \
    --threads $THREADS \
    --output-file $outfile
done

# Merge into one wide TPM matrix
python3 - <<PYEOF
import os, glob, pandas as pd
OUT="$OUT"
dfs=[]
for tf in sorted(glob.glob(f"{OUT}/coverm_*.tsv")):
    sample=os.path.basename(tf).replace('coverm_','').replace('.tsv','')
    df=pd.read_csv(tf, sep='\t')
    tpm_col=[c for c in df.columns if 'TPM' in c]
    if not tpm_col: continue
    s=df.set_index(df.columns[0])[tpm_col[0]].rename(sample)
    dfs.append(s)
if dfs:
    pd.concat(dfs, axis=1).fillna(0).to_csv(f"{OUT}/tpm_matrix.tsv", sep='\t')
    print(f"merged: {len(dfs)} samples")
PYEOF

echo "[$(date '+%F %T')] DONE"
