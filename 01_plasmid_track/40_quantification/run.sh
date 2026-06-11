#!/bin/bash
# === minimap2 + CoverM TPM quantification per sample ===
# Reads sample reads from $PROJECT/reads/<sample>/{*.fq.gz | *_1.fq.gz,_2.fq.gz}
# Maps each sample's reads to dereplicated.fna
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

REF=$PROJECT/01_plasmid_track/02_drep/dereplicated.fna
READS_DIR=${READS_DIR:-$PROJECT/reads}
PRESET=${MINIMAP2_PRESET:-map-hifi}  # use sr for short reads
OUT=$PROJECT/01_plasmid_track/40_quantification
mkdir -p $OUT $OUT/bam

[ -s $REF ] || { echo "ERROR: $REF missing"; exit 1; }

activate_env "$ENV_COVERM"

# Auto-detect samples by checking $READS_DIR/<sample>/*.fq.gz
shopt -s nullglob
for SDIR in $READS_DIR/*/; do
  sample=$(basename $SDIR)
  bam=$OUT/bam/${sample}.bam
  if [ -s $bam ]; then echo "skip $sample (bam exists)"; continue; fi
  # Determine read layout
  R1=( ${SDIR}*_1.fq.gz ${SDIR}*_R1.fq.gz ${SDIR}*_1.fastq.gz )
  R2=( ${SDIR}*_2.fq.gz ${SDIR}*_R2.fq.gz ${SDIR}*_2.fastq.gz )
  HIFI=( ${SDIR}*.fq.gz ${SDIR}*.fastq.gz )
  echo "[$(date '+%F %T')] $sample → minimap2 -ax $PRESET"
  if [ -n "${R1[0]:-}" ] && [ -e "${R1[0]}" ] && [ -n "${R2[0]:-}" ] && [ -e "${R2[0]}" ]; then
    minimap2 -ax sr -t $THREADS_MINIMAP2 $REF ${R1[0]} ${R2[0]} \
      | samtools sort -@ $THREADS_MINIMAP2 -o $bam -
  else
    # single-file fallback (HiFi or merged)
    READ=""
    for f in "${HIFI[@]}"; do
      bn=$(basename $f)
      [[ $bn == *_1.* || $bn == *_R1.* || $bn == *_2.* || $bn == *_R2.* ]] && continue
      READ=$f; break
    done
    [ -z "$READ" ] && { echo "no reads for $sample"; continue; }
    minimap2 -ax $PRESET -t $THREADS_MINIMAP2 $REF $READ \
      | samtools sort -@ $THREADS_MINIMAP2 -o $bam -
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
