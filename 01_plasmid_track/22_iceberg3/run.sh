#!/bin/bash
# === ICEberg3: BLAST plasmid vs ICE_combined ===
# Filter: pident≥80, alnlen≥500, qcov≥50
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FNA=$PROJECT/01_plasmid_track/02_drep/dereplicated.fna
OUT=$PROJECT/01_plasmid_track/22_iceberg3
mkdir -p $OUT

[ -s $OUT/iceberg3_filtered.tsv ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_DIAMOND"   # base — has BLAST+

# Build DB if .nhr missing
if [ ! -f ${ICEBERG3_DB}.nhr ]; then
  makeblastdb -in $ICEBERG3_DB -dbtype nucl -out $ICEBERG3_DB
fi

blastn -task megablast \
  -query $FNA -db $ICEBERG3_DB \
  -evalue 1e-5 -num_threads $THREADS_BLAST \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs' \
  -max_target_seqs 25 \
  -out $OUT/iceberg3.tsv

awk -F'\t' '$3>=80 && $4>=500 && $15>=50' $OUT/iceberg3.tsv > $OUT/iceberg3_filtered.tsv

echo "[$(date '+%F %T')] ICEberg3 hits: $(wc -l < $OUT/iceberg3.tsv), filtered: $(wc -l < $OUT/iceberg3_filtered.tsv)"
