#!/bin/bash
# === BacMet (metal/biocide resistance) via DIAMOND blastp ===
# Filter: id≥70, qcov≥50, scov≥50, evalue≤1e-5
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/plasmid/04_master_orf/all/master.faa
OUT=$PROJECT/plasmid/10_bacmet
mkdir -p $OUT

[ -s $OUT/bacmet.tsv ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_DIAMOND"

diamond blastp \
  --query $FAA --db $BACMET_DMND \
  --id 70 --query-cover 50 --subject-cover 50 \
  --evalue 1e-5 --max-target-seqs 1 \
  --threads $THREADS --memory-limit ${MEM_DIAMOND%G} \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
  --out $OUT/bacmet.tsv

echo "[$(date '+%F %T')] BacMet hits: $(wc -l < $OUT/bacmet.tsv)"
