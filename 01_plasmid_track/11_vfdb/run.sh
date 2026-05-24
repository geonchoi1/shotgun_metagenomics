#!/bin/bash
# === VFDB (virulence factors) via DIAMOND blastp ===
# Filter: id≥70, qcov≥50, scov≥50, evalue≤1e-5
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/plasmid/04_master_orf/all/master.faa
OUT=$PROJECT/plasmid/11_vfdb
mkdir -p $OUT

[ -s $OUT/vfdb.tsv ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_DIAMOND"

diamond blastp \
  --query $FAA --db $VFDB_DMND \
  --id 70 --query-cover 50 --subject-cover 50 \
  --evalue 1e-5 --max-target-seqs 1 \
  --threads $THREADS --memory-limit ${MEM_DIAMOND%G} \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
  --out $OUT/vfdb.tsv

echo "[$(date '+%F %T')] VFDB hits: $(wc -l < $OUT/vfdb.tsv)"
