#!/bin/bash
# === 11 TADB (toxin-antitoxin) DIAMOND blastp on UB ORFs ===
# Input:  $PROJECT/unbinned/03_master_orf/all/master.faa
# Output: $PROJECT/unbinned/11_tadb/all/tadb.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/unbinned/03_master_orf/all/master.faa
OUT=$PROJECT/unbinned/11_tadb/all
mkdir -p "$OUT"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }
[ -s "$TADB_DMND" ] || { echo "ERROR: missing $TADB_DMND" >&2; exit 1; }

if [ -s "$OUT/tadb.tsv" ]; then
    echo "[$(date '+%F %T')] UB TADB already done — skip"; exit 0
fi

activate_env "$ENV_DIAMOND"

echo "[$(date '+%F %T')] diamond blastp UB ORFs vs TADB"
diamond blastp \
    --query "$IN" \
    --db "$TADB_DMND" \
    --threads "$THREADS_BLAST" \
    --evalue 1e-10 \
    --id 40 \
    --query-cover 70 \
    --max-target-seqs 1 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
    -o "$OUT/tadb.tsv" 2> "$OUT/tadb.log"

n=$(wc -l < "$OUT/tadb.tsv")
echo "[$(date '+%F %T')] DONE — TADB hits: $n"
