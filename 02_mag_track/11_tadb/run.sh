#!/bin/bash
# === 11 TADB (toxin-antitoxin) via DIAMOND blastp ===
# Input:  $PROJECT/mag/03_master_orf/all/master.faa
# Output: $PROJECT/mag/11_tadb/tadb.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/mag/03_master_orf/all/master.faa
OUT=$PROJECT/mag/11_tadb
mkdir -p "$OUT"

[ -s "$FAA" ] || { echo "ERROR: $FAA empty" >&2; exit 1; }

activate_env "$ENV_DIAMOND"

echo "[$(date '+%F %T')] DIAMOND blastp vs TADB (e<=1e-5)"
diamond blastp \
    -q "$FAA" -d "$TADB_DMND" \
    -o "$OUT/tadb.tsv" \
    -p "$THREADS_BLAST" \
    -e 1e-5 \
    --max-target-seqs 5 \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp

n=$(wc -l < "$OUT/tadb.tsv")
echo "[$(date '+%F %T')] DONE — $n TA hits"
