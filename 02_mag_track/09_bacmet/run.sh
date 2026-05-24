#!/bin/bash
# === 09 BacMet (metal/biocide resistance) via DIAMOND blastp ===
# Filters: --id 70  --query-cover 50  --subject-cover 50  -e 1e-5
# Input:  $PROJECT/mag/03_master_orf/all/master.faa
# Output: $PROJECT/mag/09_bacmet/bacmet.tsv (filtered)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/mag/03_master_orf/all/master.faa
OUT=$PROJECT/mag/09_bacmet
mkdir -p "$OUT"

[ -s "$FAA" ] || { echo "ERROR: $FAA empty" >&2; exit 1; }

activate_env "$ENV_DIAMOND"

echo "[$(date '+%F %T')] DIAMOND blastp vs BacMet (id>=70, qcov>=50, scov>=50, e<=1e-5)"
diamond blastp \
    -q "$FAA" -d "$BACMET_DMND" \
    -o "$OUT/bacmet.tsv" \
    -p "$THREADS_BLAST" \
    -e 1e-5 --id 70 --query-cover 50 --subject-cover 50 \
    --max-target-seqs 5 \
    -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
    --memory-limit 32

n=$(wc -l < "$OUT/bacmet.tsv")
echo "[$(date '+%F %T')] DONE — $n hits"
