#!/bin/bash
# === 21 ICEberg3 (ICE/IME) DIAMOND blastp on UB ORFs ===
# Input:  $PROJECT/03_unbinned_track/03_master_orf/all/master.faa
# Output: $PROJECT/03_unbinned_track/21_iceberg3/all/iceberg3.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/03_unbinned_track/03_master_orf/all/master.faa
OUT=$PROJECT/03_unbinned_track/21_iceberg3/all
mkdir -p "$OUT"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }

# ICEberg3 DB may be DIAMOND .dmnd or FASTA — handle both
DB="$ICEBERG3_DB"
if [ -s "${DB}.dmnd" ]; then DB="${DB}.dmnd"; fi
[ -s "$DB" ] || { echo "ERROR: missing $ICEBERG3_DB[.dmnd]" >&2; exit 1; }

if [ -s "$OUT/iceberg3.tsv" ]; then
    echo "[$(date '+%F %T')] UB ICEberg3 already done — skip"; exit 0
fi

activate_env "$ENV_DIAMOND"

# If FASTA, build dmnd on the fly into project-local dir (DB itself stays under DB_ROOT)
if [[ "$DB" != *.dmnd ]]; then
    DMND_LOCAL="$OUT/iceberg3_ref.dmnd"
    if [ ! -s "$DMND_LOCAL" ]; then
        diamond makedb --in "$DB" --db "${DMND_LOCAL%.dmnd}" --threads "$THREADS" \
            > "$OUT/makedb.log" 2>&1
    fi
    DB="$DMND_LOCAL"
fi

echo "[$(date '+%F %T')] diamond blastp UB ORFs vs ICEberg3"
diamond blastp \
    --query "$IN" \
    --db "$DB" \
    --threads "$THREADS_BLAST" \
    --evalue 1e-10 \
    --id 40 \
    --query-cover 70 \
    --max-target-seqs 1 \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp scovhsp \
    -o "$OUT/iceberg3.tsv" 2> "$OUT/iceberg3.log"

n=$(wc -l < "$OUT/iceberg3.tsv")
echo "[$(date '+%F %T')] DONE — ICEberg3 hits: $n"
