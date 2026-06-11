#!/bin/bash
# === 21 ICEberg3 BLASTn (mobile element / ICE / IME) ===
# Input:  $PROJECT/02_mag_track/19_isescan/work/all_mag.fna  (per-MAG concat)
# Output: $PROJECT/02_mag_track/21_iceberg3/iceberg3.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

RAW_BASE=$PROJECT/02_mag_track/01_raw_fasta
OUT=$PROJECT/02_mag_track/21_iceberg3
mkdir -p "$OUT"

CONCAT=$PROJECT/02_mag_track/19_isescan/work/all_mag.fna
if [ ! -s "$CONCAT" ]; then
    echo "[$(date '+%F %T')] building per-MAG concat FNA"
    mkdir -p "$(dirname "$CONCAT")"
    : > "$CONCAT"
    for topo in circ frag; do
        for fa in "$RAW_BASE/$topo"/*.fa; do
            [ -f "$fa" ] || continue
            mag=$(basename "$fa" .fa)
            awk -v m="$mag" '/^>/{sub(/^>/,">"m"|"); print; next}{print}' "$fa" >> "$CONCAT"
        done
    done
fi

activate_env "$ENV_DIAMOND"  # blastn lives in base by convention

if [ ! -f "${ICEBERG3_DB}.nhr" ] && [ ! -f "${ICEBERG3_DB}.nin" ]; then
    echo "  building BLAST nucleotide DB"
    makeblastdb -in "$ICEBERG3_DB" -dbtype nucl -out "$ICEBERG3_DB"
fi

echo "[$(date '+%F %T')] blastn vs ICEberg3 (e<=1e-5, qcov>=50, pident>=70)"
blastn \
    -query "$CONCAT" \
    -db "$ICEBERG3_DB" \
    -out "$OUT/iceberg3.tsv" \
    -evalue 1e-5 \
    -perc_identity 70 \
    -qcov_hsp_perc 50 \
    -num_threads "$THREADS_BLAST" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp" \
    -max_target_seqs 5

n=$(wc -l < "$OUT/iceberg3.tsv")
echo "[$(date '+%F %T')] DONE — $n hits"
