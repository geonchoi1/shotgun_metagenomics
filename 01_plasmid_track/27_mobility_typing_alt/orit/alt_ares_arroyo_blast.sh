#!/bin/bash
# === ALT for oriT ===
# Ares-Arroyo 2023 (10.1093/nar/gkad084) 91-oriT BLAST.
# Lightweight; complementary to oriTfinder2. Strict filter id≥80% + cov≥80% (paper threshold).

set -e
source ~/anaconda3/etc/profile.d/conda.sh
conda activate base

THREADS=${THREADS:-16}
INPUT_FNA=${INPUT_FNA:-../inputs/dereplicated.fna}
ARES_DB=${ARES_DB:-/home/gchoi/tools/ares_arroyo_oriT/ares_arroyo_91_oriT.fna}
OUT_DIR=${OUT_DIR:-../outputs/03_orit_ares_arroyo}
mkdir -p $OUT_DIR

echo "[$(date '+%F %T')] BLAST plasmid vs Ares-Arroyo 91 oriT"
blastn -task blastn-short -evalue 0.01 \
  -query $INPUT_FNA \
  -db $ARES_DB \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
  -num_threads $THREADS \
  -out $OUT_DIR/orit_raw.tsv

# Filter: id > 80%, query-coverage > 80%
awk -F'\t' '$3 > 80 && ($4/$13*100) > 80' $OUT_DIR/orit_raw.tsv > $OUT_DIR/orit_filtered.tsv
cut -f1 $OUT_DIR/orit_filtered.tsv | sort -u > $OUT_DIR/contigs_with_orit.txt
echo "  raw hits: $(wc -l < $OUT_DIR/orit_raw.tsv)"
echo "  filtered hits: $(wc -l < $OUT_DIR/orit_filtered.tsv)"
echo "  unique plasmid with oriT (Ares-Arroyo): $(wc -l < $OUT_DIR/contigs_with_orit.txt)"
