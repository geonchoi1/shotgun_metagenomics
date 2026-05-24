#!/bin/bash
# === ALT for rep ===
# PlasmidFinder CLI: clinical-biased Enterobacterales replicon DB.
# Very few hits expected on environmental plasmid (we observed 7 at default, 9 at loose).
# Runs both default (0.9 identity / 0.6 coverage) and loose (0.8 / 0.6) — paper-reported pair.

set -e
source ~/anaconda3/etc/profile.d/conda.sh
conda activate plasmidfinder

INPUT_FNA=${INPUT_FNA:-../inputs/dereplicated.fna}
PF_DB=${PF_DB:-/mnt/nas/DB/geon/plasmidfinder_db2}
OUT_DIR=${OUT_DIR:-../outputs/01_rep_plasmidfinder}
mkdir -p $OUT_DIR/{full,loose80}

echo "[$(date '+%F %T')] PlasmidFinder default (-t 0.9 -l 0.6)"
plasmidfinder.py \
  -i $INPUT_FNA \
  -o $OUT_DIR/full \
  -p $PF_DB \
  -mp $(which blastn) \
  -t 0.90 -l 0.60 -x -q

echo "[$(date '+%F %T')] PlasmidFinder loose (-t 0.8 -l 0.6)"
plasmidfinder.py \
  -i $INPUT_FNA \
  -o $OUT_DIR/loose80 \
  -p $PF_DB \
  -mp $(which blastn) \
  -t 0.80 -l 0.60 -x -q

echo "DONE"
echo "  default hits: $(($(wc -l < $OUT_DIR/full/results_tab.tsv) - 1))"
echo "  loose hits: $(($(wc -l < $OUT_DIR/loose80/results_tab.tsv) - 1))"
