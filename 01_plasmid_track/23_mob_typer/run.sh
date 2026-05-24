#!/bin/bash
# === DEFAULT for rep ===
# mob_typer: replicon typing via PLSDB-based DIAMOND.
# Environmental plasmid coverage > PlasmidFinder.

set -e
source ~/anaconda3/etc/profile.d/conda.sh
conda activate mob_suite

THREADS=${THREADS:-16}
INPUT_FNA=${INPUT_FNA:-../inputs/dereplicated.fna}
MOB_DB=${MOB_DB:-/mnt/nas/DB/geon/mob_suite}
OUT_DIR=${OUT_DIR:-../outputs/01_rep_mob_typer}
mkdir -p $OUT_DIR

echo "[$(date '+%F %T')] mob_typer on $INPUT_FNA"
mob_typer \
  --infile $INPUT_FNA \
  --out_file $OUT_DIR/mobtyper_report.tsv \
  --num_threads $THREADS \
  --multi \
  -d $MOB_DB

echo "DONE — col 5=rep_type, 7=relaxase, 9=mpf, 11=orit, 13=mobility"
echo "  rep_type non-empty: $(awk -F'\t' 'NR>1 && $6!="-"' $OUT_DIR/mobtyper_report.tsv | wc -l)"
