#!/bin/bash
# === DEFAULT for MPF (T4SS) ===
# MacSyFinder CONJScan/Plasmids gembase mode.
# Detects full T4SS_typeX (full MPF) + dCONJ_typeX (degraded) + MOB systems.
# Coluzzi 2022 framework — locked default.

set -e
source ~/anaconda3/etc/profile.d/conda.sh
conda activate macsyfinder

THREADS=${THREADS:-16}
INPUT_FAA=${INPUT_FAA:-../inputs/plasmidome.master.faa}
TOPOLOGY_FILE=${TOPOLOGY_FILE:-../inputs/topology.txt}
CONJSCAN_MODELS_DIR=${CONJSCAN_MODELS_DIR:-/mnt/nas/DB/geon/conjscan_models}
OUT_DIR=${OUT_DIR:-../outputs/04_mpf_conjscan}
mkdir -p $OUT_DIR

# Topology placeholder generation (assume linear; user must edit for known circular)
if [ ! -f $TOPOLOGY_FILE ]; then
    echo "WARN: $TOPOLOGY_FILE missing — generating all-linear placeholder"
    grep "^>" $INPUT_FAA | awk -F"_" '{print $1"_"$2}' | sort -u | \
      sed 's/^>//;s/[^A-Za-z0-9_.-]/_/g' | awk '{print $0": linear"}' > $TOPOLOGY_FILE
fi

echo "[$(date '+%F %T')] CONJScan gembase"
macsyfinder \
  --db-type gembase \
  --models-dir $CONJSCAN_MODELS_DIR \
  --models CONJScan/Plasmids all \
  --sequence-db $INPUT_FAA \
  --topology-file $TOPOLOGY_FILE \
  --out-dir $OUT_DIR/result \
  -w $THREADS

echo "DONE"
echo "  Model distribution:"
awk -F'\t' 'NR>1 && !/^#/{print $5}' $OUT_DIR/result/best_solution.tsv | sort | uniq -c | sort -rn | head -20
