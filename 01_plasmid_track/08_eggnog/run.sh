#!/bin/bash
# === eggNOG-mapper (functional + ortholog assignment) ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

FAA=$PROJECT/01_plasmid_track/04_master_orf/all/master.faa
OUT=$PROJECT/01_plasmid_track/08_eggnog
mkdir -p $OUT

[ -s $OUT/eggnog.emapper.annotations ] && { echo "skip (exists)"; exit 0; }

activate_env "$ENV_EGGNOG"

emapper.py \
  -i $FAA --itype proteins \
  --data_dir $EGGNOG_DATA \
  --cpu $THREADS \
  -o eggnog --output_dir $OUT \
  --override

echo "[$(date '+%F %T')] eggNOG annotations: $(grep -vc '^#' $OUT/eggnog.emapper.annotations)"
