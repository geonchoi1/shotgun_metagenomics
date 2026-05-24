#!/bin/bash
# === BiG-SCAPE on antiSMASH BGC GBK outputs ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

AS=$PROJECT/plasmid/17_antismash
OUT=$PROJECT/plasmid/18_bigscape
mkdir -p $OUT $OUT/input_gbks

[ -d $OUT/network_files ] && { echo "skip (exists)"; exit 0; }

# Collect region GBKs from both circ and frag
find $AS -name "*.region*.gbk" -exec cp -nf {} $OUT/input_gbks/ \;
N_GBK=$(ls $OUT/input_gbks/*.gbk 2>/dev/null | wc -l)
echo "[$(date '+%F %T')] BiG-SCAPE input region GBKs: $N_GBK"
[ $N_GBK -eq 0 ] && { echo "no BGCs — skip"; exit 0; }

activate_env "$ENV_BIGSCAPE"

bigscape \
  -i $OUT/input_gbks \
  -o $OUT/network_files \
  -c $THREADS \
  --mibig \
  --mix --no_classify || \
bigscape.py \
  -i $OUT/input_gbks \
  -o $OUT/network_files \
  -c $THREADS \
  --mibig --mix --no_classify

echo "[$(date '+%F %T')] BiG-SCAPE done"
