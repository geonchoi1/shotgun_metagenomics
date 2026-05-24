#!/bin/bash
# === IntegronFinder2: --circ on circular, --linear on fragmented (default local-max) ===
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

CIRC=$PROJECT/plasmid/02_drep/circ.fna
FRAG=$PROJECT/plasmid/02_drep/frag.fna
OUT=$PROJECT/plasmid/21_integronfinder
mkdir -p $OUT/circ $OUT/frag

activate_env "$ENV_INTEGRONFINDER"

if [ -s $CIRC ] && [ ! -s $OUT/circ/Results_Integron_Finder_circ/circ.integrons ]; then
  echo "[$(date '+%F %T')] IntegronFinder circ (--circ)"
  integron_finder --local-max --circ --cpu $THREADS \
    --outdir $OUT/circ $CIRC
fi
if [ -s $FRAG ] && [ ! -s $OUT/frag/Results_Integron_Finder_frag/frag.integrons ]; then
  echo "[$(date '+%F %T')] IntegronFinder frag (linear)"
  integron_finder --local-max --linear --cpu $THREADS \
    --outdir $OUT/frag $FRAG
fi

echo "[$(date '+%F %T')] DONE"
