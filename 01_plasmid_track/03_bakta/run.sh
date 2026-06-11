#!/bin/bash
# === Bakta annotation on dereplicated plasmids ===
# - circ: --complete (Bakta treats contigs as full circular plasmids)
# - frag: default (linear, may be incomplete)
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

CIRC=$PROJECT/01_plasmid_track/02_drep/circ.fna
FRAG=$PROJECT/01_plasmid_track/02_drep/frag.fna
OUT=$PROJECT/01_plasmid_track/03_bakta
mkdir -p $OUT/circ $OUT/frag

activate_env "$ENV_BAKTA"

# ---- Circular: --complete ----
if [ ! -s $OUT/circ/plasmid_circ.faa ] && [ -s $CIRC ] && [ "$(grep -c '^>' $CIRC)" -gt 0 ]; then
  echo "[$(date '+%F %T')] Bakta circ (--complete): $(grep -c '^>' $CIRC) contigs"
  bakta --complete --force --keep-contig-headers \
        --db "$BAKTA_DB" \
        --output $OUT/circ --prefix plasmid_circ \
        --threads $THREADS_BAKTA \
        $CIRC > $OUT/circ/bakta.log 2>&1
else
  echo "[$(date '+%F %T')] circ output exists or empty — skip"
fi

# ---- Fragmented: default ----
if [ ! -s $OUT/frag/plasmid_frag.faa ] && [ -s $FRAG ] && [ "$(grep -c '^>' $FRAG)" -gt 0 ]; then
  echo "[$(date '+%F %T')] Bakta frag (default): $(grep -c '^>' $FRAG) contigs"
  bakta --force --keep-contig-headers \
        --db "$BAKTA_DB" \
        --output $OUT/frag --prefix plasmid_frag \
        --threads $THREADS_BAKTA \
        $FRAG > $OUT/frag/bakta.log 2>&1
else
  echo "[$(date '+%F %T')] frag output exists or empty — skip"
fi

echo "[$(date '+%F %T')] DONE — circ ORFs=$(grep -c '^>' $OUT/circ/plasmid_circ.faa 2>/dev/null || echo 0), frag ORFs=$(grep -c '^>' $OUT/frag/plasmid_frag.faa 2>/dev/null || echo 0)"
