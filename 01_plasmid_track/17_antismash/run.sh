#!/bin/bash
# === antiSMASH on plasmid FNA (circ + frag separately) ===
# Uses --fullhmmer --cc-mibig --genefinding-tool prodigal-m
set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

CIRC=$PROJECT/01_plasmid_track/02_drep/circ.fna
FRAG=$PROJECT/01_plasmid_track/02_drep/frag.fna
OUT=$PROJECT/01_plasmid_track/17_antismash
mkdir -p $OUT/circ $OUT/frag

activate_env "$ENV_ANTISMASH"

run_one() {
  local fna=$1; local out=$2
  [ -s $fna ] || { echo "skip — empty input: $fna"; return 0; }
  [ -s $out/index.html ] && { echo "skip $out (exists)"; return 0; }
  antismash \
    --fullhmmer --cc-mibig \
    --genefinding-tool prodigal-m \
    --output-dir $out \
    --cpus $THREADS \
    $fna
}

run_one $CIRC $OUT/circ
run_one $FRAG $OUT/frag

echo "[$(date '+%F %T')] antiSMASH circ regions: $(ls $OUT/circ/*.region*.gbk 2>/dev/null | wc -l), frag: $(ls $OUT/frag/*.region*.gbk 2>/dev/null | wc -l)"
