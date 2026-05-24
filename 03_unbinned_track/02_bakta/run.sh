#!/bin/bash
# === 02 Bakta on unbinned chromosomal contigs (default; no --complete) ===
# Input:  $PROJECT/unbinned/01_raw_fasta/all/unbinned.fna
# Output: $PROJECT/unbinned/02_bakta/unbinned/unbinned.{faa,ffn,fna,gff3,tsv}

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/unbinned/01_raw_fasta/all/unbinned.fna
OUT=$PROJECT/unbinned/02_bakta/unbinned
mkdir -p "$OUT"

[ -s "$IN" ] || { echo "ERROR: missing $IN — run 01_raw_fasta first" >&2; exit 1; }

if [ -s "$OUT/unbinned.gff3" ]; then
    echo "[$(date '+%F %T')] unbinned Bakta already done — skip"
    exit 0
fi

activate_env "$ENV_BAKTA"

echo "[$(date '+%F %T')] Bakta on unbinned (default, no --complete)"
bakta --db "$BAKTA_DB" \
      --output "$OUT" \
      --prefix unbinned \
      --threads "$THREADS_BAKTA" \
      --locus-tag UB \
      --force \
      "$IN" > "$OUT/bakta.log" 2>&1

echo "[$(date '+%F %T')] DONE"
echo "  faa : $OUT/unbinned.faa"
echo "  gff : $OUT/unbinned.gff3"
echo "  ffn : $OUT/unbinned.ffn"
