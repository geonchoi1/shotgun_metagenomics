#!/bin/bash
# === 15 DefenseFinder (anti-phage defense systems) on UB ORFs ===
# Input:  $PROJECT/unbinned/03_master_orf/all/master.faa
# Output: $PROJECT/unbinned/15_defensefinder/all/defense_finder_systems.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/unbinned/03_master_orf/all/master.faa
OUT=$PROJECT/unbinned/15_defensefinder/all
mkdir -p "$OUT"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }
[ -d "$DEFENSEFINDER_MODELS" ] || { echo "ERROR: missing $DEFENSEFINDER_MODELS" >&2; exit 1; }

if [ -s "$OUT/defense_finder_systems.tsv" ]; then
    echo "[$(date '+%F %T')] UB DefenseFinder already done — skip"; exit 0
fi

activate_env "$ENV_DEFENSEFINDER"

echo "[$(date '+%F %T')] defense-finder run on UB ORFs"
defense-finder run \
    --models-dir "$DEFENSEFINDER_MODELS" \
    --out-dir "$OUT" \
    --workers "$THREADS" \
    --preserve-raw \
    "$IN" > "$OUT/defensefinder.log" 2>&1

# defense-finder produces <basename>_defense_finder_*.tsv — symlink for consistency
base=$(basename "$IN" .faa)
for f in defense_finder_systems defense_finder_genes defense_finder_hmmer; do
    src=$(find "$OUT" -maxdepth 2 -name "${base}_${f}.tsv" | head -1)
    [ -n "$src" ] && [ -f "$src" ] && ln -sf "$src" "$OUT/${f}.tsv"
done

n=$(awk 'NR>1' "$OUT/defense_finder_systems.tsv" 2>/dev/null | wc -l || echo 0)
echo "[$(date '+%F %T')] DONE — defense systems: $n"
