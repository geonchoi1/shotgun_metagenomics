#!/bin/bash
# === 07.08 dRep species-level cross-sample dereplication (-sa 0.95) ===
# Input:  $PROJECT/00_shared/07_mag_production/04_dastool/<SAMPLE>/_DASTool_bins/*.fa  (all samples pooled)
#         $PROJECT/00_shared/07_mag_production/05_checkm2/genomeInfo.csv
# Output: $PROJECT/00_shared/07_mag_production/08_drep_species/dereplicated_genomes/*.fa

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

DAS_BASE=$PROJECT/00_shared/07_mag_production/04_dastool
GINFO=$PROJECT/00_shared/07_mag_production/05_checkm2/genomeInfo.csv
OUT_DIR=$PROJECT/00_shared/07_mag_production/08_drep_species
INPUT_DIR=$OUT_DIR/input_bins
mkdir -p "$OUT_DIR" "$INPUT_DIR"

if [ -d "$OUT_DIR/dereplicated_genomes" ] && ls "$OUT_DIR/dereplicated_genomes"/*.fa 2>/dev/null | grep -q .; then
    echo "[$(date '+%F %T')] dRep 0.95 already done"; exit 0
fi

[ -s "$GINFO" ] || { echo "ERROR: $GINFO missing — run 05_checkm2 first" >&2; exit 2; }

echo "[$(date '+%F %T')] gather MAGs -> $INPUT_DIR"
for d in "$DAS_BASE"/*/_DASTool_bins/; do
    [ -d "$d" ] || continue
    sample=$(basename "$(dirname "$d")")
    for f in "$d"/*.fa; do
        [ -f "$f" ] || continue
        ln -sf "$(readlink -f "$f")" "$INPUT_DIR/${sample}__$(basename "$f")"
    done
done

activate_env "$ENV_DREP"
echo "[$(date '+%F %T')] dRep dereplicate -sa 0.95"
dRep dereplicate "$OUT_DIR" \
    -g "$INPUT_DIR"/*.fa \
    --genomeInfo "$GINFO" \
    -sa 0.95 \
    -comp 50 -con 10 \
    -p "$THREADS" \
    > "$OUT_DIR/run.log" 2>&1

echo "[$(date '+%F %T')] DONE"
