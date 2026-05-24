#!/bin/bash
# === 07.10 GTDB-Tk classify_wf on species-dereplicated MAGs ===
# Input:  $PROJECT/07_mag/08_drep_species/dereplicated_genomes/*.fa
# Output: $PROJECT/07_mag/10_gtdbtk/classify/gtdbtk.{bac120,ar53}.summary.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN_DIR=$PROJECT/07_mag/08_drep_species/dereplicated_genomes
OUT_DIR=$PROJECT/07_mag/10_gtdbtk
mkdir -p "$OUT_DIR"

[ -d "$IN_DIR" ] && ls "$IN_DIR"/*.fa 2>/dev/null | grep -q . \
    || { echo "ERROR: no dereplicated MAGs in $IN_DIR" >&2; exit 2; }

if ls "$OUT_DIR"/classify/gtdbtk.*.summary.tsv 2>/dev/null | grep -q .; then
    echo "[$(date '+%F %T')] GTDB-Tk already done"; exit 0
fi

activate_env "$ENV_GTDBTK"
export GTDBTK_DATA_PATH

echo "[$(date '+%F %T')] gtdbtk classify_wf (pplacer_cpus=${PPLACER_CPUS:-32})"
gtdbtk classify_wf \
    --genome_dir "$IN_DIR" \
    --out_dir "$OUT_DIR" \
    --extension fa \
    --cpus "$THREADS_GTDBTK" \
    --pplacer_cpus "${PPLACER_CPUS:-32}" \
    --skip_ani_screen \
    > "$OUT_DIR/run.log" 2>&1

echo "[$(date '+%F %T')] DONE"
