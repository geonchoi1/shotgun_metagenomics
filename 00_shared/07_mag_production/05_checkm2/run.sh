#!/bin/bash
# === 07.05 CheckM2 on DAS_Tool consensus bins (all samples pooled) ===
# Input:  $PROJECT/00_shared/07_mag_production/04_dastool/<SAMPLE>/_DASTool_bins/*.fa
# Output: $PROJECT/00_shared/07_mag_production/05_checkm2/quality_report.tsv
#         $PROJECT/00_shared/07_mag_production/05_checkm2/genomeInfo.csv (drep-ready)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

DAS_BASE=$PROJECT/00_shared/07_mag_production/04_dastool
OUT_DIR=$PROJECT/00_shared/07_mag_production/05_checkm2
# input_bins must live OUTSIDE OUT_DIR — checkm2 predict --force wipes the output dir.
INPUT_DIR=${OUT_DIR}_input_bins
mkdir -p "$OUT_DIR" "$INPUT_DIR"

REPORT=$OUT_DIR/quality_report.tsv
if [ -s "$REPORT" ]; then
    echo "[$(date '+%F %T')] CheckM2 already done"
else
    echo "[$(date '+%F %T')] gather DAS_Tool bins -> $INPUT_DIR"
    for d in "$DAS_BASE"/*/_DASTool_bins/; do
        [ -d "$d" ] || continue
        sample=$(basename "$(dirname "$d")")
        for f in "$d"/*.fa; do
            [ -f "$f" ] || continue
            ln -sf "$(readlink -f "$f")" "$INPUT_DIR/${sample}__$(basename "$f")"
        done
    done

    activate_env "$ENV_CHECKM2"
    echo "[$(date '+%F %T')] CheckM2 predict"
    checkm2 predict \
        --threads "$THREADS" \
        --input "$INPUT_DIR" \
        --output-directory "$OUT_DIR" \
        --database_path "$CHECKM2_DB" \
        --force \
        > "$OUT_DIR/run.log" 2>&1
fi

# Build dRep genomeInfo.csv: genome,completeness,contamination
GINFO=$OUT_DIR/genomeInfo.csv
if [ ! -s "$GINFO" ] && [ -s "$REPORT" ]; then
    awk -F'\t' 'NR==1{
        for(i=1;i<=NF;i++) c[$i]=i; print "genome,completeness,contamination"; next
    }
    { printf "%s.fa,%s,%s\n", $c["Name"], $c["Completeness"], $c["Contamination"] }' "$REPORT" > "$GINFO"
fi

echo "[$(date '+%F %T')] DONE"
