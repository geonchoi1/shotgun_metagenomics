#!/bin/bash
# === 07.02 Per-contig depth via jgi_summarize_bam_contig_depths (metabat2 env) ===
# Input:  $PROJECT/00_shared/07_mag_production/01_mapping/<SAMPLE>.bam
# Output: $PROJECT/00_shared/07_mag_production/02_depth/<SAMPLE>_depth.txt

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

BAM_DIR=$PROJECT/00_shared/07_mag_production/01_mapping
OUT_DIR=$PROJECT/00_shared/07_mag_production/02_depth
mkdir -p "$OUT_DIR"

activate_env "$ENV_METABAT2"

echo "[$(date '+%F %T')] jgi_summarize_bam_contig_depths"
for bam in "$BAM_DIR"/*.bam; do
    [ -f "$bam" ] || continue
    sample=$(basename "$bam" .bam)
    out="$OUT_DIR/${sample}_depth.txt"
    if [ -s "$out" ]; then
        echo "  $sample already done"; continue
    fi
    echo "  $sample"
    jgi_summarize_bam_contig_depths --outputDepth "$out" "$bam"
done

echo "[$(date '+%F %T')] DONE"
