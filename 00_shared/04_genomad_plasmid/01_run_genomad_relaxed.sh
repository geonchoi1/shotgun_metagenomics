#!/bin/bash
# === 04.01 geNomad plasmid identification (relaxed -s 4.8 + score calibration) ===
# Input:  $PROJECT/00_shared/02_assembly/metaflye/<SAMPLE>/assembly_1kb.fasta | scaffolds.fasta
# Output: $PROJECT/00_shared/04_genomad_plasmid/<SAMPLE>/<SAMPLE>_summary/<SAMPLE>_plasmid_summary.tsv
#         $PROJECT/00_shared/04_genomad_plasmid/<SAMPLE>/<SAMPLE>_summary/<SAMPLE>_plasmid.fna

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

ASM_BASE=$PROJECT/00_shared/02_assembly/metaflye
OUT_BASE=$PROJECT/00_shared/04_genomad_plasmid
mkdir -p "$OUT_BASE"

ENV_GENOMAD=${ENV_GENOMAD:-genomad}
activate_env "$ENV_GENOMAD"

echo "[$(date '+%F %T')] geNomad end-to-end (plasmid, -s 4.8 --relaxed --enable-score-calibration)"
for d in "$ASM_BASE"/*/; do
    [ -d "$d" ] || continue
    sample=$(basename "$d")

    asm="$d/assembly_1kb.fasta"
    [ -f "$asm" ] || asm="$d/assembly_1kb.fasta"
    [ -f "$asm" ] || { echo "  skip $sample: no assembly"; continue; }

    out_dir="$OUT_BASE/$sample"
    summary_tsv="$out_dir/${sample}_summary/${sample}_plasmid_summary.tsv"
    if [ -s "$summary_tsv" ]; then
        echo "  $sample already done"
        continue
    fi

    mkdir -p "$out_dir"
    echo "  $sample"
    genomad end-to-end \
        -s 4.8 --relaxed --enable-score-calibration \
        "$asm" \
        "$out_dir" \
        "$GENOMAD_DB" \
        --threads "$THREADS" \
        --cleanup \
        > "$out_dir/run.log" 2>&1 || true
done

echo "[$(date '+%F %T')] DONE"
