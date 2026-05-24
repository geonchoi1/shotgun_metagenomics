#!/bin/bash
# === 03 geNomad virus identification (default params) ===
# Input:  $PROJECT/02_assembly/<SAMPLE>/assembly.fasta (HiFi) | scaffolds.fasta (Illumina)
# Output: $PROJECT/03_virus/<SAMPLE>/<SAMPLE>_summary/<SAMPLE>_virus.fna
#         $PROJECT/03_virus/all_virus.fna   (aggregated across samples)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

ASM_BASE=$PROJECT/02_assembly
OUT_BASE=$PROJECT/03_virus
mkdir -p "$OUT_BASE"

# geNomad env — fallback if not set in config.sh
ENV_GENOMAD=${ENV_GENOMAD:-genomad}
activate_env "$ENV_GENOMAD"

echo "[$(date '+%F %T')] geNomad end-to-end (virus, default params)"
for d in "$ASM_BASE"/*/; do
    [ -d "$d" ] || continue
    sample=$(basename "$d")

    asm="$d/assembly.fasta"
    [ -f "$asm" ] || asm="$d/scaffolds.fasta"
    [ -f "$asm" ] || { echo "  skip $sample: no assembly"; continue; }

    out_dir="$OUT_BASE/$sample"
    virus_fna="$out_dir/${sample}_summary/${sample}_virus.fna"
    if [ -s "$virus_fna" ]; then
        echo "  $sample already done"
        continue
    fi

    mkdir -p "$out_dir"
    echo "  $sample"
    genomad end-to-end \
        "$asm" \
        "$out_dir" \
        "$GENOMAD_DB" \
        --threads "$THREADS" \
        --cleanup \
        > "$out_dir/run.log" 2>&1 || true
done

echo "[$(date '+%F %T')] aggregate virus FNAs"
agg="$OUT_BASE/all_virus.fna"
: > "$agg"
for d in "$OUT_BASE"/*/; do
    [ -d "$d" ] || continue
    sample=$(basename "$d")
    fna="$d/${sample}_summary/${sample}_virus.fna"
    [ -s "$fna" ] || continue
    # Prefix contig IDs with sample to avoid collisions across samples
    awk -v s="$sample" '/^>/{sub(/^>/,">"s"|"); print; next} {print}' "$fna" >> "$agg"
done

n=$(grep -c '^>' "$agg" 2>/dev/null || echo 0)
echo "[$(date '+%F %T')] DONE — $n virus contigs aggregated -> $agg"
