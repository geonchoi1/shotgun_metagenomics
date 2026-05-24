#!/bin/bash
# === 06 Chromosomal extract = assembly - plasmid IDs - virus IDs (per sample) ===
# Input:  $PROJECT/02_assembly/<SAMPLE>/assembly.fasta | scaffolds.fasta
#         $PROJECT/03_virus/<SAMPLE>/<SAMPLE>_summary/<SAMPLE>_virus.fna
#         $PROJECT/04_plasmid/<SAMPLE>/F12345.ids
# Output: $PROJECT/06_chromosomal/<SAMPLE>/chromosomal.fasta
#         $PROJECT/06_chromosomal/<SAMPLE>/excluded.ids

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

ASM_BASE=$PROJECT/02_assembly
VIR_BASE=$PROJECT/03_virus
PLA_BASE=$PROJECT/04_plasmid
OUT_BASE=$PROJECT/06_chromosomal
mkdir -p "$OUT_BASE"

echo "[$(date '+%F %T')] chromosomal extract per sample"
for d in "$ASM_BASE"/*/; do
    [ -d "$d" ] || continue
    sample=$(basename "$d")

    asm="$d/assembly.fasta"
    [ -f "$asm" ] || asm="$d/scaffolds.fasta"
    [ -f "$asm" ] || { echo "  skip $sample: no assembly"; continue; }

    out_dir="$OUT_BASE/$sample"
    out_fa="$out_dir/chromosomal.fasta"
    excl="$out_dir/excluded.ids"
    if [ -s "$out_fa" ]; then
        echo "  $sample already done"
        continue
    fi
    mkdir -p "$out_dir"

    # Collect plasmid + virus contig IDs (these are raw assembly contig names)
    : > "$excl"
    [ -f "$PLA_BASE/$sample/F12345.ids" ] && cat "$PLA_BASE/$sample/F12345.ids" >> "$excl"
    vfna="$VIR_BASE/$sample/${sample}_summary/${sample}_virus.fna"
    if [ -s "$vfna" ]; then
        # geNomad provirus IDs may have |provirus suffix — strip to match assembly contigs
        grep '^>' "$vfna" | sed -e 's/^>//' -e 's/ .*//' -e 's/|provirus.*//' >> "$excl"
    fi
    sort -u "$excl" -o "$excl"

    # Filter assembly: keep contigs NOT in excl
    awk 'BEGIN{while((getline l < "'"$excl"'")>0) drop[l]=1}
         /^>/{id=$1; sub(/^>/,"",id); take = !(id in drop)}
         take{print}' "$asm" > "$out_fa"

    n_asm=$(grep -c '^>' "$asm")
    n_out=$(grep -c '^>' "$out_fa" 2>/dev/null || echo 0)
    n_excl=$(wc -l < "$excl")
    echo "    $sample: $n_asm asm - $n_excl excl = $n_out chromosomal"
done

echo "[$(date '+%F %T')] DONE"
