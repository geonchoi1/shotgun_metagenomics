#!/bin/bash
# === 06 Chromosomal extract = assembly_1kb - plasmid IDs - FREE virus IDs (per sample) ===
# Input:  $PROJECT/00_shared/02_assembly/metaflye/<SAMPLE>/assembly_1kb.fasta
#         $PROJECT/00_shared/03_genomad_virus/<SAMPLE>/*_summary/*_virus.fna
#         $PROJECT/00_shared/04_genomad_plasmid[/filter5_s48]/<SAMPLE>/F12345.ids
# Output: $PROJECT/00_shared/06_chromosomal_extract/<SAMPLE>/chromosomal.fasta
#         $PROJECT/00_shared/06_chromosomal_extract/<SAMPLE>/excluded.ids
#
# Removes: plasmid contigs (F12345) + FREE (whole-contig) viruses only.
# Keeps:   proviruses (integrated prophage lives inside a host chromosome contig, so
#          removing it would drop a multi-Mbp complete chromosome). geNomad headers:
#          a free virus is a plain contig ID; a provirus carries a "|provirus_..." suffix.

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

ASM_BASE=$PROJECT/00_shared/02_assembly/metaflye
VIR_BASE=$PROJECT/00_shared/03_genomad_virus
PLA_BASE=$PROJECT/00_shared/04_genomad_plasmid
OUT_BASE=$PROJECT/00_shared/06_chromosomal_extract
mkdir -p "$OUT_BASE"

echo "[$(date '+%F %T')] chromosomal extract (plasmid + FREE virus removed; proviruses kept)"
for d in "$ASM_BASE"/*/; do
    [ -d "$d" ] || continue
    sample=$(basename "$d")

    asm="$d/assembly_1kb.fasta"
    [ -f "$asm" ] || asm="$d/assembly.fasta"
    [ -f "$asm" ] || asm="$d/scaffolds.fasta"
    [ -f "$asm" ] || { echo "  skip $sample: no assembly"; continue; }

    out_dir="$OUT_BASE/$sample"
    out_fa="$out_dir/chromosomal.fasta"
    excl="$out_dir/excluded.ids"
    if [ -s "$out_fa" ]; then echo "  $sample already done"; continue; fi
    mkdir -p "$out_dir"

    : > "$excl"
    # plasmid F12345 (repo path or community filter5_s48 path)
    pids="$PLA_BASE/$sample/F12345.ids"
    [ -f "$pids" ] || pids="$PLA_BASE/filter5_s48/$sample/F12345.ids"
    [ -f "$pids" ] && cat "$pids" >> "$excl"
    # FREE virus only (drop |provirus rows so host chromosomes are kept)
    vfna=$(ls "$VIR_BASE/$sample"/*_summary/*_virus.fna 2>/dev/null | head -1 || true)
    if [ -n "${vfna:-}" ] && [ -s "$vfna" ]; then
        grep '^>' "$vfna" | sed -e 's/^>//' -e 's/ .*//' | grep -v '|provirus' >> "$excl" || true
    fi
    sort -u "$excl" -o "$excl"

    awk 'BEGIN{while((getline l < "'"$excl"'")>0) drop[l]=1}
         /^>/{id=$1; sub(/^>/,"",id); take = !(id in drop)}
         take{print}' "$asm" > "$out_fa"

    n_asm=$(grep -c '^>' "$asm" || true)
    n_out=$(grep -c '^>' "$out_fa" 2>/dev/null || echo 0)
    n_excl=$(wc -l < "$excl")
    echo "    $sample: $n_asm asm - $n_excl excl = $n_out chromosomal"
done
echo "[$(date '+%F %T')] DONE"
