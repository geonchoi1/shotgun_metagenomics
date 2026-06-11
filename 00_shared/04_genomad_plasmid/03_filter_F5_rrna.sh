#!/bin/bash
# === 04.03 F5 rRNA filter (Infernal cmscan vs RFAM_CM) ===
# Removes any F1234 contig that has an rRNA hit (5S/16S/23S, bac+arc).
#
# Input:  $PROJECT/00_shared/04_genomad_plasmid/<SAMPLE>/F1234.ids
#         $PROJECT/00_shared/04_genomad_plasmid/<SAMPLE>/<SAMPLE>_summary/<SAMPLE>_plasmid.fna
# Output: $PROJECT/00_shared/04_genomad_plasmid/<SAMPLE>/F12345.ids
#         $PROJECT/00_shared/04_genomad_plasmid/<SAMPLE>/F12345.fna
#         $PROJECT/00_shared/04_genomad_plasmid/all_putative.fna   (aggregated)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

BASE=$PROJECT/00_shared/04_genomad_plasmid
[ -f "$RFAM_CM" ] || { echo "ERROR: RFAM_CM not found: $RFAM_CM" >&2; exit 2; }

activate_env "$ENV_INFERNAL"

echo "[$(date '+%F %T')] F5 cmscan vs $RFAM_CM"
for d in "$BASE"/*/; do
    [ -d "$d" ] || continue
    sample=$(basename "$d")
    ids="$d/F1234.ids"
    fna_in="$d/${sample}_summary/${sample}_plasmid.fna"
    [ -s "$ids" ] || { echo "  skip $sample: no F1234.ids"; continue; }
    [ -s "$fna_in" ] || { echo "  skip $sample: no plasmid.fna"; continue; }

    out_ids="$d/F12345.ids"
    out_fna="$d/F12345.fna"
    if [ -s "$out_ids" ] && [ -s "$out_fna" ]; then
        echo "  $sample already done"
        continue
    fi

    f1234_fna="$d/F1234.fna"
    # Extract F1234 contigs (samtools is in many envs; fall back to awk)
    awk 'BEGIN{while((getline l < "'"$ids"'")>0) keep[l]=1}
         /^>/{id=$1; sub(/^>/,"",id); take = (id in keep)}
         take{print}' "$fna_in" > "$f1234_fna"

    cm_tbl="$d/F5_cmscan.tblout"
    cmscan \
        --cpu "$THREADS" \
        --rfam --cut_ga --nohmmonly \
        --tblout "$cm_tbl" \
        "$RFAM_CM" "$f1234_fna" \
        > "$d/F5_cmscan.log" 2>&1 || true

    # contigs with any rRNA hit -> reject; query name is column 3 in tblout
    rrna_ids="$d/F5_rrna_hits.ids"
    awk '!/^#/{print $3}' "$cm_tbl" | sort -u > "$rrna_ids" || true

    # F12345 = F1234 minus rrna_hits
    grep -vx -F -f "$rrna_ids" "$ids" > "$out_ids" || true

    awk 'BEGIN{while((getline l < "'"$out_ids"'")>0) keep[l]=1}
         /^>/{id=$1; sub(/^>/,"",id); take = (id in keep)}
         take{print}' "$fna_in" > "$out_fna"

    n_in=$(wc -l < "$ids")
    n_out=$(wc -l < "$out_ids")
    echo "    $sample: $n_in -> $n_out after F5"
done

echo "[$(date '+%F %T')] aggregate F12345 fasta"
agg="$BASE/all_putative.fna"
: > "$agg"
for d in "$BASE"/*/; do
    [ -d "$d" ] || continue
    sample=$(basename "$d")
    f="$d/F12345.fna"
    [ -s "$f" ] || continue
    awk -v s="$sample" '/^>/{sub(/^>/,">"s"|"); print; next} {print}' "$f" >> "$agg"
done

n=$(grep -c '^>' "$agg" 2>/dev/null || echo 0)
echo "[$(date '+%F %T')] DONE — $n putative plasmid contigs -> $agg"
