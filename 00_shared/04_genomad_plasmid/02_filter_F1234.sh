#!/bin/bash
# === 04.02 F1-F4 filter on geNomad plasmid summary TSV ===
# F1: plasmid_score        >= 0.7
# F2: FDR                  <  0.05
# F3: n_hallmarks          >= 1
# F4: conjugation_genes... not used; USCG (universal single-copy genes) <= 1
#
# Input:  $PROJECT/04_plasmid/<SAMPLE>/<SAMPLE>_summary/<SAMPLE>_plasmid_summary.tsv
# Output: $PROJECT/04_plasmid/<SAMPLE>/F1234.ids
#         $PROJECT/04_plasmid/<SAMPLE>/F1234.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

BASE=$PROJECT/04_plasmid

echo "[$(date '+%F %T')] F1-F4 filter (score>=0.7, FDR<0.05, n_hallmarks>=1, USCG<=1)"
for d in "$BASE"/*/; do
    [ -d "$d" ] || continue
    sample=$(basename "$d")
    tsv="$d/${sample}_summary/${sample}_plasmid_summary.tsv"
    [ -s "$tsv" ] || { echo "  skip $sample: no summary tsv"; continue; }

    ids="$d/F1234.ids"
    out_tsv="$d/F1234.tsv"
    if [ -s "$ids" ] && [ -s "$out_tsv" ]; then
        echo "  $sample already done"
        continue
    fi

    echo "  $sample"
    # Determine column indices from header (geNomad columns:
    # seq_name length topology n_genes genetic_code plasmid_score fdr n_hallmarks
    # marker_enrichment_plasmid conjugation_genes amr_genes)
    # USCG-like signal: use n_hallmarks counterpart — geNomad reports n_uscg via
    # marker output in genes.tsv; v1.7+ also exposes 'n_uscg' inside summary
    # if score calibration is on. Otherwise fall back to skipping F4.
    awk -F'\t' -v sample="$sample" '
        NR==1{
            for (i=1;i<=NF;i++) col[$i]=i
            need="seq_name plasmid_score fdr n_hallmarks"
            n=split(need, arr, " ")
            for (k=1;k<=n;k++) if(!(arr[k] in col)){print "ERROR: missing col "arr[k] > "/dev/stderr"; exit 2}
            has_uscg = ("n_uscg" in col) ? 1 : 0
            next
        }
        {
            score = $col["plasmid_score"]+0
            fdr   = ($col["fdr"]=="NA" || $col["fdr"]=="") ? 1 : $col["fdr"]+0
            nh    = $col["n_hallmarks"]+0
            uscg  = has_uscg ? ($col["n_uscg"]+0) : 0
            if (score>=0.7 && fdr<0.05 && nh>=1 && uscg<=1) print
        }
    ' "$tsv" > "$out_tsv"

    head -n 1 "$tsv" > "$out_tsv.tmp" 2>/dev/null || true
    cat "$out_tsv" >> "$out_tsv.tmp"; mv "$out_tsv.tmp" "$out_tsv"

    awk -F'\t' 'NR>1{print $1}' "$out_tsv" > "$ids"
    n=$(wc -l < "$ids")
    echo "    $sample: $n contigs passed F1-F4"
done

echo "[$(date '+%F %T')] DONE"
