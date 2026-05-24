#!/bin/bash
# === 22 MMseqs2 easy-taxonomy LCA on UB contigs (Priest 2025 method) ===
# Input:  $PROJECT/unbinned/01_raw_fasta/all/unbinned.fna
# DB:     $MMSEQS_GTDB_DB (GTDB amino-acid db)
# Output: $PROJECT/unbinned/22_mmseqs2_lca_taxonomy/all/
#           lca_output_lca.tsv      raw MMseqs LCA per contig
#           lca_output_report       Kraken-like report
#           lca_output_tophit_aln   top hit alignments
#           contig_taxonomy.tsv     contig<TAB>kingdom..species (7 ranks)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

IN=$PROJECT/unbinned/01_raw_fasta/all/unbinned.fna
OUT=$PROJECT/unbinned/22_mmseqs2_lca_taxonomy/all
TMP=$OUT/tmp
mkdir -p "$OUT" "$TMP"

[ -s "$IN" ] || { echo "ERROR: missing $IN" >&2; exit 1; }
[ -s "${MMSEQS_GTDB_DB}" ] || [ -s "${MMSEQS_GTDB_DB}.dbtype" ] || \
    { echo "ERROR: missing MMseqs DB $MMSEQS_GTDB_DB" >&2; exit 1; }

LCA_PREFIX=$OUT/lca_output

if [ -s "${LCA_PREFIX}_lca.tsv" ] && [ -s "$OUT/contig_taxonomy.tsv" ]; then
    echo "[$(date '+%F %T')] UB mmseqs2 LCA already done — skip"; exit 0
fi

activate_env "$ENV_MMSEQS"

echo "[$(date '+%F %T')] mmseqs easy-taxonomy on UB contigs"
mmseqs easy-taxonomy \
    "$IN" \
    "$MMSEQS_GTDB_DB" \
    "$LCA_PREFIX" \
    "$TMP" \
    --tax-lineage 1 \
    --threads "$THREADS" \
    > "$OUT/mmseqs.log" 2>&1

rm -rf "$TMP"

# Parse _lca.tsv (cols: query, taxid, rank, sciname, lineage)
# Lineage example: d_Bacteria;p_Pseudomonadota;c_Gammaproteobacteria;o_Pseudomonadales;...
echo "[$(date '+%F %T')] parse contig → kingdom..species"
out_tax="$OUT/contig_taxonomy.tsv"
{
    printf "contig\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"
    awk -F'\t' '
    {
        contig=$1; lineage=$NF
        split("",rank,"")
        n=split(lineage,parts,";")
        for(i=1;i<=n;i++){
            t=parts[i]
            if(t ~ /^d_/) rank["k"]=substr(t,3)
            else if(t ~ /^p_/) rank["p"]=substr(t,3)
            else if(t ~ /^c_/) rank["c"]=substr(t,3)
            else if(t ~ /^o_/) rank["o"]=substr(t,3)
            else if(t ~ /^f_/) rank["f"]=substr(t,3)
            else if(t ~ /^g_/) rank["g"]=substr(t,3)
            else if(t ~ /^s_/) rank["s"]=substr(t,3)
        }
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", contig,
            (("k" in rank)?rank["k"]:"NA"),
            (("p" in rank)?rank["p"]:"NA"),
            (("c" in rank)?rank["c"]:"NA"),
            (("o" in rank)?rank["o"]:"NA"),
            (("f" in rank)?rank["f"]:"NA"),
            (("g" in rank)?rank["g"]:"NA"),
            (("s" in rank)?rank["s"]:"NA")
    }' "${LCA_PREFIX}_lca.tsv"
} > "$out_tax"

n=$(awk 'NR>1' "$out_tax" | wc -l)
echo "[$(date '+%F %T')] DONE — $n contigs taxonomy-assigned"
