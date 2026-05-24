#!/bin/bash
# === 03 Master ORF (UB): wrap Bakta output as unified master files ===
# Input:  $PROJECT/unbinned/02_bakta/unbinned/unbinned.{faa,ffn,gff3,fna}
# Output: $PROJECT/unbinned/03_master_orf/all/master.{faa,ffn,fna,gff}
#         $PROJECT/unbinned/03_master_orf/all/orf2contig.tsv
# Notes: UB has no circ/frag distinction unless 01_raw_fasta produced it; we always emit "all/".

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

BAK=$PROJECT/unbinned/02_bakta/unbinned
OUT=$PROJECT/unbinned/03_master_orf/all
mkdir -p "$OUT"

for ext in faa ffn fna gff3; do
    [ -s "$BAK/unbinned.$ext" ] || { echo "ERROR: missing $BAK/unbinned.$ext" >&2; exit 1; }
done

echo "[$(date '+%F %T')] wrap Bakta output as master ORF (UB)"
ln -sf "$BAK/unbinned.faa"  "$OUT/master.faa"
ln -sf "$BAK/unbinned.ffn"  "$OUT/master.ffn"
ln -sf "$BAK/unbinned.fna"  "$OUT/master.fna"
ln -sf "$BAK/unbinned.gff3" "$OUT/master.gff"

# orf2contig.tsv: ORF_id<TAB>contig_id (parsed from GFF CDS lines)
orf2c="$OUT/orf2contig.tsv"
awk -F'\t' '$3=="CDS" {
    contig=$1
    n=split($9,a,";")
    id=""
    for(i=1;i<=n;i++) if(a[i] ~ /^ID=/){ id=substr(a[i],4); break }
    if(id!="") print id"\t"contig
}' "$BAK/unbinned.gff3" > "$orf2c"

n_orf=$(grep -c '^>' "$OUT/master.faa" || true)
n_ctg=$(cut -f2 "$orf2c" | sort -u | wc -l)
echo "[$(date '+%F %T')] DONE — ORFs=$n_orf on contigs=$n_ctg"
