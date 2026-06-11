#!/bin/bash
# === 03 Master ORF (UB): combine circ + frag Bakta outputs ===
# Input:  $PROJECT/03_unbinned_track/02_bakta/{circ/unbinned_circ,frag/unbinned_frag}.{faa,ffn,gff3}
# Output: $PROJECT/03_unbinned_track/03_master_orf/{circ,frag,all}/master.{faa,ffn,gff}
#         $PROJECT/03_unbinned_track/03_master_orf/all/orf2contig.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

BAK=$PROJECT/03_unbinned_track/02_bakta
OUT=$PROJECT/03_unbinned_track/03_master_orf
mkdir -p "$OUT"/{circ,frag,all}

# circ subset (may be empty if no circ contigs)
if [ -s "$BAK/circ/unbinned_circ.faa" ]; then
    ln -sf "$BAK/circ/unbinned_circ.faa"  "$OUT/circ/master.faa"
    ln -sf "$BAK/circ/unbinned_circ.ffn"  "$OUT/circ/master.ffn"
    ln -sf "$BAK/circ/unbinned_circ.gff3" "$OUT/circ/master.gff"
fi

# frag subset
if [ -s "$BAK/frag/unbinned_frag.faa" ]; then
    ln -sf "$BAK/frag/unbinned_frag.faa"  "$OUT/frag/master.faa"
    ln -sf "$BAK/frag/unbinned_frag.ffn"  "$OUT/frag/master.ffn"
    ln -sf "$BAK/frag/unbinned_frag.gff3" "$OUT/frag/master.gff"
fi

# combined all/
echo "[$(date '+%F %T')] Building combined master ORF (circ + frag)"
> "$OUT/all/master.faa"
> "$OUT/all/master.ffn"
> "$OUT/all/master.gff"
first=1
for src in "$BAK/circ/unbinned_circ" "$BAK/frag/unbinned_frag"; do
    [ -s "${src}.faa" ] || continue
    cat "${src}.faa" >> "$OUT/all/master.faa"
    cat "${src}.ffn" >> "$OUT/all/master.ffn"
    if [ "$first" = "1" ]; then
        cp "${src}.gff3" "$OUT/all/master.gff"
        first=0
    else
        grep -v "^#" "${src}.gff3" >> "$OUT/all/master.gff"
    fi
done

# orf2contig.tsv (parse locus_tag → contig from GFF CDS lines)
awk -F'\t' '$3=="CDS" {
    contig=$1
    n=split($9,a,";")
    id=""
    for(i=1;i<=n;i++) if(a[i] ~ /^locus_tag=/){ id=substr(a[i],11); break }
    if(id=="") for(i=1;i<=n;i++) if(a[i] ~ /^ID=/){ id=substr(a[i],4); break }
    if(id!="") print id"\t"contig
}' "$OUT/all/master.gff" > "$OUT/all/orf2contig.tsv"

n_orf=$(grep -c '^>' "$OUT/all/master.faa" || true)
n_ctg=$(cut -f2 "$OUT/all/orf2contig.tsv" | sort -u | wc -l)
echo "[$(date '+%F %T')] DONE — ORFs=$n_orf on contigs=$n_ctg"
