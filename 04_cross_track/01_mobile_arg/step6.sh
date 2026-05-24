#!/bin/bash
# === Step 6: Cross-source BLAST — all modules vs combined PL ∪ MAG ∪ UB FNA ===
# Filters: pident in {95, 99, 100} tiers (post-filter), aln/qlen >= 0.95,
# length >= 1000, self-hit removed, same-origin-region hit removed (in step7).
#
# Output:
#   $PROJECT/cross/mobile_arg/step6/combined_targets.fna  (PL+MAG+UB FNA concat with prefixes)
#   $PROJECT/cross/mobile_arg/step6/combined.{ndb,nhr,nin,nsq,not,ntf,nto}
#   $PROJECT/cross/mobile_arg/step6/modules_combined.fna
#   $PROJECT/cross/mobile_arg/step6/modules_vs_combined.blast.tsv
#   $PROJECT/cross/mobile_arg/step6/modules_vs_combined.filt.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

OUT=$PROJECT/cross/mobile_arg/step6
S5=$PROJECT/cross/mobile_arg/step5
mkdir -p "$OUT"

PL_FNA=$PROJECT/plasmid/04_master_orf/all/master.fna
MAG_FNA=$PROJECT/mag/03_master_orf/all/master.fna
UB_FNA=$PROJECT/unbinned/03_master_orf/all/master.fna

# 1) build combined target FASTA with source prefixes (PL_, MAG_, UB_)
COMB=$OUT/combined_targets.fna
if [ ! -s "$COMB" ]; then
    : > "$COMB"
    for s in plasmid mag unbinned; do
        case $s in
            plasmid)  fa=$PL_FNA;  pfx=PL ;;
            mag)      fa=$MAG_FNA; pfx=MAG ;;
            unbinned) fa=$UB_FNA;  pfx=UB ;;
        esac
        if [ -s "$fa" ]; then
            awk -v p="$pfx" '/^>/{print ">"p"_"substr($1,2); next}{print}' "$fa" >> "$COMB"
        else
            echo "  WARN: $s FNA missing: $fa"
        fi
    done
fi

# 2) makeblastdb (idempotent)
if [ ! -s "$COMB.nhr" ] && [ ! -s "$OUT/combined.nhr" ]; then
    activate_env "$ENV_DIAMOND"  # base env should have blast+
    makeblastdb -in "$COMB" -dbtype nucl -out "$OUT/combined" > "$OUT/makeblastdb.log" 2>&1
fi

# 3) concat module FASTAs (modules_combined.fna)
MOD=$OUT/modules_combined.fna
: > "$MOD"
for s in plasmid mag unbinned; do
    f=$S5/${s}_modules.fna
    [ -s "$f" ] && cat "$f" >> "$MOD"
done
n_mod=$(grep -c '^>' "$MOD" || true)
if [ "$n_mod" -eq 0 ]; then
    echo "[step6] no modules to BLAST — skip"; exit 0
fi

# 4) BLAST modules vs combined
BLAST_OUT=$OUT/modules_vs_combined.blast.tsv
if [ ! -s "$BLAST_OUT" ]; then
    activate_env "$ENV_DIAMOND"
    blastn -query "$MOD" \
           -db "$OUT/combined" \
           -num_threads "$THREADS_BLAST" \
           -evalue 1e-20 \
           -perc_identity 95 \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
           -out "$BLAST_OUT" 2> "$OUT/blastn.log"
fi

# 5) filter: aln/qlen >= 0.95, length >= 1000, pident in {>=95, >=99, ==100}; remove self-hit
FILT=$OUT/modules_vs_combined.filt.tsv
awk -F'\t' '
BEGIN{
    OFS="\t"
    print "qseqid","sseqid","pident","length","qlen","slen","aln_qcov","tier"
}
{
    qseqid=$1; sseqid=$2; pident=$3+0; length=$4+0;
    qlen=$13+0; slen=$14+0
    if(length < 1000) next
    if(qlen == 0) next
    qcov = length / qlen
    if(qcov < 0.95) next
    # remove self-hit: parse module id "src|contig|start_end|ARG:..." vs target "PFX_contig"
    n=split(qseqid,p,"|")
    qcontig = p[2]
    # target sseqid like "PL_contigname" — split on first _
    pos = index(sseqid,"_")
    if(pos>0){
        tcontig = substr(sseqid,pos+1)
        if(tcontig == qcontig) next
    }
    if(pident == 100) tier="clonal"
    else if(pident >= 99) tier="very_recent"
    else if(pident >= 95) tier="recent"
    else next
    print qseqid,sseqid,pident,length,qlen,slen,qcov,tier
}' "$BLAST_OUT" > "$FILT"

n=$(awk 'NR>1' "$FILT" | wc -l)
echo "[step6] DONE — modules=$n_mod filtered_hits=$n -> $FILT"
