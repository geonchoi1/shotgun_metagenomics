#!/bin/bash
# === 08 Split chromosomal contigs into binned vs unbinned (+ circular topology) ===
# Binned   = chromosomal contig captured in ANY DAS_Tool MAG (any sample)
# Unbinned = the rest
# Circular topology is taken from metaFlye assembly_info (col 4 "circ." = Y/N) so the
# MAG (02_mag_track) and unbinned (03_unbinned_track) tracks can run Bakta with
# --complete on circular replicons.
#
# Input:  $PROJECT/00_shared/02_assembly/metaflye/<SAMPLE>/assembly_info.txt
#         $PROJECT/00_shared/06_chromosomal_extract/<SAMPLE>/chromosomal.fasta
#         $PROJECT/00_shared/07_mag_production/04_dastool/<SAMPLE>/_DASTool_bins/*.fa
#         $PROJECT/00_shared/07_mag_production/08_drep_species/dereplicated_genomes/*.fa
# Output: $PROJECT/00_shared/08_split_binned_unbinned/{binned,unbinned}/all.fna
#         $PROJECT/00_shared/08_split_binned_unbinned/{binned,unbinned}/circ.ids
#         $PROJECT/00_shared/07_mag_production/circ_mag.tsv   (MAG<TAB>circ|frag)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

ASM_BASE=$PROJECT/00_shared/02_assembly/metaflye
CHR_BASE=$PROJECT/00_shared/06_chromosomal_extract
DAS_BASE=$PROJECT/00_shared/07_mag_production/04_dastool
DREP_DIR=$PROJECT/00_shared/07_mag_production/08_drep_species/dereplicated_genomes
OUT_BASE=$PROJECT/00_shared/08_split_binned_unbinned
mkdir -p "$OUT_BASE/binned" "$OUT_BASE/unbinned"

BINNED_FA=$OUT_BASE/binned/all.fna
UNBINNED_FA=$OUT_BASE/unbinned/all.fna
BIN_CIRC=$OUT_BASE/binned/circ.ids
UNB_CIRC=$OUT_BASE/unbinned/circ.ids
: > "$BINNED_FA"; : > "$UNBINNED_FA"; : > "$BIN_CIRC"; : > "$UNB_CIRC"

echo "[$(date '+%F %T')] split per sample (binned/unbinned + circular flag)"
for d in "$CHR_BASE"/*/; do
    [ -d "$d" ] || continue
    sample=$(basename "$d")
    chr="$d/chromosomal.fasta"
    [ -s "$chr" ] || { echo "  skip $sample"; continue; }

    # contigs captured in any MAG for this sample
    bin_ids="$OUT_BASE/binned/${sample}.ids"
    : > "$bin_ids"
    das_dir="$DAS_BASE/$sample/_DASTool_bins"
    if [ -d "$das_dir" ]; then
        for f in "$das_dir"/*.fa; do
            [ -f "$f" ] || continue
            grep '^>' "$f" | sed -e 's/^>//' -e 's/ .*//' >> "$bin_ids"
        done
    fi
    sort -u "$bin_ids" -o "$bin_ids"

    # circular contig IDs for this sample (assembly_info col4 == Y)
    info="$ASM_BASE/$sample/assembly_info.txt"
    circ_set="$OUT_BASE/.circ_${sample}.ids"
    if [ -s "$info" ]; then awk -F'\t' 'NR>1 && $4=="Y"{print $1}' "$info" > "$circ_set"; else : > "$circ_set"; fi

    # append (>>) binned/unbinned FASTAs (prefixed) and circular IDs
    awk -v s="$sample" -v ids="$bin_ids" -v circ="$circ_set" \
        -v bo="$BINNED_FA" -v uo="$UNBINNED_FA" -v bc="$BIN_CIRC" -v uc="$UNB_CIRC" '
        BEGIN{ while((getline l < ids)>0) keep[l]=1; while((getline c < circ)>0) iscirc[c]=1 }
        /^>/{
            id=$1; sub(/^>/,"",id)
            binned=(id in keep)
            out = binned ? bo : uo
            print ">" s "|" id >> out
            if(id in iscirc){ if(binned) print s"|"id >> bc; else print s"|"id >> uc }
            next
        }
        { print >> out }
    ' "$chr"
    rm -f "$circ_set"
    echo "    $sample done (binned=$(grep -c '^>' $BINNED_FA), unbinned=$(grep -c '^>' $UNBINNED_FA))"
done

# Per-MAG topology: a MAG is "circ" iff it is a SINGLE circular contig (complete replicon).
# MAG files are named <sample>__<bin>.fa; contigs inside are unprefixed (e.g. contig_5).
CIRC_MAG=$PROJECT/00_shared/07_mag_production/circ_mag.tsv
if [ -d "$DREP_DIR" ] && ls "$DREP_DIR"/*.fa >/dev/null 2>&1; then
    : > "$CIRC_MAG"
    for mag in "$DREP_DIR"/*.fa; do
        [ -f "$mag" ] || continue
        name=$(basename "$mag" .fa)
        sample=${name%%__*}
        info="$ASM_BASE/$sample/assembly_info.txt"
        n_ctg=$(grep -c '^>' "$mag")
        label=frag
        if [ "$n_ctg" -eq 1 ] && [ -s "$info" ]; then
            ctg=$(grep '^>' "$mag" | head -1 | sed -e 's/^>//' -e 's/ .*//')
            awk -F'\t' -v c="$ctg" 'NR>1 && $1==c && $4=="Y"{f=1} END{exit !f}' "$info" && label=circ
        fi
        printf '%s\t%s\n' "$name" "$label" >> "$CIRC_MAG"
    done
    echo "  circ_mag.tsv: $(awk -F'\t' '$2=="circ"' "$CIRC_MAG" | wc -l) circ / $(wc -l < "$CIRC_MAG") MAGs"
fi

echo "[$(date '+%F %T')] DONE"
echo "  binned   : $(grep -c '^>' $BINNED_FA) contigs ($(wc -l < $BIN_CIRC) circular)"
echo "  unbinned : $(grep -c '^>' $UNBINNED_FA) contigs ($(wc -l < $UNB_CIRC) circular)"
