#!/bin/bash
# === 08 Split chromosomal contigs into binned vs unbinned (+ circular topology) ===
# Binned   = chromosomal contig captured in ANY DAS_Tool MAG (any sample)
# Unbinned = the rest
# Circular topology is read from the MASTER contig_topology.tsv (auto-generated from
# metaFlye assembly_info, the same single source 05 uses for plasmids) so the MAG
# (02_mag_track) and unbinned (03_unbinned_track) tracks can run Bakta --complete on
# circular replicons.
#
# Input:  $PROJECT/00_shared/02_assembly/contig_topology.tsv  (auto; <sample>|<contig> TAB circ|frag)
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
MASTER=$PROJECT/00_shared/02_assembly/contig_topology.tsv
mkdir -p "$OUT_BASE/binned" "$OUT_BASE/unbinned"

# Auto-generate the master topology from metaFlye assembly_info if missing.
if [ ! -s "$MASTER" ]; then
    echo "[$(date '+%F %T')] generating master contig_topology.tsv from assembly_info"
    : > "$MASTER"
    for info in "$ASM_BASE"/*/assembly_info.txt; do
        [ -s "$info" ] || continue
        s=$(basename "$(dirname "$info")")
        awk -F'\t' -v s="$s" 'NR>1{print s"|"$1"\t"($4=="Y"?"circ":"frag")}' "$info" >> "$MASTER"
    done
fi
# circular <sample>|<contig> keys (master is the single source of truth)
CIRC_MASTER=$OUT_BASE/.circ_master.ids
awk -F'\t' '$2=="circ"{print $1}' "$MASTER" > "$CIRC_MASTER" 2>/dev/null || : > "$CIRC_MASTER"

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

    # append (>>) binned/unbinned FASTAs (prefixed) and circular IDs (from master)
    awk -v s="$sample" -v ids="$bin_ids" -v circ="$CIRC_MASTER" \
        -v bo="$BINNED_FA" -v uo="$UNBINNED_FA" -v bc="$BIN_CIRC" -v uc="$UNB_CIRC" '
        BEGIN{ while((getline l < ids)>0) keep[l]=1; while((getline c < circ)>0) iscirc[c]=1 }
        /^>/{
            id=$1; sub(/^>/,"",id)
            binned=(id in keep)
            out = binned ? bo : uo
            key = s "|" id
            print ">" key >> out
            if(key in iscirc){ if(binned) print key >> bc; else print key >> uc }
            next
        }
        { print >> out }
    ' "$chr"
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
        n_ctg=$(grep -c '^>' "$mag")
        label=frag
        if [ "$n_ctg" -eq 1 ]; then
            ctg=$(grep '^>' "$mag" | head -1 | sed -e 's/^>//' -e 's/ .*//')
            grep -qxF "${sample}|${ctg}" "$CIRC_MASTER" && label=circ
        fi
        printf '%s\t%s\n' "$name" "$label" >> "$CIRC_MAG"
    done
    echo "  circ_mag.tsv: $(awk -F'\t' '$2=="circ"' "$CIRC_MAG" | wc -l) circ / $(wc -l < "$CIRC_MAG") MAGs"
fi
rm -f "$CIRC_MASTER"

echo "[$(date '+%F %T')] DONE"
echo "  binned   : $(grep -c '^>' $BINNED_FA) contigs ($(wc -l < $BIN_CIRC) circular)"
echo "  unbinned : $(grep -c '^>' $UNBINNED_FA) contigs ($(wc -l < $UNB_CIRC) circular)"
