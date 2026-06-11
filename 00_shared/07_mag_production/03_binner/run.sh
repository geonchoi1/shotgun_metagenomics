#!/bin/bash
# === 07.03 Binner dispatcher ===
# HiFi      -> MetaBinner + MetaDecoder + SemiBin2 (--sequencing-type long_read)
# Illumina  -> MetaBAT2 (default)     + SemiBin2 (--sequencing-type short_read)
#
# Inputs:  $PROJECT/00_shared/06_chromosomal_extract/<SAMPLE>/chromosomal.fasta
#          $PROJECT/00_shared/07_mag_production/02_depth/<SAMPLE>_depth.txt    (jgi format)
#          $PROJECT/00_shared/07_mag_production/01_mapping/<SAMPLE>.bam
# Outputs: $PROJECT/00_shared/07_mag_production/03_binner/{metabinner,metadecoder,semibin,metabat2}/<SAMPLE>/

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}
: ${READ_TYPE:?ERROR: export READ_TYPE=hifi or illumina}

REF_BASE=$PROJECT/00_shared/06_chromosomal_extract
DEPTH_DIR=$PROJECT/00_shared/07_mag_production/02_depth
BAM_DIR=$PROJECT/00_shared/07_mag_production/01_mapping
BIN_BASE=$PROJECT/00_shared/07_mag_production/03_binner
mkdir -p "$BIN_BASE"

list_samples() {
    for d in "$REF_BASE"/*/; do [ -d "$d" ] && basename "$d"; done
}

# Convert jgi depth (MetaBAT2 format) -> MetaBinner coverage_profile.tsv
# jgi cols: contigName contigLen totalAvgDepth <SAMPLE>.bam <SAMPLE>.bam-var ...
make_metabinner_cov() {
    local jgi=$1 out=$2
    awk 'BEGIN{OFS="\t"}
        NR==1{
            # build header: contigName + every non-variance sample column
            printf "contigName"
            for(i=4;i<=NF;i+=2) printf "\t%s", $i
            printf "\n"; next
        }
        {
            printf "%s", $1
            for(i=4;i<=NF;i+=2) printf "\t%s", $i
            printf "\n"
        }' "$jgi" > "$out"
}

run_metabinner() {
    local sample=$1
    local ref="$REF_BASE/$sample/chromosomal.fasta"
    local depth="$DEPTH_DIR/${sample}_depth.txt"
    local out="$BIN_BASE/metabinner/$sample"
    [ -d "$out/metabinner_res" ] && [ "$(ls -A "$out/metabinner_res" 2>/dev/null)" ] && { echo "    metabinner $sample done"; return 0; }
    mkdir -p "$out"
    # kmer profile (>=1kb) via MetaBinner's gen_kmer.py — generated FIRST so the coverage
    # below can be built on exactly the same contig set/order.
    local kmer="$out/kmer_4_f1000.csv"
    local mb_dir
    mb_dir=$(dirname "$(dirname "$(command -v run_metabinner.sh 2>/dev/null || echo /dev/null)")")
    if [ ! -s "$kmer" ] && [ -f "$mb_dir/scripts/gen_kmer.py" ]; then
        python "$mb_dir/scripts/gen_kmer.py" "$ref" 1000 4
        mv "$(dirname "$ref")"/$(basename "$ref" .fasta)_kmer_4_f1000.csv "$kmer" 2>/dev/null || true
    fi

    # Coverage profile built FROM the kmer contig list (identical set AND order), with depth
    # from the jgi totalAvgDepth column (0 if a contig is absent from the depth file).
    # MetaBinner requires coverage and kmer on exactly the same contigs in the same order;
    # a jgi-ordered/filtered coverage triggers a split_hhbins.py KeyError, so build coverage
    # directly from the kmer file.
    local cov="$out/coverage_profile.tsv"
    if [ -s "$kmer" ] && [ -s "$depth" ]; then
        tail -n +2 "$kmer" | cut -d',' -f1 > "$out/.kmer_contigs.txt"
        awk -F'\t' 'NR>1{print $1"\t"$3}' "$depth" > "$out/.depth_lookup.tsv"
        { printf 'contigName\t%s_depth\n' "$sample"
          awk -F'\t' 'NR==FNR{d[$1]=$2;next}{print $1"\t"(($1 in d)?d[$1]:0)}' "$out/.depth_lookup.tsv" "$out/.kmer_contigs.txt"; } > "$cov"
        rm -f "$out/.kmer_contigs.txt" "$out/.depth_lookup.tsv"
    fi

    # MetaBinner's KMeans (n_jobs=-1) leaves contigs unbinned at very high thread counts,
    # which then crashes split_hhbins.py (KeyError). Cap MetaBinner threads (other binners
    # keep $THREADS). 24 matches the original working run.
    local mb_t=$(( THREADS > 24 ? 24 : THREADS ))
    run_metabinner.sh \
        -a "$(readlink -f "$ref")" \
        -d "$(readlink -f "$cov")" \
        -k "$(readlink -f "$kmer")" \
        -o "$(readlink -f "$out")" \
        -p "$mb_dir" \
        -t "$mb_t" \
        -s small \
        > "$out/run.log" 2>&1 || true
}

run_metadecoder() {
    local sample=$1
    local ref="$REF_BASE/$sample/chromosomal.fasta"
    local bam="$BAM_DIR/${sample}.bam"
    local out="$BIN_BASE/metadecoder/$sample"
    [ -d "$out/bins" ] && [ "$(ls -A "$out/bins" 2>/dev/null)" ] && { echo "    metadecoder $sample done"; return 0; }
    mkdir -p "$out/bins"
    ( cd "$out" && \
        metadecoder coverage --threads "$THREADS" -b "$bam" -o coverage.tsv && \
        metadecoder seed --threads "$THREADS" -f "$ref" -o seed.txt && \
        metadecoder cluster -f "$ref" -c coverage.tsv -s seed.txt -o bins/bin \
    ) > "$out/run.log" 2>&1 || true
}

run_semibin_long() {
    local sample=$1
    local ref="$REF_BASE/$sample/chromosomal.fasta"
    local bam="$BAM_DIR/${sample}.bam"
    local out="$BIN_BASE/semibin/$sample"
    [ -d "$out/output_bins" ] && [ "$(ls -A "$out/output_bins" 2>/dev/null)" ] && { echo "    semibin $sample done"; return 0; }
    mkdir -p "$out"
    SemiBin2 single_easy_bin \
        --sequencing-type long_read \
        --self-supervised \
        -i "$ref" -b "$bam" -o "$out" \
        --threads "$THREADS" \
        > "$out/run.log" 2>&1 || true
}

run_semibin_short() {
    local sample=$1
    local ref="$REF_BASE/$sample/chromosomal.fasta"
    local bam="$BAM_DIR/${sample}.bam"
    local out="$BIN_BASE/semibin/$sample"
    [ -d "$out/output_bins" ] && [ "$(ls -A "$out/output_bins" 2>/dev/null)" ] && { echo "    semibin $sample done"; return 0; }
    mkdir -p "$out"
    SemiBin2 single_easy_bin \
        --sequencing-type short_read \
        --self-supervised \
        -i "$ref" -b "$bam" -o "$out" \
        --threads "$THREADS" \
        > "$out/run.log" 2>&1 || true
}

run_metabat2() {
    local sample=$1
    local ref="$REF_BASE/$sample/chromosomal.fasta"
    local depth="$DEPTH_DIR/${sample}_depth.txt"
    local out="$BIN_BASE/metabat2/$sample"
    [ -d "$out" ] && ls "$out"/bin.*.fa 2>/dev/null | grep -q . && { echo "    metabat2 $sample done"; return 0; }
    mkdir -p "$out"
    metabat2 -i "$ref" -a "$depth" -o "$out/bin" -t "$THREADS" \
        > "$out/run.log" 2>&1 || true
}

echo "[$(date '+%F %T')] binner dispatch (READ_TYPE=$READ_TYPE)"
SAMPLES=( $(list_samples) )

if [ "$READ_TYPE" = "hifi" ]; then
    for s in "${SAMPLES[@]}"; do
        echo "  $s — MetaBinner";   activate_env "$ENV_METABINNER";  run_metabinner   "$s"
        echo "  $s — MetaDecoder";  activate_env "$ENV_METADECODER"; run_metadecoder  "$s"
        echo "  $s — SemiBin2 (long)"; activate_env "$ENV_SEMIBIN";  run_semibin_long "$s"
    done
elif [ "$READ_TYPE" = "illumina" ]; then
    for s in "${SAMPLES[@]}"; do
        echo "  $s — MetaBAT2";     activate_env "$ENV_METABAT2";    run_metabat2     "$s"
        echo "  $s — SemiBin2 (short)"; activate_env "$ENV_SEMIBIN"; run_semibin_short "$s"
    done
else
    echo "ERROR: READ_TYPE must be hifi or illumina" >&2; exit 2
fi

echo "[$(date '+%F %T')] DONE"
