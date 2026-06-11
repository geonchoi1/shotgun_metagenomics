#!/bin/bash
# === 09 Kraken2 + Bracken community-level taxonomy ===
# Read-level classification (independent of MAG/bin assembly).
# Answers "what's in the community" — complements MAG GTDB-Tk
# which answers "what was recovered as a genome".
#
# Input:  $PROJECT/00_shared/01_read_qc/dehuman/<SAMPLE>_clean.fastq.gz       (HiFi, single)
#         $PROJECT/00_shared/01_read_qc/dehuman/<SAMPLE>_clean_R{1,2}.fastq.gz (Illumina, paired)
# Output: $PROJECT/00_shared/09_kraken2_community/<SAMPLE>/
#           kraken2.report, kraken2.output, bracken.S.report, bracken.G.report
#         $PROJECT/00_shared/09_kraken2_community/merged_bracken_species.tsv
#         $PROJECT/00_shared/09_kraken2_community/merged_bracken_genus.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}
: ${READ_TYPE:?ERROR: export READ_TYPE=hifi or illumina}

IN_DIR=$PROJECT/00_shared/01_read_qc/dehuman
OUT_BASE=$PROJECT/00_shared/09_kraken2_community
mkdir -p "$OUT_BASE"

activate_env "$ENV_KRAKEN2"

# Discover samples
if [ "$READ_TYPE" = "hifi" ]; then
    SAMPLES=$(ls "$IN_DIR"/*_clean.fastq.gz 2>/dev/null | xargs -n1 basename | sed 's/_clean\.fastq\.gz$//' | sort -u)
else
    SAMPLES=$(ls "$IN_DIR"/*_clean_R1.fastq.gz 2>/dev/null | xargs -n1 basename | sed 's/_clean_R1\.fastq\.gz$//' | sort -u)
fi
[ -n "$SAMPLES" ] || { echo "ERROR: no clean reads in $IN_DIR" >&2; exit 1; }

# === Kraken2 + Bracken per sample ===
for SAMPLE in $SAMPLES; do
    OUT=$OUT_BASE/$SAMPLE
    mkdir -p "$OUT"
    if [ -s "$OUT/bracken.S.report" ]; then
        echo "[$(date '+%F %T')] $SAMPLE already done — skip"
        continue
    fi

    echo "[$(date '+%F %T')] Kraken2 — $SAMPLE"
    if [ "$READ_TYPE" = "hifi" ]; then
        kraken2 \
            --db "$KRAKEN2_DB" \
            --threads "$THREADS" \
            --gzip-compressed \
            --report "$OUT/kraken2.report" \
            --output "$OUT/kraken2.output" \
            "$IN_DIR/${SAMPLE}_clean.fastq.gz"
    else
        kraken2 \
            --db "$KRAKEN2_DB" \
            --threads "$THREADS" \
            --paired \
            --gzip-compressed \
            --report "$OUT/kraken2.report" \
            --output "$OUT/kraken2.output" \
            "$IN_DIR/${SAMPLE}_clean_R1.fastq.gz" \
            "$IN_DIR/${SAMPLE}_clean_R2.fastq.gz"
    fi

    echo "[$(date '+%F %T')] Bracken (species + genus) — $SAMPLE"
    bracken -d "$KRAKEN2_DB" -i "$OUT/kraken2.report" \
            -o "$OUT/bracken.S.tsv" -w "$OUT/bracken.S.report" \
            -r "$BRACKEN_READ_LEN" -l S -t 10
    bracken -d "$KRAKEN2_DB" -i "$OUT/kraken2.report" \
            -o "$OUT/bracken.G.tsv" -w "$OUT/bracken.G.report" \
            -r "$BRACKEN_READ_LEN" -l G -t 10
done

# === Merge per-sample Bracken tables into sample × taxon matrix ===
echo "[$(date '+%F %T')] Merging Bracken outputs"
for L in S G; do
    OUT_TSV=$OUT_BASE/merged_bracken_${L,,}.tsv   # species → s, genus → g
    [ "$L" = "S" ] && OUT_TSV=$OUT_BASE/merged_bracken_species.tsv
    [ "$L" = "G" ] && OUT_TSV=$OUT_BASE/merged_bracken_genus.tsv

    combine_bracken_outputs.py \
        --files $(ls $OUT_BASE/*/bracken.${L}.tsv | sort) \
        --output "$OUT_TSV"
done

echo "[$(date '+%F %T')] DONE — see $OUT_BASE/"
echo "  merged_bracken_species.tsv, merged_bracken_genus.tsv"
echo "  per-sample: <SAMPLE>/{kraken2.report, kraken2.output, bracken.{S,G}.report}"
