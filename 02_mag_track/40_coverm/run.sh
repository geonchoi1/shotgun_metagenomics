#!/bin/bash
# === 40 CoverM genome — relative abundance / trimmed_mean / covered_fraction / TPM ===
# Input:  $PROJECT/07_mag/08_drep_species/dereplicated_genomes/*.fa  (genome dir)
#         $PROJECT/01_qc/02_dehuman/*_clean.fastq.gz                  (HiFi reads)
# Output: $PROJECT/mag/40_coverm/<sample>.tsv
#         $PROJECT/mag/40_coverm/coverm_tpm_matrix.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

GENOME_DIR=$PROJECT/07_mag/08_drep_species/dereplicated_genomes
READS_DIR=$PROJECT/01_qc/02_dehuman
OUT=$PROJECT/mag/40_coverm
mkdir -p "$OUT"

[ -d "$GENOME_DIR" ] || { echo "ERROR: $GENOME_DIR missing" >&2; exit 1; }
[ -d "$READS_DIR" ]  || { echo "ERROR: $READS_DIR missing"  >&2; exit 1; }

activate_env "$ENV_COVERM"

# HiFi via minimap2-hifi; toggle MAPPER if Illumina (set MAPPER=bwa-mem)
MAPPER=${MAPPER:-minimap2-hifi}

echo "[$(date '+%F %T')] coverm genome (mapper=$MAPPER)"
for fq in "$READS_DIR"/*_clean.fastq.gz; do
    [ -f "$fq" ] || continue
    sample=$(basename "$fq" _clean.fastq.gz)
    out="$OUT/${sample}.tsv"
    if [ -s "$out" ]; then
        echo "  $sample already done"
        continue
    fi
    echo "  $sample"
    coverm genome \
        --single "$fq" \
        --genome-fasta-directory "$GENOME_DIR" \
        --genome-fasta-extension fa \
        -p "$MAPPER" \
        -m relative_abundance trimmed_mean covered_fraction tpm \
        --min-read-aligned-percent 75 \
        --min-read-percent-identity 95 \
        --min-covered-fraction 0 \
        -t "$THREADS" \
        -o "$out"
done

# Combined TPM matrix
echo "[$(date '+%F %T')] build TPM matrix"
python3 - <<PYEOF
import os, glob
OUT = "$OUT"
files = sorted(glob.glob(os.path.join(OUT, "*.tsv")))
files = [f for f in files if not f.endswith("coverm_tpm_matrix.tsv")]
data = {}; samples = []
for f in files:
    s = os.path.basename(f)[:-4]
    samples.append(s)
    with open(f) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        # find TPM column (case-insensitive)
        tpm_idx = next((i for i,h in enumerate(header) if h.strip().lower().endswith("tpm")), None)
        if tpm_idx is None:
            print(f"  WARN no TPM col in {f}"); continue
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            g = parts[0]
            try: v = float(parts[tpm_idx])
            except: v = 0.0
            data.setdefault(g, {})[s] = v
with open(os.path.join(OUT, "coverm_tpm_matrix.tsv"), "w") as fh:
    fh.write("genome\t" + "\t".join(samples) + "\n")
    for g, row in sorted(data.items()):
        fh.write(g + "\t" + "\t".join(f"{row.get(s,0.0):.4f}" for s in samples) + "\n")
print(f"  TPM matrix: {len(data)} genomes x {len(samples)} samples")
PYEOF

echo "[$(date '+%F %T')] DONE"
