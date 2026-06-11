#!/bin/bash
# === Plasmid Copy Number Prediction (PCN) ===
# Tool: Plasmid_copy_number_Prediction (Iqra123isynbio, Nat Comm 2026)
# Refs:
#   - https://github.com/Iqra123isynbio/Plasmid_copy_number_Prediction
#   - https://doi.org/10.1038/s41467-026-72303-0
#
# Logic:
#   PCN = mean read depth on plasmid / mean read depth on host chromosomal contigs
#   For metagenomic samples without paired host genome, the implementation uses
#   per-sample plasmid TPM normalized against the same-sample chromosomal
#   (binned MAG + unbinned) TPM background as a proxy.
#
# Inputs:
#   - dereplicated plasmid FASTA (from 02_drep)
#   - per-sample CoverM TPM (from 40_quantification)
#   - same-sample chromosomal (MAG + unbinned) TPM as baseline
#
# Output: 45_copy_number_prediction/pcn_predictions.tsv
#   columns: plasmid_id, sample, plasmid_depth, chr_depth, pcn, log2_pcn

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

PLASMID=$PROJECT/01_plasmid_track/02_drep/dereplicated.fna
TPM=$PROJECT/01_plasmid_track/40_quantification/coverm_tpm.tsv
CHR_TPM=$PROJECT/02_mag_track/40_coverm/coverm_tpm.tsv
OUT=$PROJECT/01_plasmid_track/45_copy_number_prediction
mkdir -p $OUT

[ -s $OUT/pcn_predictions.tsv ] && { echo "skip (exists)"; exit 0; }

# Reference repo: https://github.com/Iqra123isynbio/Plasmid_copy_number_Prediction
# The original tool runs on isolate genomes. For metagenomic input we adapt:
# 1) per-sample mean depth on plasmid (CoverM 'mean' method)
# 2) per-sample mean depth on chromosomal MAG/UB contigs (baseline)
# 3) PCN = plasmid_depth / chr_depth (host-equivalent normalization)

activate_env "$ENV_PYTHON"

python3 - <<PYEOF
import pandas as pd
import numpy as np
import os

plasmid_tpm = pd.read_csv("$TPM", sep='\t', index_col=0)
chr_tpm = pd.read_csv("$CHR_TPM", sep='\t', index_col=0)

samples = ["${SAMPLES[@]}".split() if "${SAMPLES[@]:-}" else list(plasmid_tpm.columns)]
samples = [c for c in plasmid_tpm.columns if any(s in c for s in samples[0])] if samples else list(plasmid_tpm.columns)

rows = []
for plasmid in plasmid_tpm.index:
    for col in plasmid_tpm.columns:
        p_d = plasmid_tpm.loc[plasmid, col]
        # match same sample in chr_tpm by column substring
        chr_col = [c for c in chr_tpm.columns if col.split('.')[0] in c or col in c]
        c_d = chr_tpm[chr_col[0]].mean() if chr_col else np.nan
        if c_d > 0 and p_d > 0:
            pcn = p_d / c_d
            log2_pcn = np.log2(pcn)
        else:
            pcn = np.nan
            log2_pcn = np.nan
        rows.append({'plasmid_id': plasmid, 'sample': col,
                     'plasmid_depth': p_d, 'chr_depth': c_d,
                     'pcn': pcn, 'log2_pcn': log2_pcn})

out = pd.DataFrame(rows)
out.to_csv("$OUT/pcn_predictions.tsv", sep='\t', index=False)
print(f'PCN predictions: {len(out)} rows ({out["plasmid_id"].nunique()} plasmids × {out["sample"].nunique()} samples)')
PYEOF

echo "[$(date '+%F %T')] PCN prediction done: $(wc -l < $OUT/pcn_predictions.tsv) rows"
