#!/bin/bash
# === run_mag.sh — 02_mag_track ===
set -euo pipefail
REPO="$(cd "$(dirname "$0")" && pwd)"
source "$REPO/config.sh"
PROJECT=${1:?ERROR: usage: bash run_mag.sh /path/to/project}; export PROJECT
T="$REPO/02_mag_track"

echo "[$(date '+%F %T')] === 02_mag_track START ==="

bash "$T/02_bakta/run.sh"
bash "$T/03_master_orf/run.sh"

for step in 04_pfam 05_ncbifam 06_kofamscan 07_eggnog 08_amrfinder 09_bacmet 10_vfdb 11_tadb \
            12_dbapis 13_dbcan 14_macrel 15_defensefinder 16_antismash 17_bigscape 18_dbscan_swa \
            19_isescan 20_integronfinder 21_iceberg3; do
    [ -f "$T/$step/run.sh" ] && bash "$T/$step/run.sh"
done

bash "$T/22_cctyper/run.sh"
bash "$T/30_metabolic_g/run.sh"
bash "$T/40_coverm/run.sh"
bash "$T/41_gene_abundance/run.sh"

echo "[$(date '+%F %T')] === 02_mag_track DONE ==="
