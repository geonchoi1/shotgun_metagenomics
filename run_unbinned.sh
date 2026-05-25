#!/bin/bash
# === run_unbinned.sh — 03_unbinned_track ===
set -euo pipefail
REPO="$(cd "$(dirname "$0")" && pwd)"
source "$REPO/config.sh"
PROJECT=${1:?ERROR: usage: bash run_unbinned.sh /path/to/project}; export PROJECT
T="$REPO/03_unbinned_track"

echo "[$(date '+%F %T')] === 03_unbinned_track START ==="

bash "$T/02_bakta/run.sh"
bash "$T/03_master_orf/run.sh"

for step in 04_pfam 05_ncbifam 06_kofamscan 07_eggnog 08_amrfinder 13_dbcan 18_dbscan_swa \
            19_isescan 20_integronfinder 21_iceberg3; do
    [ -f "$T/$step/run.sh" ] && bash "$T/$step/run.sh"
done

bash "$T/22_cctyper/run.sh"

echo "[$(date '+%F %T')] === 03_unbinned_track DONE ==="
