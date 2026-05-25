#!/bin/bash
# === run_plasmid.sh — 01_plasmid_track ===
set -euo pipefail
REPO="$(cd "$(dirname "$0")" && pwd)"
source "$REPO/config.sh"
PROJECT=${1:?ERROR: usage: bash run_plasmid.sh /path/to/project}; export PROJECT
T="$REPO/01_plasmid_track"

echo "[$(date '+%F %T')] === 01_plasmid_track START ==="

bash "$T/02_drep/run.sh"
bash "$T/03_bakta/run.sh"
bash "$T/04_master_orf/run.sh"

for step in 05_pfam 06_ncbifam 07_kofamscan 08_eggnog 09_amrfinder 10_bacmet 11_vfdb 12_tadb \
            14_dbcan 15_macrel 16_defensefinder 17_antismash 18_bigscape 19_dbscan_swa \
            20_isescan 21_integronfinder 22_iceberg3; do
    [ -f "$T/$step/run.sh" ] && bash "$T/$step/run.sh"
done

bash "$T/23_mob_typer/run.sh"
bash "$T/24_mobscan/run.sh"
python3 "$T/25_oritfinder2/01_split_gbff.py"
bash "$T/25_oritfinder2/02_parallel.sh"
bash "$T/25_oritfinder2/03_finalize.sh"
bash "$T/26_conjscan/run.sh"
# 27_mobility_typing_alt: optional — run individually as needed
python3 "$T/28_5tier_classification/classify_5tier.py"

bash "$T/30_clustering/run.sh"
bash "$T/31_accnet/run.sh"
bash "$T/32_functional_comparison/run.sh"
bash "$T/40_quantification/run.sh"
bash "$T/50_host_prediction/run.sh"
bash "$T/60_plsdb_lookup/run.sh"

echo "[$(date '+%F %T')] === 01_plasmid_track DONE ==="
