#!/bin/bash
# === 03_unbinned_track/run.sh ===
# Unbinned chromosomal contig annotation (gene-level only, no genome-level analysis)
set -e
BASE=$(cd "$(dirname "$0")" && pwd)

# 02 Bakta
bash $BASE/02_bakta/run.sh

# 03 Master ORF
bash $BASE/03_master_orf/run.sh

# 04-21 Annotation (same set as MAG)
for step in 04_pfam 05_ncbifam 06_kofamscan 07_eggnog 08_amrfinder 09_bacmet 10_vfdb 11_tadb \
            12_dbapis 13_dbcan 14_macrel 15_defensefinder 16_antismash 17_bigscape 18_dbscan_swa \
            19_isescan 20_integronfinder 21_iceberg3; do
    [ -f $BASE/$step/run.sh ] && bash $BASE/$step/run.sh || echo "  skip $step (no run.sh)"
done

# 22 UB-specific: MMseqs2 LCA contig-level taxonomy
bash $BASE/22_mmseqs2_lca_taxonomy/run.sh

echo "[$(date '+%F %T')] 03_unbinned_track DONE"
