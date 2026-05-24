#!/bin/bash
# === 01_plasmid_track/run.sh ===
# Plasmid annotation + mobility typing + clustering + host prediction + PLSDB
set -e
BASE=$(cd "$(dirname "$0")" && pwd)

# 02 dRep (RBH 100/100 ANI + AF)
bash $BASE/02_drep/run.sh

# 03 Bakta (circ --complete + frag default)
bash $BASE/03_bakta/run.sh

# 04 Master ORF (circ/frag/all)
bash $BASE/04_master_orf/run.sh

# 05-22 Annotation (all run on master.faa or master.fna)
for step in 05_pfam 06_ncbifam 07_kofamscan 08_eggnog 09_amrfinder 10_bacmet 11_vfdb 12_tadb \
            13_dbapis 14_dbcan 15_macrel 16_defensefinder 17_antismash 18_bigscape 19_dbscan_swa \
            20_isescan 21_integronfinder 22_iceberg3; do
    [ -f $BASE/$step/run.sh ] && bash $BASE/$step/run.sh || echo "  skip $step (no run.sh)"
done

# 23-26 Mobility typing (4 defaults)
bash $BASE/23_mob_typer/run.sh
bash $BASE/24_mobscan/run.sh
python3 $BASE/25_oritfinder2/01_split_gbff.py
bash $BASE/25_oritfinder2/02_parallel.sh
bash $BASE/25_oritfinder2/03_finalize.sh
bash $BASE/26_conjscan/run.sh

# 27 Alternative mobility tools (optional supplementary)
# Run individually as needed; not part of default pipeline:
#   bash 27_mobility_typing_alt/rep/*.sh
#   python 27_mobility_typing_alt/mob/*.py; bash 27_mobility_typing_alt/mob/*.sh
#   ...

# 28 5-tier classification (consumes 23-26 outputs)
python3 $BASE/28_5tier_classification/classify_5tier.py

# 30 Clustering (Track 1 + 2 + 3 + validation)
bash $BASE/30_clustering/run.sh

# 31 AcCNET (internal + external NMI)
bash $BASE/31_accnet/run.sh

# 32 Functional comparison (PCA + GSEA + Fisher × richness + TPM-weighted)
bash $BASE/32_functional_comparison/run.sh

# 40 Quantification (minimap2 + CoverM TPM/mean/RPKM)
bash $BASE/40_quantification/run.sh

# 50 Host prediction (Track A + C + D + union)
bash $BASE/50_host_prediction/run.sh

# 60 PLSDB lookup + ecosystem positioning
bash $BASE/60_plsdb_lookup/run.sh

echo "[$(date '+%F %T')] 01_plasmid_track DONE"
