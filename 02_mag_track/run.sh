#!/bin/bash
# === 02_mag_track/run.sh ===
# MAG annotation + MAG-specific tools (cctyper, METABOLIC-G, CoverM, gene-level abundance)
set -e
BASE=$(cd "$(dirname "$0")" && pwd)

# 02 Bakta (circ --complete + frag default for complete-circular MAGs)
bash $BASE/02_bakta/run.sh

# 03 Master ORF (circ/frag/all + per_mag)
bash $BASE/03_master_orf/run.sh

# 04-21 Annotation (same set as plasmid: Pfam ~ ICEberg3)
for step in 04_pfam 05_ncbifam 06_kofamscan 07_eggnog 08_amrfinder 09_bacmet 10_vfdb 11_tadb \
            12_dbapis 13_dbcan 14_macrel 15_defensefinder 16_antismash 17_bigscape 18_dbscan_swa \
            19_isescan 20_integronfinder 21_iceberg3; do
    [ -f $BASE/$step/run.sh ] && bash $BASE/$step/run.sh || echo "  skip $step (no run.sh)"
done

# 22 MAG-specific: cctyper (CRISPR-Cas + spacer)
bash $BASE/22_cctyper/run.sh

# 30 MAG-specific: METABOLIC-G (KEGG module + biogeochemical pathway)
bash $BASE/30_metabolic_g/run.sh

# 40 Abundance: CoverM at genome level (per MAG per sample TPM)
bash $BASE/40_coverm/run.sh

# 41 Gene-level abundance pipeline (minimap2 + featureCounts + TPM + pathway aggregation)
bash $BASE/41_gene_abundance/run.sh

echo "[$(date '+%F %T')] 02_mag_track DONE"
