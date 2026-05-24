#!/bin/bash
# === 00_shared/run.sh ===
# Reads → QC → assembly → classify → bin (MAG production).
set -e
BASE=$(cd "$(dirname "$0")" && pwd)
: ${READ_TYPE:?ERROR: export READ_TYPE=illumina or hifi}

# 01 Read QC
if [ "$READ_TYPE" = "illumina" ]; then
    bash $BASE/01_read_qc/illumina_fastp.sh
    bash $BASE/01_read_qc/illumina_dehuman_bowtie2.sh
elif [ "$READ_TYPE" = "hifi" ]; then
    bash $BASE/01_read_qc/hifi_dehuman_minimap2.sh
fi

# 02 Assembly
if [ "$READ_TYPE" = "illumina" ]; then
    bash $BASE/02_assembly/illumina_metaspades.sh
else
    bash $BASE/02_assembly/hifi_metaflye.sh
fi

# 03 geNomad virus (default)
bash $BASE/03_genomad_virus/run.sh

# 04 geNomad plasmid (-s 4.8 --relaxed + 5-filter F1-F5)
bash $BASE/04_genomad_plasmid/01_run_genomad_relaxed.sh
bash $BASE/04_genomad_plasmid/02_filter_F1234.sh
bash $BASE/04_genomad_plasmid/03_filter_F5_rrna.sh

# 05 Topology split (uses user circ_frag_map.tsv)
bash $BASE/05_topology_split/split_from_user_list.sh

# 06 Chromosomal extract (assembly − plasmid − virus)
bash $BASE/06_chromosomal_extract/run.sh

# 07 MAG production
bash $BASE/07_mag_production/01_mapping/run.sh
bash $BASE/07_mag_production/02_depth/run.sh
bash $BASE/07_mag_production/03_binner/run.sh
bash $BASE/07_mag_production/04_dastool/run.sh
bash $BASE/07_mag_production/05_checkm2/run.sh
bash $BASE/07_mag_production/06_rrna_barrnap/run.sh
bash $BASE/07_mag_production/07_trna_trnascan/run.sh
bash $BASE/07_mag_production/08_drep_cross_sample/run.sh      # species (95% ANI, main)
bash $BASE/07_mag_production/08b_drep_strain/run.sh           # strain (99% ANI, supp)
bash $BASE/07_mag_production/09_mimag_classify/run.sh
bash $BASE/07_mag_production/10_gtdbtk/run.sh

# 08 Split chromosomal into binned vs unbinned
bash $BASE/08_split_binned_unbinned/run.sh

echo "[$(date '+%F %T')] 00_shared DONE"
