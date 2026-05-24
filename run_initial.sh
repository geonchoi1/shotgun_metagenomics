#!/bin/bash
# === run_initial.sh ===
# 00_shared: reads → assembly → genomad → topology split → chromosomal → MAG production → split UB
set -euo pipefail
REPO="$(cd "$(dirname "$0")" && pwd)"
source "$REPO/config.sh"
PROJECT=${1:?ERROR: usage: bash run_initial.sh /path/to/project}; export PROJECT
S="$REPO/00_shared"
: ${READ_TYPE:?ERROR: export READ_TYPE=illumina or hifi}

echo "[$(date '+%F %T')] === 00_shared START ==="

# 01 Read QC
if [ "$READ_TYPE" = "illumina" ]; then
    bash "$S/01_read_qc/illumina_fastp.sh"
    bash "$S/01_read_qc/illumina_dehuman_bowtie2.sh"
elif [ "$READ_TYPE" = "hifi" ]; then
    bash "$S/01_read_qc/hifi_dehuman_minimap2.sh"
fi

# 02 Assembly
if [ "$READ_TYPE" = "illumina" ]; then
    bash "$S/02_assembly/illumina_metaspades.sh"
else
    bash "$S/02_assembly/hifi_metaflye.sh"
fi

# 03 Virus identification (geNomad default)
bash "$S/03_genomad_virus/run.sh"

# 04 Plasmid identification (-s 4.8 --relaxed + 5-filter)
bash "$S/04_genomad_plasmid/01_run_genomad_relaxed.sh"
bash "$S/04_genomad_plasmid/02_filter_F1234.sh"
bash "$S/04_genomad_plasmid/03_filter_F5_rrna.sh"

# 05 Topology split (uses user-provided circ_frag_map.tsv)
bash "$S/05_topology_split/split_from_user_list.sh"

# 06 Chromosomal extract
bash "$S/06_chromosomal_extract/run.sh"

# 07 MAG production
bash "$S/07_mag_production/01_mapping/run.sh"
bash "$S/07_mag_production/02_depth/run.sh"
bash "$S/07_mag_production/03_binner/run.sh"
bash "$S/07_mag_production/04_dastool/run.sh"
bash "$S/07_mag_production/05_checkm2/run.sh"
bash "$S/07_mag_production/06_rrna_barrnap/run.sh"
bash "$S/07_mag_production/07_trna_trnascan/run.sh"
bash "$S/07_mag_production/08_drep_cross_sample/run.sh"
bash "$S/07_mag_production/08b_drep_strain/run.sh"
bash "$S/07_mag_production/09_mimag_classify/run.sh"
bash "$S/07_mag_production/10_gtdbtk/run.sh"

# 08 Split chromosomal contigs → binned + unbinned
bash "$S/08_split_binned_unbinned/run.sh"

echo "[$(date '+%F %T')] === 00_shared DONE ==="
