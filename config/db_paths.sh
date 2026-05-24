#!/bin/bash
# === Database paths for shotgun_metagenomics pipeline ===
# All scripts `source` this file to resolve DB locations.
# Override with environment variables before running each script:
#   PFAM_HMM=/custom/path/Pfam-A.hmm bash 05_pfam/run.sh

# Root of all reference DBs
export DB_ROOT=${DB_ROOT:-/mnt/nas/DB/geon}

# Read QC / human removal
export GRCH38=${GRCH38:-$DB_ROOT/GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz}
export GRCH38_BT2_IDX=${GRCH38_BT2_IDX:-$DB_ROOT/GRCh38.p14/bt2_index/GRCh38}

# Assembly / classification
export GENOMAD_DB=${GENOMAD_DB:-$DB_ROOT/genomad_db}
export RFAM_DB=${RFAM_DB:-$DB_ROOT/rfam}             # for 5-filter F5 rRNA cmscan
export RFAM_CM=${RFAM_CM:-$RFAM_DB/rRNA_5.cm}        # combined 5 rRNA CMs (RF00001/00177/01959/02540/02541)

# Bakta
export BAKTA_DB=${BAKTA_DB:-$DB_ROOT/bakta_db/db}

# Functional HMM DBs
export PFAM_HMM=${PFAM_HMM:-$DB_ROOT/pfam/Pfam-A.hmm}
export PFAM_HALLMARK=${PFAM_HALLMARK:-$DB_ROOT/pfam_hallmark}    # PF03090/PF03432/PF07042/PF12696/TIGR00929
export NCBIFAM_HMM=${NCBIFAM_HMM:-$DB_ROOT/ncbifam/hmm_PGAP.LIB}
export KOFAM_PROFILES=${KOFAM_PROFILES:-$DB_ROOT/kofam_db/profiles}
export KOFAM_KO_LIST=${KOFAM_KO_LIST:-$DB_ROOT/kofam_db/ko_list}
export EGGNOG_DATA=${EGGNOG_DATA:-$DB_ROOT/eggnog_db}
export DBAPIS_HMM=${DBAPIS_HMM:-$DB_ROOT/dbAPIS/dbAPIS.hmm}

# Specific annotation DBs (DIAMOND-based)
export AMRFINDER_DB=${AMRFINDER_DB:-$DB_ROOT/amrfinder_db/latest}
export BACMET_DMND=${BACMET_DMND:-$DB_ROOT/bacmet_db/BacMet.dmnd}
export VFDB_DMND=${VFDB_DMND:-$DB_ROOT/vfdb_db/VFDB.dmnd}
export TADB_DMND=${TADB_DMND:-$DB_ROOT/tadb_db/TADB.dmnd}
export DBCAN_DB=${DBCAN_DB:-$DB_ROOT/dbcan_db}
export DEFENSEFINDER_MODELS=${DEFENSEFINDER_MODELS:-$DB_ROOT/defense-finder-models-v3.1}

# MGE detection
export ISESCAN_DB=${ISESCAN_DB:-$DB_ROOT/isescan_db}
export ICEBERG3_DB=${ICEBERG3_DB:-$DB_ROOT/iceberg3/ICE_combined}

# Plasmid typing
export MOB_SUITE_DB=${MOB_SUITE_DB:-$DB_ROOT/mob_suite}
export MOBSCAN_HMM=${MOBSCAN_HMM:-$DB_ROOT/mobscan_db/MOBfamDB}
export PLASMIDFINDER_DB=${PLASMIDFINDER_DB:-$DB_ROOT/plasmidfinder_db2}
export CONJSCAN_MODELS_DIR=${CONJSCAN_MODELS_DIR:-$DB_ROOT/conjscan_models}
export ORITFINDER2_DIR=${ORITFINDER2_DIR:-$HOME/tools/oriTfinder2_linux}
export ARES_ARROYO_DB=${ARES_ARROYO_DB:-$HOME/tools/ares_arroyo_oriT/ares_arroyo_91_oriT.fna}

# MAG production
export CHECKM2_DB=${CHECKM2_DB:-$DB_ROOT/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd}
export GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH:-$DB_ROOT/gtdbtk_db/release232}
export GUNC_DB=${GUNC_DB:-$DB_ROOT/gunc_db/gunc_db_progenomes2.1.dmnd}

# Host prediction
export IPHOP_DB=${IPHOP_DB:-$DB_ROOT/iphop_db/Jun_2025_pub_rw}
export IPHOP_SPACER_DB=${IPHOP_SPACER_DB:-$DB_ROOT/iphop_db/Jun_2025_pub_rw/db/All_CRISPR_spacers_nr_clean}

# PLSDB
export PLSDB_DIR=${PLSDB_DIR:-$DB_ROOT/plsdb}
export PLSDB_FASTA=${PLSDB_FASTA:-$PLSDB_DIR/sequences.fasta}
export PLSDB_MASH_SKETCH=${PLSDB_MASH_SKETCH:-$PLSDB_DIR/metadata/plsdb_sketch.msh}

# MMseqs2 GTDB (UB contig taxonomy)
export MMSEQS_GTDB_DB=${MMSEQS_GTDB_DB:-/mnt/nas/DB/nj/mmseqs2_db/gtdbAA_DB}

# COPLA (Track 2 plasmid PTU)
export COPLA_DIR=${COPLA_DIR:-$DB_ROOT/copla_install/COPLA}

# Optional METABOLIC / DRAM (MAG-only)
export METABOLIC_DB=${METABOLIC_DB:-$DB_ROOT/metabolic_db}
export DRAM_CONFIG=${DRAM_CONFIG:-$HOME/.dram_config.json}
