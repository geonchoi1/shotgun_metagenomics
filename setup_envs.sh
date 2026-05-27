#!/bin/bash
# ==============================================================================
# Create the conda envs used by this pipeline (idempotent — skips existing)
# ==============================================================================
# Usage:
#   bash setup_envs.sh                # create all
#   bash setup_envs.sh bakta mob_suite # create only listed envs
# ==============================================================================
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || \
source "$HOME/anaconda3/etc/profile.d/conda.sh"

have_env() { conda env list | awk '{print $1}' | grep -qx "$1"; }

create_yml() {
    if have_env "$1"; then echo "[skip] $1"; return; fi
    echo "[create] $1 from environment.yml"
    conda env create -f "$(dirname "$0")/environment.yml" -n "$1"
}

create_one() {
    local name=$1; shift
    if have_env "$name"; then echo "[skip] $name"; return; fi
    echo "[create] $name <- $*"
    conda create -y -n "$name" -c conda-forge -c bioconda "$@"
}

# Restrict to specific envs if listed on CLI
WANT=("$@")
should() {
    [ ${#WANT[@]} -eq 0 ] && return 0
    for w in "${WANT[@]}"; do [ "$w" = "$1" ] && return 0; done
    return 1
}

# --- Analysis + visualization env -------------------------------------------
should shotgun-analysis && create_yml shotgun-analysis

# --- Read QC / assembly -----------------------------------------------------
should fastp           && create_one fastp           fastp
should bowtie2         && create_one bowtie2         bowtie2 samtools
should coverm          && create_one coverm          coverm minimap2 samtools
should spades          && create_one spades          spades
should metaflye        && create_one metaflye        flye

# --- Annotation core --------------------------------------------------------
should bakta           && create_one bakta           bakta
should infernal        && create_one infernal        infernal
should kofamscan       && create_one kofamscan       kofamscan
should eggnog-mapper   && create_one eggnog-mapper   eggnog-mapper
should amrfinderplus   && create_one amrfinderplus   ncbi-amrfinderplus
should dbcan           && create_one dbcan           dbcan
should macrel          && create_one macrel          macrel
should defense-finder  && create_one defense-finder  defense-finder
should antismash       && create_one antismash       antismash
should bigscape        && create_one bigscape        bigscape
should phage           && create_one phage           barrnap pyrodigal
should isescan         && create_one isescan         isescan
should integronfinder  && create_one integronfinder  integron_finder

# --- MAG production ---------------------------------------------------------
should metabinner      && create_one metabinner      metabinner
should metadecoder     && create_one metadecoder     pip
should semibin         && create_one semibin         semibin
should metabat2        && create_one metabat2        metabat2
should das_tool        && create_one das_tool        das_tool diamond
should checkm2         && create_one checkm2         checkm2
should drep            && create_one drep            drep
should gtdbtk          && create_one gtdbtk          gtdbtk

# --- Plasmid typing (defaults) ----------------------------------------------
should mob_suite       && create_one mob_suite       mob_suite
should plasmidfinder   && create_one plasmidfinder   plasmidfinder
should macsyfinder     && create_one macsyfinder     macsyfinder

# --- Plasmid clustering / host / PLSDB --------------------------------------
should copla           && create_one copla           pip
should mummer4         && create_one mummer4         mummer4
should simka           && create_one simka           simka
should iphop           && create_one iphop           iphop
should cctyper         && create_one cctyper         cctyper

# --- MAG-specific metabolism ------------------------------------------------
should metabolic       && create_one metabolic       pip perl gawk hmmer
should dram            && create_one dram            dram

# --- UB taxonomy ------------------------------------------------------------
should mmseqs2         && create_one mmseqs2         mmseqs2

# --- Community-level read taxonomy ------------------------------------------
should kraken2         && create_one kraken2         kraken2 bracken

# --- Plasmid clustering (Camargo pipeline) ----------------------------------
should snakemake-camargo && create_one snakemake-camargo snakemake blast python-igraph leidenalg numpy seqkit networkx

# --- Cross-track (Mobile ARG R) ---------------------------------------------
should r-base          && create_one r-base          r-base r-essentials

echo ""
echo "DONE. Envs created/verified."
echo "Note: some tools (oriTfinder2, Ares-Arroyo, MetaDecoder) need extra install — see tools/install_*.sh + per-env post-install steps."
