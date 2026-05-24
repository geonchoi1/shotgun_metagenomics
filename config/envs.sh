#!/bin/bash
# === Conda environment names ===
# Each script `source`s this file to standardize env names.
# If your conda env names differ, override here.

# Read QC + assembly
export ENV_FASTP=${ENV_FASTP:-fastp}
export ENV_BOWTIE2=${ENV_BOWTIE2:-bowtie2}
export ENV_MINIMAP2=${ENV_MINIMAP2:-coverm}             # minimap2 inside coverm env
export ENV_METASPADES=${ENV_METASPADES:-spades}
export ENV_METAFLYE=${ENV_METAFLYE:-metaflye}

# Annotation core
export ENV_BAKTA=${ENV_BAKTA:-bakta}
export ENV_HMMER=${ENV_HMMER:-base}                     # hmmsearch/hmmpress
export ENV_DIAMOND=${ENV_DIAMOND:-base}
export ENV_INFERNAL=${ENV_INFERNAL:-infernal}           # for 5-filter F5

# Functional tools
export ENV_KOFAMSCAN=${ENV_KOFAMSCAN:-kofamscan}
export ENV_EGGNOG=${ENV_EGGNOG:-eggnog-mapper}
export ENV_AMRFINDER=${ENV_AMRFINDER:-amrfinderplus}
export ENV_DBCAN=${ENV_DBCAN:-dbcan}
export ENV_MACREL=${ENV_MACREL:-macrel}
export ENV_DEFENSEFINDER=${ENV_DEFENSEFINDER:-defense-finder}
export ENV_ANTISMASH=${ENV_ANTISMASH:-antismash}
export ENV_BIGSCAPE=${ENV_BIGSCAPE:-bigscape}
export ENV_DBSCANSWA=${ENV_DBSCANSWA:-phage}

# MGE detection
export ENV_ISESCAN=${ENV_ISESCAN:-isescan}
export ENV_INTEGRONFINDER=${ENV_INTEGRONFINDER:-integronfinder}

# MAG production
export ENV_METABINNER=${ENV_METABINNER:-metabinner}
export ENV_METADECODER=${ENV_METADECODER:-metadecoder}
export ENV_SEMIBIN=${ENV_SEMIBIN:-semibin}
export ENV_METABAT2=${ENV_METABAT2:-metabat2}
export ENV_DASTOOL=${ENV_DASTOOL:-das_tool}
export ENV_CHECKM2=${ENV_CHECKM2:-checkm2}
export ENV_BARRNAP=${ENV_BARRNAP:-phage}                # barrnap inside phage env
export ENV_TRNASCAN=${ENV_TRNASCAN:-bakta}              # tRNAscan-SE inside bakta env
export ENV_DREP=${ENV_DREP:-drep}
export ENV_GTDBTK=${ENV_GTDBTK:-gtdbtk}
export ENV_COVERM=${ENV_COVERM:-coverm}

# Plasmid typing (4 defaults)
export ENV_MOB_SUITE=${ENV_MOB_SUITE:-mob_suite}        # mob_typer
export ENV_PLASMIDFINDER=${ENV_PLASMIDFINDER:-plasmidfinder}
export ENV_MACSYFINDER=${ENV_MACSYFINDER:-macsyfinder}  # CONJScan

# Plasmid clustering / host / PLSDB
export ENV_COPLA=${ENV_COPLA:-copla}
export ENV_MASH=${ENV_MASH:-base}
export ENV_NUCMER=${ENV_NUCMER:-mummer4}
export ENV_SIMKA=${ENV_SIMKA:-simka}
export ENV_IPHOP=${ENV_IPHOP:-iphop}
export ENV_CCTYPER=${ENV_CCTYPER:-cctyper}

# Cross-track (Mobile ARG)
export ENV_R=${ENV_R:-r-base}                           # CooccurrenceAffinity

# MAG metabolism
export ENV_METABOLIC=${ENV_METABOLIC:-metabolic}
export ENV_DRAM=${ENV_DRAM:-dram}

# UB taxonomy
export ENV_MMSEQS=${ENV_MMSEQS:-mmseqs2}

# Activation helper — `source` this in each script then call `activate_env $ENV_NAME`
activate_env() {
    if [ -z "$1" ]; then echo "activate_env: env name required" >&2; return 1; fi
    if [ -z "$CONDA_EXE" ]; then
        source ~/anaconda3/etc/profile.d/conda.sh
    fi
    conda activate "$1"
}
