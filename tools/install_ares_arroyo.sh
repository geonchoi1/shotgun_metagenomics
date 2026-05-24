#!/bin/bash
# === Ares-Arroyo 91 oriT sequences installer ===
# Downloads 91 oriT sequences (Ares-Arroyo et al. 2023 NAR, 10.1093/nar/gkad084)
# from supplementary Table S3, builds BLAST DB.
#
# The 91 sequences are normally distributed as a supplementary table — user
# may need to download manually and place at $INSTALL_DIR/ares_arroyo_91_oriT.fna
# (this script builds the BLAST DB from that FASTA).

set -e
INSTALL_DIR=${INSTALL_DIR:-$HOME/tools/ares_arroyo_oriT}
FASTA=$INSTALL_DIR/ares_arroyo_91_oriT.fna

mkdir -p $INSTALL_DIR

if [ ! -f $FASTA ]; then
    echo "ERROR: $FASTA missing."
    echo "  Get the 91-oriT FASTA from:"
    echo "    Ares-Arroyo et al. 2023 NAR, Table S3 (DOI 10.1093/nar/gkad084)"
    echo "  Place it at $FASTA, then re-run this script."
    exit 1
fi

# Build BLAST DB
if command -v makeblastdb &>/dev/null; then
    makeblastdb -in $FASTA -dbtype nucl
    echo "✓ BLAST DB built at $FASTA.{ndb,nhr,nin,...}"
else
    echo "WARN: makeblastdb not found. Install BLAST+ and re-run, or run:"
    echo "    makeblastdb -in $FASTA -dbtype nucl"
fi
