#!/bin/bash
# === oriTfinder2 installer ===
# Downloads oriTfinder2 from SJTU + EMBOSS dependency.
# Installs into $HOME/tools/oriTfinder2_linux/ (overridable).
#
# Note: oriTfinder2 download URL has been unreliable in the past — if the
# primary SJTU mirror fails, manual install from GitHub fork or local backup
# may be required. After install, validate with the example dataset.

set -e
INSTALL_DIR=${INSTALL_DIR:-$HOME/tools/oriTfinder2_linux}

if [ -d "$INSTALL_DIR" ] && [ -x "$INSTALL_DIR/run_oriTfinder2_local_batch.pl" ]; then
    echo "Already installed at $INSTALL_DIR"
    exit 0
fi

mkdir -p $(dirname $INSTALL_DIR)
cd $(dirname $INSTALL_DIR)

# Primary download — SJTU mirror (subject to change)
URL=${URL:-http://bioinfo-mml.sjtu.edu.cn/oriTfinder/oriTfinder2_linux.tar.gz}
echo "Downloading oriTfinder2 from $URL ..."
if ! wget -q --show-progress "$URL" -O oriTfinder2_linux.tar.gz; then
    echo "Primary download failed. Try one of:"
    echo "  1. Get URL from https://tool-mml.sjtu.edu.cn/oriTfinder/"
    echo "  2. Provide local copy at $(dirname $INSTALL_DIR)/oriTfinder2_linux.tar.gz then re-run"
    exit 1
fi

tar xzf oriTfinder2_linux.tar.gz
rm oriTfinder2_linux.tar.gz

# EMBOSS dependency (extractseq + revseq)
echo "Checking EMBOSS extractseq + revseq ..."
if ! command -v extractseq &>/dev/null; then
    echo "EMBOSS not found in PATH. Install via:"
    echo "  conda install -c bioconda emboss"
    echo "Then copy extractseq + revseq into $INSTALL_DIR/tools/"
else
    cp $(which extractseq) $INSTALL_DIR/tools/ 2>/dev/null || true
    cp $(which revseq) $INSTALL_DIR/tools/ 2>/dev/null || true
fi

# Sanity check
if [ ! -f $INSTALL_DIR/run_oriTfinder2_local_batch.pl ]; then
    echo "ERROR: install did not produce run_oriTfinder2_local_batch.pl"
    exit 1
fi

echo "✓ oriTfinder2 installed at $INSTALL_DIR"
echo "  Test: cd $INSTALL_DIR && perl run_oriTfinder2_local_batch.pl list_demo.txt"
