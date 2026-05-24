#!/bin/bash
# ==============================================================================
# Shotgun metagenomics pipeline — full master
# Runs all 5 sections: initial → plasmid + mag + unbinned (parallel-safe) → cross
# ==============================================================================
# Usage:
#   READ_TYPE=hifi bash run.sh /path/to/project
#   READ_TYPE=illumina bash run.sh /path/to/project
#
# Skips:
#   SKIP_INITIAL=1   skip 00_shared
#   SKIP_PLASMID=1
#   SKIP_MAG=1
#   SKIP_UNBINNED=1
#   SKIP_CROSS=1
#
# Env overrides (see config.sh for full list):
#   THREADS=32 bash run.sh /path/to/project
# ==============================================================================

set -euo pipefail
REPO="$(cd "$(dirname "$0")" && pwd)"
source "$REPO/config.sh"

PROJECT=${1:?ERROR: usage: bash run.sh /path/to/project}
PROJECT=$(readlink -f "$PROJECT")
export PROJECT REPO

: ${READ_TYPE:?ERROR: export READ_TYPE=illumina or hifi}

echo "=========================================="
echo " shotgun_metagenomics — full pipeline"
echo "   PROJECT  : $PROJECT"
echo "   REPO     : $REPO"
echo "   READ_TYPE: $READ_TYPE"
echo "   THREADS  : $THREADS"
echo "=========================================="

[ "${SKIP_INITIAL:-0}" = "1" ] || bash "$REPO/run_initial.sh"  "$PROJECT"
[ "${SKIP_PLASMID:-0}" = "1" ] || bash "$REPO/run_plasmid.sh"  "$PROJECT"
[ "${SKIP_MAG:-0}"     = "1" ] || bash "$REPO/run_mag.sh"      "$PROJECT"
[ "${SKIP_UNBINNED:-0}"= "1" ] || bash "$REPO/run_unbinned.sh" "$PROJECT"
[ "${SKIP_CROSS:-0}"   = "1" ] || bash "$REPO/run_cross.sh"    "$PROJECT"

echo "[$(date '+%F %T')] ALL DONE"
