#!/bin/bash
# === 04_cross_track/01_mobile_arg : Mobile ARG pipeline orchestrator ===
# Spans plasmid (PL) + MAG (MAG) + unbinned (UB) tracks.
# Output root: $PROJECT/04_cross_track/mobile_arg/{step1,..,step8}/
# Each step has its own script in this dir. Re-run is idempotent.

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

export PROJECT REPO SCRIPT_DIR
OUT_ROOT=$PROJECT/04_cross_track/mobile_arg
mkdir -p "$OUT_ROOT"

run_step() {
    local label=$1
    local script=$2
    echo
    echo "############################################################"
    echo "## $label — $script"
    echo "############################################################"
    [ -x "$SCRIPT_DIR/$script" ] || chmod +x "$SCRIPT_DIR/$script" 2>/dev/null || true
    case "$script" in
        *.sh)  bash    "$SCRIPT_DIR/$script" ;;
        *.py)  python  "$SCRIPT_DIR/$script" ;;
        *.R)   Rscript "$SCRIPT_DIR/$script" ;;
        *)     echo "ERROR: unknown script type: $script" >&2; exit 1 ;;
    esac
}

run_step "Step 1: MGE annotation presence check"               step1.sh
run_step "Step 2: ARG × MGE presence/absence matrix per source" step2.py
run_step "Step 3: CooccurrenceAffinity (R)"                     step3.R
run_step "Step 4: ARG-MGE coordinate distance tiers (Bakta GFF)" step4.py
run_step "Step 5: Module extraction (ARG ± 5kb with MGE)"        step5.py
run_step "Step 6: Cross-source BLAST (modules vs PL∪MAG∪UB)"     step6.sh
run_step "Step 7: Mobility pathway classification"               step7.py
run_step "Step 8: ARG mobility network (nodes.tsv + edges.tsv)"  step8.py

echo
echo "[$(date '+%F %T')] Mobile ARG pipeline DONE — output: $OUT_ROOT"
