#!/bin/bash
# === Master orchestrator ===
# Runs the entire pipeline: shared upstream → 3 tracks → cross-track.
# Requires: inputs staged + circ_frag_map.tsv provided + config/db_paths.sh + config/envs.sh edited for your environment.

set -e
BASE=$(cd "$(dirname "$0")" && pwd)
source $BASE/config/db_paths.sh
source $BASE/config/envs.sh
source $BASE/config/threads.sh

# --- Read type ---
: ${READ_TYPE:?ERROR: export READ_TYPE=illumina or hifi before running}
echo "[$(date '+%F %T')] === Master pipeline START (READ_TYPE=$READ_TYPE) ==="

# --- 00_shared: reads → bin ---
bash $BASE/00_shared/run.sh

# --- 01-03 tracks (parallel-safe; run sequentially here for clarity) ---
bash $BASE/01_plasmid_track/run.sh
bash $BASE/02_mag_track/run.sh
bash $BASE/03_unbinned_track/run.sh

# --- 04 cross-track ---
bash $BASE/04_cross_track/run.sh

echo "[$(date '+%F %T')] === Master pipeline DONE ==="
