#!/bin/bash
# === Thread budget config ===
# Override via env vars before invoking any script.

# Default thread budget per single tool invocation
export THREADS=${THREADS:-16}

# Parallel-job count (used for GNU parallel orchestrations: oriTfinder2, per-MAG Bakta, etc.)
export PARALLEL=${PARALLEL:-16}

# Heavy-tool thread overrides
export THREADS_BLAST=${THREADS_BLAST:-$THREADS}
export THREADS_MINIMAP2=${THREADS_MINIMAP2:-$THREADS}
export THREADS_HMMSEARCH=${THREADS_HMMSEARCH:-$THREADS}
export THREADS_BAKTA=${THREADS_BAKTA:-$THREADS}
export THREADS_GTDBTK=${THREADS_GTDBTK:-32}     # pplacer is memory-heavy
export PPLACER_CPUS=${PPLACER_CPUS:-16}

# Memory caps (some tools have flags)
export MEM_KCLUST=${MEM_KCLUST:-35000MB}
export MEM_DIAMOND=${MEM_DIAMOND:-32G}

# Reproducibility
export SEED=${SEED:-42}
