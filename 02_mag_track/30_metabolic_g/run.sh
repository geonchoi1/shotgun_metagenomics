#!/bin/bash
# === 30 METABOLIC-G on dereplicated MAGs ===
# Input:  $PROJECT/00_shared/07_mag_production/08_drep_species/dereplicated_genomes/*.fa
#         (re-linked as .fasta because METABOLIC requires .fasta extension)
# Output: $PROJECT/02_mag_track/30_metabolic_g/

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

SRC=$PROJECT/00_shared/07_mag_production/08_drep_species/dereplicated_genomes
OUT=$PROJECT/02_mag_track/30_metabolic_g
STAGE=$OUT/_input_fasta
mkdir -p "$OUT" "$STAGE"

[ -d "$SRC" ] || { echo "ERROR: $SRC missing" >&2; exit 1; }

# Stage .fasta symlinks (METABOLIC won't accept .fa)
echo "[$(date '+%F %T')] staging .fasta symlinks"
for fa in "$SRC"/*.fa "$SRC"/*.fna "$SRC"/*.fasta; do
    [ -f "$fa" ] || continue
    base=$(basename "$fa"); name="${base%.fa}"; name="${name%.fna}"; name="${name%.fasta}"
    ln -sf "$fa" "$STAGE/${name}.fasta"
done

activate_env "$ENV_METABOLIC"

# Discover the perl script (METABOLIC install layout: $CONDA_PREFIX/share/METABOLIC*/METABOLIC-G.pl)
METABOLIC_PL=${METABOLIC_PL:-$(command -v METABOLIC-G.pl 2>/dev/null || true)}
if [ -z "$METABOLIC_PL" ]; then
    METABOLIC_PL=$(find "${CONDA_PREFIX:-/opt}" -maxdepth 6 -name 'METABOLIC-G.pl' 2>/dev/null | head -1)
fi
[ -n "$METABOLIC_PL" ] || { echo "ERROR: METABOLIC-G.pl not found (set METABOLIC_PL=)" >&2; exit 1; }

echo "[$(date '+%F %T')] METABOLIC-G.pl"
perl "$METABOLIC_PL" \
    -in-gn "$STAGE" \
    -t "$THREADS" \
    -o "$OUT"

echo "[$(date '+%F %T')] DONE — output in $OUT"
