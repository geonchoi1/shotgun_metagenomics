#!/bin/bash
# === 17 BiG-SCAPE (BGC clustering) on UB antiSMASH gbk outputs ===
# Input:  $PROJECT/unbinned/16_antismash/all/ (recursive .region*.gbk)
# Output: $PROJECT/unbinned/17_bigscape/all/

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

ASM_DIR=$PROJECT/unbinned/16_antismash/all
OUT=$PROJECT/unbinned/17_bigscape/all
GBK_DIR=$OUT/gbk_input
mkdir -p "$OUT" "$GBK_DIR"

n_gbk=$(find "$ASM_DIR" -name "*.region*.gbk" 2>/dev/null | wc -l)
if [ "$n_gbk" -eq 0 ]; then
    echo "[$(date '+%F %T')] no antiSMASH region gbk found in $ASM_DIR — skip"
    exit 0
fi

if [ -d "$OUT/network_files" ] && [ "$(ls -A "$OUT/network_files" 2>/dev/null)" ]; then
    echo "[$(date '+%F %T')] UB BiG-SCAPE already done — skip"; exit 0
fi

# Stage region gbks into single directory (BiG-SCAPE prefers flat input)
find "$ASM_DIR" -name "*.region*.gbk" -exec ln -sf {} "$GBK_DIR"/ \;

activate_env "$ENV_BIGSCAPE"

echo "[$(date '+%F %T')] bigscape on $n_gbk UB BGCs"
bigscape \
    -i "$GBK_DIR" \
    -o "$OUT" \
    --cpus "$THREADS" \
    --mibig \
    --mix \
    --include_singletons \
    --hybrids-off > "$OUT/bigscape.log" 2>&1

echo "[$(date '+%F %T')] DONE"
