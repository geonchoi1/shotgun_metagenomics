#!/bin/bash
# === 17 BiG-SCAPE on antiSMASH GBKs ===
# Input:  $PROJECT/mag/16_antismash/*/  (region GBKs)
# Output: $PROJECT/mag/17_bigscape/network_files, html_content

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

ANT_BASE=$PROJECT/mag/16_antismash
GBK_STAGE=$PROJECT/mag/17_bigscape/gbk_in
OUT=$PROJECT/mag/17_bigscape
mkdir -p "$GBK_STAGE" "$OUT"

[ -d "$ANT_BASE" ] || { echo "ERROR: $ANT_BASE missing — run 16_antismash first" >&2; exit 1; }

# Stage region GBKs (region*.gbk) into single dir
echo "[$(date '+%F %T')] staging region GBKs"
n_gbk=0
for d in "$ANT_BASE"/*/; do
    mag=$(basename "$d")
    for gbk in "$d"/*region*.gbk; do
        [ -f "$gbk" ] || continue
        ln -sf "$gbk" "$GBK_STAGE/${mag}__$(basename "$gbk")"
        n_gbk=$((n_gbk+1))
    done
done
echo "  staged $n_gbk region GBKs"

activate_env "$ENV_BIGSCAPE"

echo "[$(date '+%F %T')] bigscape"
bigscape cluster \
    -i "$GBK_STAGE" \
    -o "$OUT" \
    --cores "$THREADS" \
    --mibig

echo "[$(date '+%F %T')] DONE"
