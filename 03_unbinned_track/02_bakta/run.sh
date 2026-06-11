#!/bin/bash
# === 02 Bakta on unbinned chromosomal contigs ===
# Splits into circ/linear (matches plasmid + MAG pattern). Linear can be very large
# (~100k+ contigs); pyrodigal.train() OOMs above ~100k. Linear is auto-chunked +
# parallelized.
#
# Input:  $PROJECT/03_unbinned_track/01_raw_fasta/{circ,frag}/unbinned.fna
# Output: $PROJECT/03_unbinned_track/02_bakta/{circ,frag}/unbinned_{circ,frag}.{faa,ffn,fna,gff3,tsv}
#
# Env overrides:
#   CHUNK_SIZE=10000          # contigs per chunk for linear
#   BAKTA_CHUNK_PARALLEL=16   # parallel chunks for linear
#   THREADS_BAKTA_CHUNK=8     # threads per chunk
#   Total cores = PARALLEL × THREADS_BAKTA_CHUNK (default 16×8 = 128)

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

CIRC_IN=$PROJECT/03_unbinned_track/01_raw_fasta/circ/unbinned.fna
FRAG_IN=$PROJECT/03_unbinned_track/01_raw_fasta/frag/unbinned.fna
OUT=$PROJECT/03_unbinned_track/02_bakta
mkdir -p "$OUT"/{circ,frag,frag_chunks}

CHUNK_SIZE=${CHUNK_SIZE:-10000}
PARALLEL=${BAKTA_CHUNK_PARALLEL:-16}
T_CHUNK=${THREADS_BAKTA_CHUNK:-8}

activate_env "$ENV_BAKTA"

# === circ Bakta (--complete) ===
if [ -s "$CIRC_IN" ] && [ ! -s "$OUT/circ/unbinned_circ.gff3" ]; then
    echo "[$(date '+%F %T')] Bakta UB circ (--complete)"
    bakta --db "$BAKTA_DB" \
          --output "$OUT/circ" \
          --prefix unbinned_circ \
          --keep-contig-headers \
          --threads "$THREADS_BAKTA" \
          --complete \
          --skip-crispr --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region --skip-gap \
          --force \
          "$CIRC_IN" > "$OUT/circ/bakta.log" 2>&1
    echo "  circ ORF: $(grep -c '^>' $OUT/circ/unbinned_circ.faa)"
else
    echo "[$(date '+%F %T')] UB circ Bakta — skip (already done or no circ input)"
fi

# === linear Bakta — CHUNKED parallel ===
if [ -s "$FRAG_IN" ] && [ ! -s "$OUT/frag/unbinned_frag.faa" ]; then
    echo "[$(date '+%F %T')] UB linear Bakta — chunking (size=$CHUNK_SIZE)"

    # Step 1: Split linear/unbinned.fna into chunks
    rm -f "$OUT/frag_chunks"/chunk_*.fna
    python3 - <<PYEOF
total = sum(1 for line in open("$FRAG_IN") if line.startswith(">"))
print(f"  Total linear contigs: {total}")
N = max(1, (total + $CHUNK_SIZE - 1) // $CHUNK_SIZE)
chunk_size = (total + N - 1) // N
idx = 0; chunk = 0
out = open(f"$OUT/frag_chunks/chunk_{chunk:02d}.fna", "w")
for line in open("$FRAG_IN"):
    if line.startswith(">"):
        if idx >= chunk_size:
            out.close(); chunk += 1; idx = 0
            out = open(f"$OUT/frag_chunks/chunk_{chunk:02d}.fna", "w")
        idx += 1
    out.write(line)
out.close()
print(f"  Created {chunk+1} chunks of ~{chunk_size} contig each")
PYEOF

    # Step 2: Run Bakta in parallel on chunks
    echo "[$(date '+%F %T')] Running ${PARALLEL}-way × ${T_CHUNK}-thread Bakta"
    BAKTA_RUN() {
        local C=$1
        local CHUNK_OUT="$OUT/frag_chunks/out_$C"
        [ -s "$CHUNK_OUT/${C}.faa" ] && return 0
        bakta --db "$BAKTA_DB" \
              --output "$CHUNK_OUT" \
              --prefix "$C" \
              --keep-contig-headers \
              --threads $T_CHUNK \
              --skip-crispr --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region --skip-gap \
              --force \
              "$OUT/frag_chunks/${C}.fna" > "$OUT/frag_chunks/${C}.runlog" 2>&1
    }
    export -f BAKTA_RUN
    export OUT T_CHUNK BAKTA_DB

    ls "$OUT/frag_chunks"/chunk_*.fna | xargs -n1 basename | sed 's/\.fna$//' | \
      parallel -j $PARALLEL --joblog "$OUT/frag_chunks/parallel.joblog" BAKTA_RUN {}

    # Step 3: Merge chunk outputs → frag/unbinned_frag.*
    echo "[$(date '+%F %T')] Merging chunk outputs"
    > "$OUT/frag/unbinned_frag.faa"
    > "$OUT/frag/unbinned_frag.ffn"
    : > "$OUT/frag/unbinned_frag.gff3"
    first=1
    for C in $(ls "$OUT/frag_chunks"/chunk_*.fna | xargs -n1 basename | sed 's/\.fna$//' | sort); do
        CHUNK_OUT="$OUT/frag_chunks/out_$C"
        [ -s "$CHUNK_OUT/${C}.faa" ] || { echo "  WARN: $C missing output"; continue; }
        cat "$CHUNK_OUT/${C}.faa" >> "$OUT/frag/unbinned_frag.faa"
        cat "$CHUNK_OUT/${C}.ffn" >> "$OUT/frag/unbinned_frag.ffn"
        if [ "$first" = "1" ]; then
            cp "$CHUNK_OUT/${C}.gff3" "$OUT/frag/unbinned_frag.gff3"
            first=0
        else
            grep -v "^#" "$CHUNK_OUT/${C}.gff3" >> "$OUT/frag/unbinned_frag.gff3"
        fi
    done
    echo "  linear ORF: $(grep -c '^>' $OUT/frag/unbinned_frag.faa)"
else
    echo "[$(date '+%F %T')] UB linear Bakta — skip (already done or no input)"
fi

echo "[$(date '+%F %T')] DONE"
echo "  circ:   $OUT/circ/unbinned_circ.{faa,gff3,ffn}"
echo "  linear: $OUT/frag/unbinned_frag.{faa,gff3,ffn}"
