#!/bin/bash
# === KO annotation via direct hmmsearch against KOfam profiles (03_unbinned_track) ===
# Replaces the previous KOfamScan-based step. Annotates ORFs by direct HMMER
# hmmsearch against the concatenated KOfam HMM database (genome.jp/ftp/db/kofam/).
#
# Cutoff modes (KOFAM_CUTOFF_MODE):
#   adaptive (default) — per-KO bit-score threshold from KOfam ko_list (KOfamScan-equivalent).
#                         KOs without an adaptive threshold fall back to KOFAM_EVALUE.
#   uniform            — single user-provided E-value (KOFAM_EVALUE) and optional
#                         bit-score (KOFAM_SCORE) applied to ALL KOs (lenient broad scan).
#
# Environment variables:
#   KOFAM_DB          = /mnt/nas/DB/geon/kofam_db
#   QUERY_FAA         = $PROJECT/03_unbinned_track/03_master_orf/unbinned.master.faa
#   OUT               = $PROJECT/03_unbinned_track/06_kofam
#   THREADS           = 112
#   KOFAM_CUTOFF_MODE = adaptive | uniform     (default adaptive)
#   KOFAM_EVALUE      = 1e-5
#   KOFAM_SCORE       = (empty; uniform mode only)
#   KOFAM_PRE_EVALUE  = 1e-2 (broad pre-filter; per-KO threshold applied below)

set -eo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

KOFAM_DB=${KOFAM_DB:-/mnt/nas/DB/geon/kofam_db}
QUERY_FAA=${QUERY_FAA:-$PROJECT/03_unbinned_track/03_master_orf/unbinned.master.faa}
OUT=${OUT:-$PROJECT/03_unbinned_track/06_kofam}
THREADS=${THREADS:-112}

KOFAM_CUTOFF_MODE=${KOFAM_CUTOFF_MODE:-adaptive}
KOFAM_EVALUE=${KOFAM_EVALUE:-1e-5}
KOFAM_SCORE=${KOFAM_SCORE:-}
KOFAM_PRE_EVALUE=${KOFAM_PRE_EVALUE:-1e-2}

mkdir -p "$OUT"
activate_env "$ENV_DIAMOND"   # base; HMMER from anaconda3 PATH

ALL_HMM="$OUT/kofam_all.hmm"
if [ ! -s "$ALL_HMM.h3p" ]; then
    echo "[$(date +%T)] concat KOfam profiles -> $ALL_HMM"
    find "$KOFAM_DB/profiles" -maxdepth 1 -name "K*.hmm" -print0 | xargs -0 cat > "$ALL_HMM"
    echo "[$(date +%T)] hmmpress..."
    hmmpress -f "$ALL_HMM" > /dev/null
fi

TBLOUT="$OUT/kofam_search.tbl"
if [ ! -s "$TBLOUT" ]; then
    echo "[$(date +%T)] hmmsearch --cpu $THREADS -E $KOFAM_PRE_EVALUE"
    hmmsearch --cpu "$THREADS" -E "$KOFAM_PRE_EVALUE" --noali \
        --tblout "$TBLOUT" "$ALL_HMM" "$QUERY_FAA" > /dev/null
fi

echo "[$(date +%T)] mode=$KOFAM_CUTOFF_MODE   E=$KOFAM_EVALUE   score=${KOFAM_SCORE:-NA}"
python3 - <<PYEOF
import sys
mode   = "$KOFAM_CUTOFF_MODE".lower()
ev_thr = float("$KOFAM_EVALUE")
sc_s   = "$KOFAM_SCORE".strip()
sc_thr = float(sc_s) if sc_s else None

ko_thr = {}
for line in open("$KOFAM_DB/ko_list"):
    p = line.rstrip("\n").split("\t")
    if not p or p[0] == "knum" or len(p) < 3: continue
    knum, thr_s, stype = p[0], p[1], p[2]
    if thr_s == "-" or not thr_s: continue
    try: ko_thr[knum] = (float(thr_s), stype)
    except ValueError: pass

best = {}
n_seen = 0; n_kept = 0
for line in open("$TBLOUT"):
    if line.startswith("#") or not line.strip(): continue
    parts = line.split()
    if len(parts) < 9: continue
    n_seen += 1
    orf, ko = parts[0], parts[2]
    try:
        evalue_full = float(parts[4]); score_full = float(parts[5])
        score_dom   = float(parts[8])
    except ValueError: continue

    if mode == "adaptive":
        if ko in ko_thr:
            thr, stype = ko_thr[ko]
            sc = score_full if stype == "full" else score_dom
            keep = (sc >= thr)
        else:
            keep = (evalue_full <= ev_thr)
    elif mode == "uniform":
        keep = (evalue_full <= ev_thr) and (sc_thr is None or score_full >= sc_thr)
    else:
        sys.exit(f"unknown KOFAM_CUTOFF_MODE: {mode}")

    if not keep: continue
    n_kept += 1
    if orf not in best or evalue_full < best[orf][0]:
        best[orf] = (evalue_full, score_full, ko)

with open("$OUT/kofam_mapper.tsv", "w") as o:
    for orf in sorted(best):
        o.write(f"{orf}\t{best[orf][2]}\n")

kc = {}
for _, (_, _, k) in best.items(): kc[k] = kc.get(k, 0) + 1
print(f"  pre-filter rows: {n_seen}")
print(f"  passed cutoff  : {n_kept}")
print(f"  ORFs annotated : {len(best)}")
print(f"  unique KOs     : {len(kc)}")
PYEOF
echo "[$(date +%T)] DONE -> $OUT/kofam_mapper.tsv"
