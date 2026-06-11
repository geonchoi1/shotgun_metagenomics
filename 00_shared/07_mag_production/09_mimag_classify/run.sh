#!/bin/bash
# === 07.09 MIMAG classification per MAG ===
# HQ: completeness>=90, contamination<5, has 5S+16S+23S (bac OR arc kingdom),
#     >=18 distinct tRNA types
# NC: completeness>=90, contamination<5, otherwise
# MQ: completeness>=50, contamination<10
# LQ: else
#
# Input:  $PROJECT/00_shared/07_mag_production/05_checkm2/quality_report.tsv
#         $PROJECT/00_shared/07_mag_production/06_rrna/<MAG>.{bac,arc}.gff
#         $PROJECT/00_shared/07_mag_production/07_trna/<MAG>.{bac,arc}.tsv
# Output: $PROJECT/00_shared/07_mag_production/09_mimag.tsv

set -euo pipefail
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$SCRIPT_DIR/../../.." && pwd)
source "$REPO/config.sh"
: ${PROJECT:?ERROR: export PROJECT=/path/to/project}

CHECKM=$PROJECT/00_shared/07_mag_production/05_checkm2/quality_report.tsv
RRNA=$PROJECT/00_shared/07_mag_production/06_rrna
TRNA=$PROJECT/00_shared/07_mag_production/07_trna
OUT=$PROJECT/00_shared/07_mag_production/09_mimag.tsv

[ -s "$CHECKM" ] || { echo "ERROR: $CHECKM missing" >&2; exit 2; }

if [ -s "$OUT" ]; then
    echo "[$(date '+%F %T')] 09 MIMAG already done"; exit 0
fi

echo "[$(date '+%F %T')] classify MIMAG"
python - <<PY
import os, re, csv, sys

checkm = "${CHECKM}"
rrna_dir = "${RRNA}"
trna_dir = "${TRNA}"
out = "${OUT}"

def parse_rrna_gff(path):
    """Return set of {5S,16S,23S} present (across bac+arc unioned by caller)."""
    found = set()
    if not (os.path.exists(path) and os.path.getsize(path) > 0):
        return found
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attr = cols[8]
            m = re.search(r"Name=([0-9.]+S)_rRNA", attr)
            if m:
                rrna = m.group(1)
                rrna = rrna.split(".")[0] + "S"  # normalize "16S" / "5S" / "23S"
                if rrna in ("5S", "16S", "23S"):
                    found.add(rrna)
    return found

def parse_trna_tsv(path):
    """Return set of distinct anticodon-derived amino-acid types from tRNAscan-SE."""
    types = set()
    if not (os.path.exists(path) and os.path.getsize(path) > 0):
        return types
    with open(path) as fh:
        # tRNAscan-SE output has 3 header lines, then columns:
        # Sequence_Name tRNA# Begin End Type Codon Begin End Score Note
        for i, line in enumerate(fh):
            if i < 3:
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            aa = parts[4]
            if aa and aa not in ("Undet", "Sup", "Pseudo"):
                types.add(aa)
    return types

# Read CheckM2 report
mags = {}
with open(checkm) as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    for row in rdr:
        name = row["Name"]
        try:
            comp = float(row["Completeness"])
            cont = float(row["Contamination"])
        except Exception:
            continue
        mags[name] = {"completeness": comp, "contamination": cont}

with open(out, "w") as out_fh:
    w = csv.writer(out_fh, delimiter="\t")
    w.writerow(["mag", "completeness", "contamination", "rrna_5S_16S_23S",
                "n_trna_types", "mimag_quality"])
    for mag, info in sorted(mags.items()):
        rrna_bac = parse_rrna_gff(os.path.join(rrna_dir, f"{mag}.bac.gff"))
        rrna_arc = parse_rrna_gff(os.path.join(rrna_dir, f"{mag}.arc.gff"))
        rrna = rrna_bac | rrna_arc
        trna_bac = parse_trna_tsv(os.path.join(trna_dir, f"{mag}.bac.tsv"))
        trna_arc = parse_trna_tsv(os.path.join(trna_dir, f"{mag}.arc.tsv"))
        trna = trna_bac | trna_arc

        rrna_complete = {"5S","16S","23S"}.issubset(rrna)
        n_trna = len(trna)
        comp = info["completeness"]
        cont = info["contamination"]

        if comp >= 90 and cont < 5 and rrna_complete and n_trna >= 18:
            q = "HQ"
        elif comp >= 90 and cont < 5:
            q = "NC"
        elif comp >= 50 and cont < 10:
            q = "MQ"
        else:
            q = "LQ"

        w.writerow([mag, f"{comp:.2f}", f"{cont:.2f}",
                    ",".join(sorted(rrna)) or "-", n_trna, q])

print(f"wrote {out}")
PY

echo "[$(date '+%F %T')] DONE"
