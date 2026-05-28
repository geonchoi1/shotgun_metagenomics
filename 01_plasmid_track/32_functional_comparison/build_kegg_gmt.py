#!/usr/bin/env python3
"""Build KEGG pathway GMT file from KEGG REST API.

Output format (GMT — gene set MAtrix Transposed):
    pathway_id <TAB> description <TAB> KO1 <TAB> KO2 <TAB> ... <TAB> KOn

Used by gseapy.prerank in run.sh.

REST endpoints:
    https://rest.kegg.jp/list/pathway        → pathway_id <TAB> description (per line)
    https://rest.kegg.jp/link/ko/<map_id>    → map_id <TAB> KO (one line per KO)

Usage:
    python build_kegg_gmt.py --out $DB_ROOT/gsea/kegg_pathway.gmt
    # or filtered to bacterial / prokaryotic only:
    python build_kegg_gmt.py --prefix map --out kegg_pathway.gmt

The script caches pathway→KO mappings in a single bulk REST call
(https://rest.kegg.jp/link/ko/pathway) for efficiency (~2 MB, ~30 sec).
"""
import argparse
import sys
import os
from collections import defaultdict
from urllib.request import urlopen
from urllib.error import URLError


def fetch(url, timeout=120):
    """GET from KEGG REST. Returns text body."""
    try:
        with urlopen(url, timeout=timeout) as r:
            return r.read().decode('utf-8')
    except URLError as e:
        print(f"ERROR fetching {url}: {e}", file=sys.stderr)
        sys.exit(1)


def build_gmt(out_path, prefix_filter='map', skip_global=True, skip_brite=True):
    """Build KEGG pathway → KO mapping GMT.

    Args:
        out_path: output .gmt path
        prefix_filter: only keep pathway IDs starting with this prefix
                       'map' = reference pathway (no organism prefix)
                       'ko'  = same as map but KO-centric
        skip_global: skip "global" maps (map01100, map01110, etc.) — too generic
        skip_brite: skip BRITE hierarchies (br...)
    """
    print(f"[1/3] Fetching pathway list from KEGG REST API...")
    pathway_list = fetch('https://rest.kegg.jp/list/pathway')
    # Format: map00010 <TAB> Glycolysis / Gluconeogenesis
    pathways = {}
    for line in pathway_list.strip().split('\n'):
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) < 2:
            continue
        pid, desc = parts[0].strip(), parts[1].strip()
        # Filter
        if not pid.startswith(prefix_filter):
            continue
        if skip_brite and pid.startswith('br'):
            continue
        # Skip global maps (overview)
        if skip_global and pid.startswith(prefix_filter + '0110'):
            continue
        if skip_global and pid.startswith(prefix_filter + '01100'):
            continue
        pathways[pid] = desc
    print(f"      pathways retained: {len(pathways)}")

    print(f"[2/3] Fetching pathway → KO links (bulk)...")
    link_text = fetch('https://rest.kegg.jp/link/ko/pathway')
    # Format: path:map00010 <TAB> ko:K00844
    pathway_to_ko = defaultdict(set)
    for line in link_text.strip().split('\n'):
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) < 2:
            continue
        pid = parts[0].replace('path:', '').strip()
        ko = parts[1].replace('ko:', '').strip()
        if pid in pathways:
            pathway_to_ko[pid].add(ko)

    n_with_kos = sum(1 for v in pathway_to_ko.values() if v)
    print(f"      pathways with ≥1 KO: {n_with_kos}")

    print(f"[3/3] Writing GMT → {out_path}")
    os.makedirs(os.path.dirname(out_path) if os.path.dirname(out_path) else '.', exist_ok=True)
    n_written = 0
    with open(out_path, 'w') as out:
        for pid, desc in sorted(pathways.items()):
            kos = sorted(pathway_to_ko.get(pid, []))
            if not kos:
                continue
            # GMT: pathway_id <TAB> description <TAB> KO1 <TAB> KO2 ...
            out.write(f"{pid}\t{desc}\t" + "\t".join(kos) + "\n")
            n_written += 1
    print(f"DONE. Wrote {n_written} pathways × KO sets to {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True, help='Output GMT path (e.g. /path/to/kegg_pathway.gmt)')
    ap.add_argument('--prefix', default='map', choices=['map', 'ko'],
                    help='Pathway ID prefix to keep ("map" = reference [default], "ko" = ko-centric duplicate)')
    ap.add_argument('--keep-global', action='store_true',
                    help='Keep global maps (map01100 etc.) — default skip')
    args = ap.parse_args()

    build_gmt(args.out,
              prefix_filter=args.prefix,
              skip_global=not args.keep_global)


if __name__ == '__main__':
    main()
