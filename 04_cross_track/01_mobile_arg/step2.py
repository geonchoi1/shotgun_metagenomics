#!/usr/bin/env python3
"""
Step 2: ARG × MGE presence/absence matrix per source.

For each source (plasmid / mag / unbinned), produces a contig-level matrix
where rows = contigs and columns = ARG genes (from AMRFinder) plus MGE classes
(IS_*, integron, ICE). Cell = 1 if that feature is present on the contig.

Output: $PROJECT/04_cross_track/mobile_arg/step2/{plasmid,mag,unbinned}_arg_mge_matrix.tsv
"""
import os
import sys
from pathlib import Path
from collections import defaultdict
import csv

PROJECT = os.environ.get("PROJECT")
if not PROJECT:
    sys.exit("ERROR: export PROJECT=...")

OUT = Path(PROJECT) / "cross/mobile_arg/step2"
OUT.mkdir(parents=True, exist_ok=True)

SOURCES = {
    "plasmid":  Path(PROJECT) / "plasmid",
    "mag":      Path(PROJECT) / "mag",
    "unbinned": Path(PROJECT) / "unbinned",
}

# Step-index per track (plasmid uses 09/20/21/22; MAG/UB use 08/19/20/21)
PATHS = {
    "plasmid":  dict(amr="09_amrfinder/all/amrfinder.tsv",
                     ise="20_isescan/all",
                     intf="21_integronfinder/all",
                     ice="22_iceberg3/all/iceberg3.tsv",
                     orf2c="04_master_orf/all/orf2contig.tsv"),
    "mag":      dict(amr="08_amrfinder/all/amrfinder.tsv",
                     ise="19_isescan/all",
                     intf="20_integronfinder/all",
                     ice="21_iceberg3/all/iceberg3.tsv",
                     orf2c="03_master_orf/all/orf2contig.tsv"),
    "unbinned": dict(amr="08_amrfinder/all/amrfinder.tsv",
                     ise="19_isescan/all",
                     intf="20_integronfinder/all",
                     ice="21_iceberg3/all/iceberg3.tsv",
                     orf2c="03_master_orf/all/orf2contig.tsv"),
}


def load_orf2contig(p):
    o2c = {}
    if not p.exists():
        return o2c
    with p.open() as fh:
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if len(a) >= 2:
                o2c[a[0]] = a[1]
    return o2c


def parse_amrfinder(p, o2c):
    """Return {contig: set(gene_symbol)}"""
    out = defaultdict(set)
    if not p.exists() or p.stat().st_size == 0:
        return out
    with p.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_orf = header.index("Protein identifier")
        except ValueError:
            i_orf = 0
        try:
            i_sym = header.index("Gene symbol")
        except ValueError:
            i_sym = 5
        try:
            i_contig = header.index("Contig id")
        except ValueError:
            i_contig = -1
        try:
            i_class = header.index("Class")
        except ValueError:
            i_class = -1
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if len(a) <= max(i_orf, i_sym):
                continue
            contig = a[i_contig] if i_contig >= 0 and i_contig < len(a) and a[i_contig] else o2c.get(a[i_orf], "")
            sym = a[i_sym].strip()
            if not contig or not sym:
                continue
            # Restrict to AMR class if column exists
            if i_class >= 0 and i_class < len(a):
                cls = a[i_class].upper()
                if cls and cls not in ("AMR",):
                    continue
            out[contig].add(f"ARG:{sym}")
    return out


def parse_isescan(d):
    """Each *.tsv in isescan output dir has columns: seqID, family, ... (header)"""
    out = defaultdict(set)
    if not d.exists():
        return out
    for f in d.glob("*.tsv"):
        with f.open() as fh:
            head = fh.readline().rstrip("\n").split("\t")
            try:
                i_seq = head.index("seqID")
            except ValueError:
                i_seq = 0
            try:
                i_fam = head.index("family")
            except ValueError:
                i_fam = 1
            for line in fh:
                a = line.rstrip("\n").split("\t")
                if len(a) <= max(i_seq, i_fam):
                    continue
                out[a[i_seq]].add(f"IS:{a[i_fam]}")
    return out


def parse_integronfinder(d):
    out = defaultdict(set)
    if not d.exists():
        return out
    for summ in d.rglob("*.summary"):
        with summ.open() as fh:
            head = fh.readline().rstrip("\n").split("\t")
            try:
                i_id = head.index("ID_replicon")
            except ValueError:
                i_id = 0
            for line in fh:
                a = line.rstrip("\n").split("\t")
                if len(a) <= i_id:
                    continue
                out[a[i_id]].add("INT:integron")
    return out


def parse_iceberg(p, o2c):
    out = defaultdict(set)
    if not p.exists() or p.stat().st_size == 0:
        return out
    with p.open() as fh:
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if len(a) < 2:
                continue
            contig = o2c.get(a[0], "")
            if contig:
                out[contig].add("ICE:iceberg3")
    return out


def build(src, root):
    cfg = PATHS[src]
    o2c = load_orf2contig(root / cfg["orf2c"])
    arg = parse_amrfinder(root / cfg["amr"], o2c)
    ise = parse_isescan(root / cfg["ise"])
    intf = parse_integronfinder(root / cfg["intf"])
    ice = parse_iceberg(root / cfg["ice"], o2c)

    feats = defaultdict(set)
    for d in (arg, ise, intf, ice):
        for c, sset in d.items():
            feats[c] |= sset

    all_feats = sorted({f for s in feats.values() for f in s})
    out_tsv = OUT / f"{src}_arg_mge_matrix.tsv"
    with out_tsv.open("w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["contig"] + all_feats)
        for c in sorted(feats):
            w.writerow([c] + [1 if f in feats[c] else 0 for f in all_feats])
    print(f"  {src}: contigs={len(feats)} features={len(all_feats)} -> {out_tsv}")


def main():
    for src, root in SOURCES.items():
        if not root.exists():
            print(f"  WARN: {src} root missing: {root}")
            continue
        build(src, root)
    print("[step2] DONE")


if __name__ == "__main__":
    main()
