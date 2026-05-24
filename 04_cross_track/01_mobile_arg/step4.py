#!/usr/bin/env python3
"""
Step 4: ARG–MGE coordinate distance tiers from Bakta GFF.

For each source, build per-contig coordinate intervals of ARGs (AMRFinder)
and MGEs (ISEScan + IntegronFinder + ICEberg3 mapped via GFF/orf2contig).
Compute pairwise distances and classify into 3 tiers:
  strict   <= 2 kb
  standard <= 5 kb
  broad    <= 10 kb
Output: $PROJECT/cross/mobile_arg/step4/{src}_arg_mge_distance.tsv
"""
import os
import sys
from pathlib import Path
from collections import defaultdict

PROJECT = os.environ.get("PROJECT")
if not PROJECT:
    sys.exit("ERROR: export PROJECT=...")

OUT = Path(PROJECT) / "cross/mobile_arg/step4"
OUT.mkdir(parents=True, exist_ok=True)

SOURCES = ["plasmid", "mag", "unbinned"]
ROOTS = {s: Path(PROJECT) / (s if s != "plasmid" else "plasmid") for s in SOURCES}

PATHS = {
    "plasmid":  dict(gff="03_bakta",       # per-genome dirs *.gff3 — fall back to master
                     master_gff="04_master_orf/all/master.gff",
                     amr="09_amrfinder/all/amrfinder.tsv",
                     ise="20_isescan/all",
                     intf="21_integronfinder/all",
                     ice="22_iceberg3/all/iceberg3.tsv",
                     orf2c="04_master_orf/all/orf2contig.tsv"),
    "mag":      dict(master_gff="03_master_orf/all/master.gff",
                     amr="08_amrfinder/all/amrfinder.tsv",
                     ise="19_isescan/all",
                     intf="20_integronfinder/all",
                     ice="21_iceberg3/all/iceberg3.tsv",
                     orf2c="03_master_orf/all/orf2contig.tsv"),
    "unbinned": dict(master_gff="03_master_orf/all/master.gff",
                     amr="08_amrfinder/all/amrfinder.tsv",
                     ise="19_isescan/all",
                     intf="20_integronfinder/all",
                     ice="21_iceberg3/all/iceberg3.tsv",
                     orf2c="03_master_orf/all/orf2contig.tsv"),
}


def parse_gff_cds(p):
    """Return {orf_id: (contig, start, end, strand)}"""
    out = {}
    if not p.exists():
        return out
    with p.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            a = line.rstrip("\n").split("\t")
            if len(a) < 9 or a[2] != "CDS":
                continue
            attrs = {kv.split("=", 1)[0]: kv.split("=", 1)[1]
                     for kv in a[8].split(";") if "=" in kv}
            oid = attrs.get("ID", "")
            if not oid:
                continue
            out[oid] = (a[0], int(a[3]), int(a[4]), a[6])
    return out


def load_orf2contig(p):
    o2c = {}
    if p.exists():
        with p.open() as fh:
            for line in fh:
                a = line.rstrip("\n").split("\t")
                if len(a) >= 2:
                    o2c[a[0]] = a[1]
    return o2c


def collect_args(amr_p, gff_index, o2c):
    """Return [(contig, start, end, label)]"""
    out = []
    if not amr_p.exists() or amr_p.stat().st_size == 0:
        return out
    with amr_p.open() as fh:
        head = fh.readline().rstrip("\n").split("\t")
        i_orf = head.index("Protein identifier") if "Protein identifier" in head else 0
        i_sym = head.index("Gene symbol") if "Gene symbol" in head else 5
        i_class = head.index("Class") if "Class" in head else -1
        i_contig = head.index("Contig id") if "Contig id" in head else -1
        i_start = head.index("Start") if "Start" in head else -1
        i_stop  = head.index("Stop")  if "Stop"  in head else -1
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if len(a) <= i_sym:
                continue
            if i_class >= 0 and a[i_class] and a[i_class].upper() not in ("AMR",):
                continue
            sym = a[i_sym]
            if i_contig >= 0 and i_start >= 0 and i_stop >= 0 \
                    and a[i_contig] and a[i_start] and a[i_stop]:
                try:
                    out.append((a[i_contig], int(a[i_start]),
                                int(a[i_stop]), f"ARG:{sym}"))
                    continue
                except ValueError:
                    pass
            # fallback: ORF coordinates from GFF
            orf = a[i_orf]
            if orf in gff_index:
                c, s, e, _ = gff_index[orf]
                out.append((c, s, e, f"ARG:{sym}"))
    return out


def collect_isescan(d):
    out = []
    if not d.exists():
        return out
    for f in d.glob("*.tsv"):
        with f.open() as fh:
            head = fh.readline().rstrip("\n").split("\t")
            try:
                i_seq = head.index("seqID")
                i_fam = head.index("family")
                i_s = head.index("isBegin")
                i_e = head.index("isEnd")
            except ValueError:
                continue
            for line in fh:
                a = line.rstrip("\n").split("\t")
                if len(a) <= max(i_seq, i_fam, i_s, i_e):
                    continue
                try:
                    out.append((a[i_seq], int(a[i_s]), int(a[i_e]),
                                f"IS:{a[i_fam]}"))
                except ValueError:
                    continue
    return out


def collect_integron(d):
    out = []
    if not d.exists():
        return out
    for it in d.rglob("*.integrons"):
        with it.open() as fh:
            head_line = fh.readline()
            while head_line.startswith("#"):
                head_line = fh.readline()
            head = head_line.rstrip("\n").split("\t")
            try:
                i_id = head.index("ID_replicon")
                i_s = head.index("pos_beg")
                i_e = head.index("pos_end")
            except ValueError:
                continue
            for line in fh:
                a = line.rstrip("\n").split("\t")
                if len(a) <= max(i_id, i_s, i_e):
                    continue
                try:
                    out.append((a[i_id], int(a[i_s]), int(a[i_e]),
                                "INT:integron"))
                except ValueError:
                    continue
    return out


def collect_iceberg(p, gff_index, o2c):
    out = []
    if not p.exists() or p.stat().st_size == 0:
        return out
    with p.open() as fh:
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if len(a) < 2:
                continue
            orf = a[0]
            if orf in gff_index:
                c, s, e, _ = gff_index[orf]
                out.append((c, s, e, "ICE:iceberg3"))
    return out


def interval_distance(a_s, a_e, b_s, b_e):
    if a_e < b_s:
        return b_s - a_e
    if b_e < a_s:
        return a_s - b_e
    return 0  # overlap


def tier(d):
    if d <= 2000:
        return "strict"
    if d <= 5000:
        return "standard"
    if d <= 10000:
        return "broad"
    return None


def run_source(src):
    cfg = PATHS[src]
    root = ROOTS[src]
    gff = parse_gff_cds(root / cfg["master_gff"])
    o2c = load_orf2contig(root / cfg["orf2c"])
    args_ = collect_args(root / cfg["amr"], gff, o2c)
    isvs = collect_isescan(root / cfg["ise"])
    ints = collect_integron(root / cfg["intf"])
    ices = collect_iceberg(root / cfg["ice"], gff, o2c)
    mges = isvs + ints + ices

    by_contig_arg = defaultdict(list)
    by_contig_mge = defaultdict(list)
    for c, s, e, lab in args_:
        by_contig_arg[c].append((s, e, lab))
    for c, s, e, lab in mges:
        by_contig_mge[c].append((s, e, lab))

    out_tsv = OUT / f"{src}_arg_mge_distance.tsv"
    n_pairs = 0
    with out_tsv.open("w") as out:
        out.write("source\tcontig\targ_label\targ_start\targ_end\t"
                  "mge_label\tmge_start\tmge_end\tdistance_bp\ttier\n")
        for c, arg_list in by_contig_arg.items():
            mge_list = by_contig_mge.get(c, [])
            if not mge_list:
                continue
            for a_s, a_e, a_lab in arg_list:
                for m_s, m_e, m_lab in mge_list:
                    d = interval_distance(a_s, a_e, m_s, m_e)
                    t = tier(d)
                    if t is None:
                        continue
                    out.write(f"{src}\t{c}\t{a_lab}\t{a_s}\t{a_e}\t"
                              f"{m_lab}\t{m_s}\t{m_e}\t{d}\t{t}\n")
                    n_pairs += 1
    print(f"  {src}: ARG={len(args_)} MGE={len(mges)} pairs<=10kb={n_pairs}")


def main():
    for src in SOURCES:
        if not ROOTS[src].exists():
            print(f"  WARN: {src} root missing"); continue
        run_source(src)
    print("[step4] DONE")


if __name__ == "__main__":
    main()
