#!/usr/bin/env python3
"""
Step 5: Module extraction — for each ARG with an MGE within +-5 kb (uses Step 4
'strict' or 'standard' tier rows), extract the contig region spanning
ARG ± 5 kb flanking. Modules contain at least one MGE.

Input:  $PROJECT/04_cross_track/mobile_arg/step4/{src}_arg_mge_distance.tsv
        per-source contig FASTAs (PL: 01_raw_fasta/all/*.fna ; MAG: master.fna ; UB: master.fna)
Output: $PROJECT/04_cross_track/mobile_arg/step5/{src}_modules.fna
        $PROJECT/04_cross_track/mobile_arg/step5/{src}_modules.tsv  (module_id, contig, start, end, ARGs, MGEs)
"""
import os
import sys
from pathlib import Path
from collections import defaultdict

PROJECT = os.environ.get("PROJECT")
if not PROJECT:
    sys.exit("ERROR: export PROJECT=...")

OUT = Path(PROJECT) / "cross/mobile_arg/step5"
OUT.mkdir(parents=True, exist_ok=True)
STEP4 = Path(PROJECT) / "cross/mobile_arg/step4"

FLANK = 5000
SOURCES = ["plasmid", "mag", "unbinned"]

CONTIG_FASTA = {
    "plasmid":  Path(PROJECT) / "plasmid/04_master_orf/all/master.fna",
    "mag":      Path(PROJECT) / "mag/03_master_orf/all/master.fna",
    "unbinned": Path(PROJECT) / "unbinned/03_master_orf/all/master.fna",
}


def load_fasta(p):
    seqs = {}
    if not p.exists():
        return seqs
    name, chunks = None, []
    with p.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs


def run_source(src):
    dist_tsv = STEP4 / f"{src}_arg_mge_distance.tsv"
    if not dist_tsv.exists():
        print(f"  skip {src}: no step4 output"); return
    seqs = load_fasta(CONTIG_FASTA[src])
    if not seqs:
        print(f"  WARN {src}: no contig FASTA at {CONTIG_FASTA[src]}"); return

    # Group: per contig collect ARGs (with coords) and MGEs (with coords) tier<=standard
    arg_per_contig = defaultdict(list)   # (s,e,label)
    mge_per_contig = defaultdict(list)
    seen_arg = set(); seen_mge = set()
    with dist_tsv.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if a[ix["tier"]] not in ("strict", "standard"):
                continue
            c = a[ix["contig"]]
            arg_s = int(a[ix["arg_start"]]); arg_e = int(a[ix["arg_end"]])
            arg_l = a[ix["arg_label"]]
            mge_s = int(a[ix["mge_start"]]); mge_e = int(a[ix["mge_end"]])
            mge_l = a[ix["mge_label"]]
            k1 = (c, arg_s, arg_e, arg_l)
            if k1 not in seen_arg:
                arg_per_contig[c].append((arg_s, arg_e, arg_l))
                seen_arg.add(k1)
            k2 = (c, mge_s, mge_e, mge_l)
            if k2 not in seen_mge:
                mge_per_contig[c].append((mge_s, mge_e, mge_l))
                seen_mge.add(k2)

    out_fna = OUT / f"{src}_modules.fna"
    out_tsv = OUT / f"{src}_modules.tsv"
    n_mod = 0
    with out_fna.open("w") as ofa, out_tsv.open("w") as otsv:
        otsv.write("module_id\tsource\tcontig\tstart\tend\tlength\t"
                   "ARG_labels\tMGE_labels\n")
        for c, args_ in arg_per_contig.items():
            seq = seqs.get(c)
            if seq is None:
                continue
            L = len(seq)
            for (a_s, a_e, a_l) in args_:
                ms = max(1, a_s - FLANK)
                me = min(L, a_e + FLANK)
                # MGEs that fall inside [ms, me]
                inside_mge = [m for m in mge_per_contig.get(c, [])
                              if m[1] >= ms and m[0] <= me]
                if not inside_mge:
                    continue
                mod_id = f"{src}|{c}|{ms}_{me}|{a_l}"
                n_mod += 1
                ofa.write(f">{mod_id}\n{seq[ms-1:me]}\n")
                otsv.write(f"{mod_id}\t{src}\t{c}\t{ms}\t{me}\t{me-ms+1}\t"
                           f"{a_l}\t{','.join(m[2] for m in inside_mge)}\n")
    print(f"  {src}: modules={n_mod} -> {out_fna}")


def main():
    for src in SOURCES:
        run_source(src)
    print("[step5] DONE")


if __name__ == "__main__":
    main()
