#!/usr/bin/env python3
"""
Step 7: Mobility pathway classification.

Reads step6 filtered BLAST hits and assigns each (module, hit) to a pathway:
  plasmid <-> plasmid  (further: same-pOTU vs cross-pOTU using PL pOTU map if available)
  plasmid <-> mag
  plasmid <-> unbinned
  mag <-> mag
  mag <-> unbinned
  unbinned <-> unbinned

Same-origin-region hits (already excluded in step6 by self-contig) are skipped.

Source of module = first field of "src|contig|start_end|ARG:..."
Source of target = prefix of sseqid (PL_/MAG_/UB_)

Optional pOTU map: $PROJECT/plasmid/30_clustering/potu_map.tsv  (contig<TAB>pOTU_id)

Output: $PROJECT/cross/mobile_arg/step7/mobility_pathways.tsv
"""
import os
import sys
from pathlib import Path

PROJECT = os.environ.get("PROJECT")
if not PROJECT:
    sys.exit("ERROR: export PROJECT=...")

OUT = Path(PROJECT) / "cross/mobile_arg/step7"
OUT.mkdir(parents=True, exist_ok=True)

FILT = Path(PROJECT) / "cross/mobile_arg/step6/modules_vs_combined.filt.tsv"
if not FILT.exists():
    sys.exit(f"ERROR: missing {FILT} — run step6 first")

POTU_MAP = Path(PROJECT) / "plasmid/30_clustering/potu_map.tsv"
potu = {}
if POTU_MAP.exists():
    with POTU_MAP.open() as fh:
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if len(a) >= 2:
                potu[a[0]] = a[1]

PFX2SRC = {"PL": "plasmid", "MAG": "mag", "UB": "unbinned"}

out_tsv = OUT / "mobility_pathways.tsv"
n_keep = 0
with FILT.open() as fh, out_tsv.open("w") as out:
    header = fh.readline().rstrip("\n").split("\t")
    out.write("\t".join(header) + "\tquery_source\ttarget_source\t"
              "pathway\tquery_contig\ttarget_contig\tpotu_relation\n")
    for line in fh:
        a = line.rstrip("\n").split("\t")
        d = dict(zip(header, a))
        qseqid = d["qseqid"]; sseqid = d["sseqid"]
        # query module id: "src|contig|start_end|ARG:..."
        qparts = qseqid.split("|")
        if len(qparts) < 4:
            continue
        q_src = qparts[0]; q_contig = qparts[1]
        # target sseqid: "PFX_contig"
        pos = sseqid.find("_")
        if pos <= 0:
            continue
        t_pfx = sseqid[:pos]
        t_contig = sseqid[pos+1:]
        t_src = PFX2SRC.get(t_pfx)
        if t_src is None:
            continue
        if q_contig == t_contig:
            continue  # belt-and-suspenders self-hit filter
        # Pathway label canonicalized alphabetically
        endpoints = tuple(sorted([q_src, t_src]))
        pathway = f"{endpoints[0]}<->{endpoints[1]}"
        # pOTU relation only meaningful when both ends are plasmid
        if q_src == "plasmid" and t_src == "plasmid" and potu:
            qp = potu.get(q_contig, "NA"); tp = potu.get(t_contig, "NA")
            if qp == "NA" or tp == "NA":
                rel = "potu_unknown"
            elif qp == tp:
                rel = "same_pOTU"
            else:
                rel = "cross_pOTU"
        else:
            rel = "NA"
        out.write("\t".join(a) + f"\t{q_src}\t{t_src}\t{pathway}\t"
                  f"{q_contig}\t{t_contig}\t{rel}\n")
        n_keep += 1

# Quick pathway counts
counts = {}
with out_tsv.open() as fh:
    h = fh.readline().rstrip("\n").split("\t")
    i_pw = h.index("pathway")
    for line in fh:
        p = line.rstrip("\n").split("\t")[i_pw]
        counts[p] = counts.get(p, 0) + 1

print(f"[step7] DONE — classified {n_keep} hits -> {out_tsv}")
for k, v in sorted(counts.items(), key=lambda kv: -kv[1]):
    print(f"  {k}: {v}")
