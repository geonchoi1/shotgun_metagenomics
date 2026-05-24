#!/usr/bin/env python3
"""
Step 8: ARG mobility network — Cytoscape/Gephi compatible nodes + edges.

Node = contig (annotated with source PL/MAG/UB and optional pOTU/taxonomy).
Edge = mobility link from step7 (one edge per (query_contig, target_contig, ARG, tier)).

Outputs (under $PROJECT/cross/mobile_arg/step8/):
  nodes.tsv:  id, source, n_args, attrs(potu, taxonomy)
  edges.tsv:  source, target, ARG, identity_tier, pathway, potu_relation, pident, length
"""
import os
import sys
from pathlib import Path
from collections import defaultdict

PROJECT = os.environ.get("PROJECT")
if not PROJECT:
    sys.exit("ERROR: export PROJECT=...")

OUT = Path(PROJECT) / "cross/mobile_arg/step8"
OUT.mkdir(parents=True, exist_ok=True)

PATH7 = Path(PROJECT) / "cross/mobile_arg/step7/mobility_pathways.tsv"
if not PATH7.exists():
    sys.exit(f"ERROR: missing {PATH7} — run step7 first")

# Optional metadata maps
POTU = Path(PROJECT) / "plasmid/30_clustering/potu_map.tsv"
UB_TAX = Path(PROJECT) / "unbinned/22_mmseqs2_lca_taxonomy/all/contig_taxonomy.tsv"
MAG_TAX = Path(PROJECT) / "mag/10_gtdbtk/gtdbtk.summary.tsv"

potu = {}
if POTU.exists():
    for line in POTU.open():
        a = line.rstrip("\n").split("\t")
        if len(a) >= 2:
            potu[a[0]] = a[1]

ub_tax = {}
if UB_TAX.exists():
    with UB_TAX.open() as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            a = line.rstrip("\n").split("\t")
            if len(a) >= 8:
                # join class..species short
                ub_tax[a[0]] = f"d_{a[1]};p_{a[2]};c_{a[3]};g_{a[6]}"

mag_tax = {}
if MAG_TAX.exists():
    try:
        with MAG_TAX.open() as fh:
            head = fh.readline().rstrip("\n").split("\t")
            i_id = head.index("user_genome") if "user_genome" in head else 0
            i_cl = head.index("classification") if "classification" in head else -1
            for line in fh:
                a = line.rstrip("\n").split("\t")
                if i_cl >= 0 and len(a) > max(i_id, i_cl):
                    mag_tax[a[i_id]] = a[i_cl]
    except Exception:
        pass

# Aggregate edges
edges = []     # list of dicts
node_arg = defaultdict(set)   # contig -> set of ARG labels
node_src = {}                  # contig -> source

with PATH7.open() as fh:
    head = fh.readline().rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(head)}
    for line in fh:
        a = line.rstrip("\n").split("\t")
        q_contig = a[ix["query_contig"]]
        t_contig = a[ix["target_contig"]]
        q_src = a[ix["query_source"]]
        t_src = a[ix["target_source"]]
        pident = a[ix["pident"]]
        length = a[ix["length"]]
        tier   = a[ix["tier"]]
        pw     = a[ix["pathway"]]
        rel    = a[ix["potu_relation"]]
        # ARG label is embedded in qseqid module id: "src|contig|s_e|ARG:gene"
        qseqid = a[ix["qseqid"]]
        qp = qseqid.split("|")
        arg = qp[3] if len(qp) >= 4 else "ARG:?"
        edges.append(dict(source=q_contig, target=t_contig, ARG=arg,
                          identity_tier=tier, pathway=pw,
                          potu_relation=rel, pident=pident, length=length))
        node_src[q_contig] = q_src
        node_src[t_contig] = t_src
        node_arg[q_contig].add(arg)
        node_arg[t_contig].add(arg)

# Write nodes.tsv
nodes_tsv = OUT / "nodes.tsv"
with nodes_tsv.open("w") as out:
    out.write("id\tsource\tn_args\tARGs\tpOTU\ttaxonomy\n")
    for c in sorted(node_src):
        src = node_src[c]
        args_str = ",".join(sorted(node_arg.get(c, set())))
        p = potu.get(c, "NA")
        if src == "unbinned":
            tax = ub_tax.get(c, "NA")
        elif src == "mag":
            tax = mag_tax.get(c, "NA")
        else:
            tax = "NA"
        out.write(f"{c}\t{src}\t{len(node_arg.get(c, set()))}\t{args_str}\t{p}\t{tax}\n")

# Write edges.tsv (dedupe identical 5-tuples keeping max pident)
edges_tsv = OUT / "edges.tsv"
key2best = {}
for e in edges:
    k = (e["source"], e["target"], e["ARG"], e["identity_tier"])
    cur = key2best.get(k)
    if cur is None or float(e["pident"]) > float(cur["pident"]):
        key2best[k] = e
with edges_tsv.open("w") as out:
    out.write("source\ttarget\tARG\tidentity_tier\tpathway\tpotu_relation\tpident\tlength\n")
    for k, e in sorted(key2best.items()):
        out.write(f"{e['source']}\t{e['target']}\t{e['ARG']}\t"
                  f"{e['identity_tier']}\t{e['pathway']}\t"
                  f"{e['potu_relation']}\t{e['pident']}\t{e['length']}\n")

print(f"[step8] DONE — nodes={len(node_src)} edges={len(key2best)}")
print(f"  {nodes_tsv}")
print(f"  {edges_tsv}")
