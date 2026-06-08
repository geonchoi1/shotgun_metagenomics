#!/usr/bin/env python
"""Build & verify ARG-MGE cassette-sharing pools.
 (1) all_candidates.tsv : every module(ARG±5kb)-vs-contig hit (>=2kb, >=95% id, non-self) = candidate pool.
 (2) For each distinct (query_plasmid, target) pair, blast the WHOLE query plasmid vs the target and measure
     query coverage -> verifies query and target are DIFFERENT replicons (low cov = cassette-only) vs same/near-same.
 (3) cassette_only.tsv : verified pairs where the two replicons share ONLY the cassette (whole-plasmid cov < 30%).
Source raw blast: /tmp/xsource/{pl,mag,ub}.blast.tsv (modules vs PL/MAG/UB, pident>=90)."""
import subprocess,tempfile,os,collections
from pyfaidx import Fasta
W="/home/gchoi/wwtp_plasmidome"; OUT=f"{W}/analysis/04_arg_mge/cassette_sharing"
pl=Fasta(f"{W}/01_plasmid_track/01_raw_fasta/dereplicated_1869.fna")
ub=Fasta(f"{W}/03_unbinned_track/01_raw_fasta/all.fna")
def seq(c): return str(pl[c][:].seq) if c in pl else (str(ub[c][:].seq) if c in ub else None)
def clen(c): return len(pl[c]) if c in pl else (len(ub[c]) if c in ub else 0)
tier={}
for i,ln in enumerate(open(f"{W}/analysis/plasmid_cargo/plasmid_cargo.tsv")):
    p=ln.rstrip("\n").split("\t")
    if i==0: ci=p.index("contig"); ti=p.index("tier"); continue
    tier[p[ci]]=p[ti]
# strict ARG: AMRFinder Type=AMR symbols only (exclude STRESS metal/biocide e.g. mer,qacEdelta1)
AMRG=set()
for i,ln in enumerate(open(f"{W}/01_plasmid_track/09_amrfinder/amr.tsv")):
    p=ln.rstrip("\n").split("\t")
    if i==0: HH={k:j for j,k in enumerate(p)}; continue
    if p[HH["Type"]]=="AMR": AMRG.add(p[HH["Element symbol"]])
# (1) candidate pool from raw blast
SRCFILE={"PL":"pl","MAG":"mag","UB":"ub"}
pool=[]  # (src,qcontig,qtier,arg,target,pident,alen)
pair=collections.defaultdict(lambda:{"alen":0,"args":set(),"srcs":set()})  # (q,t)->agg
for src,tag in SRCFILE.items():
    f=f"/tmp/xsource/{tag}.blast.tsv"
    if not os.path.exists(f): continue
    for ln in open(f):
        c=ln.rstrip("\n").split("\t")
        qid,sid,pid,alen=c[0],c[1],float(c[2]),int(c[3])
        a=qid.split("|"); qc=a[0]+"|"+a[1]; arg=a[3]
        if arg not in AMRG: continue                          # strict ARG only (drop metal/biocide)
        if pid<95 or alen<2000 or sid==qc: continue
        pool.append((src,qc,tier.get(qc,""),arg,sid,pid,alen))
        k=(qc,sid); pair[k]["alen"]=max(pair[k]["alen"],alen); pair[k]["args"].add(arg); pair[k]["srcs"].add(src)
# (2) whole-plasmid coverage per query (verify different replicon)
q2t=collections.defaultdict(set)
for (qc,sid) in pair: q2t[qc].add(sid)
cov={}  # (q,t)->coverage fraction of query
for q in q2t:
    qs=seq(q)
    if not qs: continue
    ts_list=[t for t in q2t[q] if seq(t)]
    tf=tempfile.NamedTemporaryFile("w",suffix=".fa",delete=False)
    for t in ts_list: tf.write(f">{t.replace('|','_')}\n{seq(t)}\n")
    tf.close(); qf=tempfile.NamedTemporaryFile("w",suffix=".fa",delete=False); qf.write(f">q\n{qs}\n"); qf.close()
    out=subprocess.run(["blastn","-query",qf.name,"-subject",tf.name,"-perc_identity","95",
        "-outfmt","6 sseqid qstart qend length"],capture_output=True,text=True).stdout
    cv=collections.defaultdict(lambda:[0]*(len(qs)+1))
    for l in out.splitlines():
        s,a,b,L=l.split("\t")
        if int(L)<500: continue
        lo,hi=sorted((int(a),int(b)))
        for i in range(lo,hi+1): cv[s][i]=1
    os.unlink(tf.name); os.unlink(qf.name)
    for t in ts_list:
        sk=t.replace("|","_"); cov[(q,t)]=sum(cv[sk])/len(qs) if sk in cv else 0
def ttype(t): return "plasmid" if t in pl else "chromosomal(unbinned)"
def verdict(c): return "cassette-only" if c<0.30 else ("partial" if c<0.80 else "near-whole/same")
# (3) write
with open(f"{OUT}/all_candidates.tsv","w") as fo:
    fo.write("src\tquery_contig\tquery_tier\tquery_len\tARG\ttarget\ttarget_type\ttarget_len\tcassette_align_bp\tpident\twhole_plasmid_cov_pct\tverdict\n")
    for src,qc,qt,arg,t,pid,alen in sorted(pool):
        c=cov.get((qc,t),0)
        fo.write(f"{src}\t{qc}\t{qt}\t{clen(qc)}\t{arg}\t{t}\t{ttype(t)}\t{clen(t)}\t{alen}\t{pid:.1f}\t{c*100:.0f}\t{verdict(c)}\n")
with open(f"{OUT}/cassette_only.tsv","w") as fo:
    fo.write("query_contig\tquery_tier\tquery_len\tARG\ttarget\ttarget_type\ttarget_len\tcassette_align_bp\twhole_plasmid_cov_pct\n")
    seen=set()
    for (qc,t),pa in pair.items():
        c=cov.get((qc,t),1)
        if c<0.30 and (qc,t) not in seen:
            seen.add((qc,t))
            fo.write(f"{qc}\t{tier.get(qc,'')}\t{clen(qc)}\t{','.join(sorted(pa['args']))}\t{t}\t{ttype(t)}\t{clen(t)}\t{pa['alen']}\t{c*100:.0f}\n")
# summary
allp=len(pair); cof=sum(1 for k in pair if cov.get(k,1)<0.30)
print(f"candidate pool rows (module hits)={len(pool)}")
print(f"distinct (query,target) pairs={allp}")
print(f"cassette-only verified pairs (cov<30%)={cof}")
print(f"-> {OUT}/all_candidates.tsv , {OUT}/cassette_only.tsv")
