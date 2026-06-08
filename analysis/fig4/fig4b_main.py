#!/usr/bin/env python
"""Fig4b (main) — three shared ARG-MGE cassettes across mobility classes, clinical diversity:
  A pdCONJ : multi-ARG ICE  erm(F)+tet(X1/X2)+tet(Q)   plasmid Anaerobic|contig_37029 ~ CHROMOSOMAL IN|contig_1562
  D pMOB   : class-1 integron w/ carbapenemase blaGES-5 (+aac(6')-IIa,blaOXA-17,sul1)  plasmid IN|contig_7336 ~ plasmid IN|contig_8285
  F pCONJ  : sul2 + IS91     plasmid Oxic|contig_16581 ~ plasmid RAS|contig_34541 (cross-zone)
Plasmid genes annotated from master.gff + AMRFinder(ARG/Metal) + ISEScan(IS family) + ICEberg3(ICE);
chromosomal copy from bakta .gbff. Same-named resistance genes linked between copies; tier badge per plasmid."""
import re, csv, numpy as np
from collections import defaultdict, Counter
from Bio import SeqIO
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Polygon, FancyBboxPatch, Patch
from matplotlib.lines import Line2D
W="/home/gchoi/wwtp_plasmidome"; OUT=f"{W}/analysis/04_arg_mge"; PT=f"{W}/01_plasmid_track"
PGFF=f"{PT}/04_master_orf/plasmidome.master.gff"; AMR=f"{PT}/09_amrfinder/amr.tsv"
COL={"ARG":"#d62728","Metal":"#8c564b","IS":"#2166ac","ICE":"#1b7837","Integron":"#ff7f0e","other":"#cfcfcf"}
TIERCOL={"pdCONJ":"#6a51a3","pMOB":"#2171b5","pCONJ":"#238b45","pNT":"#737373","pOriT":"#d94801"}
ARGRE=re.compile(r"^(tet|erm|bla|sul|aac|aph|mcr|qac|cat|dfr|mph|aad|ant|van|qnr|flo|cml|ble|fos|ges|ere|mef|msr|cfx)",re.I)
def fam(n): return re.sub(r"\d+$","",n)   # cfxA3 -> cfxA (ARG variant family for homology link)
METALRE=re.compile(r"^(mer|ars|cop|czc|pco|sil|cad|znt)",re.I)
def shortlab(lab,c):
    if c in("ARG","Metal"): return lab
    l=lab.lower()
    for k,t in [("rtec","rteC"),("inti","intI1"),("transpos","transposase"),("integrase","integrase"),
                ("recombinase","recombinase"),("relaxase","relaxase"),("conjug","conjugation"),
                ("mobiliz","mob relaxase"),("excision","excisionase"),("helicase","helicase")]:
        if k in l: return t
    return lab if len(lab)<=13 else lab[:12]+"…"
def catUB(label):
    if not label: return "other"
    if len(label)<=14 and METALRE.search(label): return "Metal"
    if len(label)<=14 and ARGRE.search(label): return "ARG"
    if re.search(r"intI|integron integrase",label,re.I): return "Integron"
    if re.search(r"transpos|insertion sequence|\bIS[0-9]|\btnp\b",label,re.I): return "IS"
    if re.search(r"integrase|recombinase|conjug|relaxase|mobiliz|excision|\brte[abc]\b|\btra[a-z]\b|pilus",label,re.I): return "ICE"
    return "other"
# AMRFinder ARG/metal locus -> (symbol, ARG|Metal)
amr={}
for i,ln in enumerate(open(AMR)):
    p=ln.rstrip("\n").split("\t")
    if i==0: H={k:j for j,k in enumerate(p)}; continue
    if p[H["Type"]]=="AMR":    amr[p[H["Protein id"]]]=(p[H["Element symbol"]],"ARG")   # AMR only (metal dropped per request)
# plasmid contig -> mobility tier (5-class)
PTIER={}
for i,ln in enumerate(open(f"{W}/analysis/plasmid_cargo/plasmid_cargo.tsv")):
    p=ln.rstrip("\n").split("\t")
    if i==0: _ci=p.index("contig"); _ti=p.index("tier"); continue
    PTIER[p[_ci]]=p[_ti]
def plabel(contig,is_chrom):
    if is_chrom: return f"Chromosome: {contig}"
    t=PTIER.get(contig,""); return f"Plasmid ({t}): {contig}" if t else f"Plasmid: {contig}"
# ISEScan IS regions + ICEberg ICE regions per contig
ISREG=defaultdict(list)
for r in csv.DictReader(open(f"{PT}/20_isescan/is_summary.csv")):
    fam_lbl=r["family"] if r["family"].startswith("IS") else "IS"+r["family"]; ISREG[r["seqID"]].append((int(r["isBegin"]),int(r["isEnd"]),fam_lbl))
ICEREG=defaultdict(list)
try:
    for ln in open(f"{PT}/22_iceberg3/iceberg3_filtered.tsv"):
        f=ln.split("\t");
        if len(f)>8: ICEREG[f[0]].append((int(f[6]),int(f[7])))
except FileNotFoundError: pass
# class-1 integron intI1 + attC sites (IntegronFinder --local-max re-run on cassette contigs)
INTI={"IN|contig_7336":[(996,2009)],"IN|contig_8285":[(23641,24654)]}
ATTC={"IN|contig_7336":[(3062,3171),(3740,3799)],"IN|contig_8285":[(21851,21910),(22479,22588)]}
def ov(a0,a1,b0,b1): return max(0,min(a1,b1)-max(a0,b0))
def genes_pl(contig,reg,use_ice=False):
    s0,s1=reg; g=[]
    for ln in open(PGFF):
        if "\tCDS\t" not in ln: continue
        f=ln.split("\t")
        if f[0]!=contig: continue
        st,en=int(f[3]),int(f[4])
        if en<s0 or st>s1: continue
        lt=re.search(r"ID=([^;]+)",f[8]); lt=lt.group(1) if lt else ""
        c,lab="other",""
        if lt in amr: lab,c=amr[lt]
        else:
            isf=[fam for (a,b,fam) in ISREG.get(contig,[]) if ov(st,en,a,b)>0.3*(en-st)]
            if isf: c,lab="IS",isf[0]                        # specific IS transposase ORF (ISEScan)
            elif any(ov(st,en,a,b)>0.5*(en-st) for (a,b) in INTI.get(contig,[])): c,lab="Integron","intI1"  # integron integrase
        g.append([st,en,1 if f[6]=="+" else -1,lab,c])
    return g
def genes_ub(gbff,contig,reg):
    s0,s1=reg; g=[]
    for r in SeqIO.parse(gbff,"genbank"):
        if r.id!=contig: continue
        for ft in r.features:
            if ft.type!="CDS": continue
            st=int(ft.location.start)+1; en=int(ft.location.end)
            if en<s0 or st>s1: continue
            lab=ft.qualifiers.get("gene",ft.qualifiers.get("product",[""]))[0]
            g.append([st,en,ft.location.strand or 1,lab,catUB(lab)])
        break
    return g
def norm(g,reg,rev):
    s0,s1=reg; L=s1-s0; out=[]
    for st,en,sd,lab,c in g:
        if not rev: x0,x1,s=st-s0,en-s0,sd
        else:       x0,x1,s=s1-en,s1-st,-sd
        out.append([max(0,x0),min(L,x1),s,lab,c])
    return out,L
CHK=f"{W}/03_unbinned_track/02_bakta/linear_chunks"
WIN=5000   # +-5 kb around the shared cassette anchor (10 kb total, fixed for every row)
CASS=[  # verified cassette-only: plasmid <-> chromosome (low genome-wide id) sharing ONLY the cassette
 dict(tier="pMOB",title="aph(3')-Ia composite transposon (IS1182)",rid="17% genome-wide id",note=None,anchor="aph(3')-Ia",
      top=("Anaerobic|contig_36490",32624,None), bot=("IN|contig_92",4491,f"{CHK}/out_chunk_03/chunk_03.gbff")),
 dict(tier="pMOB",title="mef(A)/msr(D) macrolide efflux + IS5",rid="10% genome-wide id",note=None,anchor="mef(A)",
      top=("Anaerobic|contig_17698",6351,None), bot=("IN|contig_15622",6993,f"{CHK}/out_chunk_00/chunk_00.gbff")),
 dict(tier="pCONJ",title="sul2 + IS91",rid="5% genome-wide id",note=None,anchor="sul2",
      top=("RAS|contig_34541",28233,None), bot=("Oxic|contig_16151",10111,f"{CHK}/out_chunk_08/chunk_08.gbff")),
]
def fetch_anchored(side,anchor):
    contig,apos,gbff=side
    genes=genes_ub(gbff,contig,(max(1,apos-12000),apos+12000)) if gbff else genes_pl(contig,(max(1,apos-12000),apos+12000),False)
    ag=[g for g in genes if g[4]in("ARG","Metal") and fam(g[3])==fam(anchor)]
    if not ag: return [],contig,apos,1,bool(gbff)
    a=min(ag,key=lambda g:abs((g[0]+g[1])/2-apos)); ac=(a[0]+a[1])/2; asd=a[2] or 1
    out=[]
    for st,en,sd,lab,c in genes:
        x0,x1,s=(st-ac,en-ac,(sd or 1)) if asd>0 else (ac-en,ac-st,-(sd or 1))
        if x1<-WIN or x0>WIN: continue
        out.append([max(-WIN,x0),min(WIN,x1),s,lab,c])
    return out,contig,ac,asd,bool(gbff)
rows=[]
for cs in CASS:
    tg,tcontig,tac,tasd,tchrom=fetch_anchored(cs["top"],cs["anchor"])
    bg,bcontig,bac,basd,bchrom=fetch_anchored(cs["bot"],cs["anchor"])
    if bchrom:                                               # bakta target -> transfer IS/ICE/Integron to prodigal plasmid top
        for gt in tg:
            if gt[4]=="other":
                best=max(bg,key=lambda u:ov(gt[0],gt[1],u[0],u[1]),default=None)
                if best and ov(gt[0],gt[1],best[0],best[1])>0.4*(gt[1]-gt[0]) and best[4] in("IS","ICE","Integron"):
                    gt[4]=best[4]; gt[3]=best[3]
    rows.append((cs,tg,bg,(tcontig,tac,tasd),(bcontig,bac,basd)))
    print(cs["tier"],cs["title"],"| top",[g[3] for g in tg if g[4]in("ARG","Metal")],"| bot",[g[3] for g in bg if g[4]in("ARG","Metal")])

from pyfaidx import Fasta as _Fa
import subprocess as _sp, tempfile as _tf, os as _os
_PLF=_Fa(f"{W}/01_plasmid_track/01_raw_fasta/dereplicated_1869.fna"); _UBF=_Fa(f"{W}/03_unbinned_track/01_raw_fasta/all.fna")
def clen(c): return len(_PLF[c]) if c in _PLF else (len(_UBF[c]) if c in _UBF else 10**9)
def _seqwin(c,center):
    s=_PLF[c] if c in _PLF else _UBF[c]; ws=max(1,int(center)-12000); we=min(len(s),int(center)+12000)
    return str(s[ws-1:we].seq),ws
def hsp_bands(tc,tac,tasd,bc,bac,basd):   # true sequence-level homology (incl non-coding) between the two copies
    qseq,qws=_seqwin(tc,tac); sseq,sws=_seqwin(bc,bac)
    qf=_tf.NamedTemporaryFile("w",suffix=".fa",delete=False); qf.write(">q\n"+qseq+"\n"); qf.close()
    sf=_tf.NamedTemporaryFile("w",suffix=".fa",delete=False); sf.write(">s\n"+sseq+"\n"); sf.close()
    out=_sp.run(["blastn","-query",qf.name,"-subject",sf.name,"-perc_identity","90","-outfmt","6 qstart qend sstart send length"],capture_output=True,text=True).stdout
    _os.unlink(qf.name); _os.unlink(sf.name)
    tx=lambda p:(p-tac) if tasd>0 else (tac-p); bx=lambda p:(p-bac) if basd>0 else (bac-p)
    bands=[]
    for l in out.splitlines():
        qa,qb,sa,sb,L=map(int,l.split("\t"))
        if L<300: continue
        qx0,qx1=sorted((tx(qws+qa-1),tx(qws+qb-1))); sx0,sx1=sorted((bx(sws+sa-1),bx(sws+sb-1)))
        if qx1<-WIN or qx0>WIN: continue
        if abs((qx0+qx1)/2-(sx0+sx1)/2)>2500: continue       # drop off-diagonal repeat matches (keep cassette block)
        bands.append((max(-WIN,qx0),min(WIN,qx1),max(-WIN,sx0),min(WIN,sx1)))
    return bands
H=0.34
fig,axes=plt.subplots(len(CASS),1,figsize=(10.5,7.8),gridspec_kw=dict(hspace=0.45))
def arrow(ax,x0,x1,y,sd,c):
    w=x1-x0
    if w<=0: return
    hl=min(330,w*0.5)
    xx,ww=(x0,w) if sd>=0 else (x1,-w)
    ax.add_patch(FancyArrow(xx,y,ww,0,width=H*1.3,head_width=H*1.85,head_length=hl,length_includes_head=True,facecolor=COL[c],edgecolor="black",lw=0.5,zorder=4))
for ax,(cs,tg,bg,tanc,banc) in zip(axes,rows):
    yT,yB=0.8,-0.8
    for qx0,qx1,sx0,sx1 in hsp_bands(*tanc,*banc):   # sequence-homology shading (coding + non-coding)
        ax.add_patch(Polygon([(qx0,yT-H*0.7),(qx1,yT-H*0.7),(sx1,yB+H*0.7),(sx0,yB+H*0.7)],
                     closed=True,facecolor="#cf9090",alpha=0.32,edgecolor="none",zorder=1))
    for (gl,y,anc) in [(tg,yT,tanc),(bg,yB,banc)]:
        _c,_ac,_asd=anc; _L=clen(_c)                          # clip dashed baseline to the actual contig extent
        _xs,_xe=(1-_ac,_L-_ac) if _asd>0 else (_ac-_L,_ac-1)
        ax.plot([max(-WIN,_xs),min(WIN,_xe)],[y,y],color="0.55",lw=1.0,ls=(0,(5,3)),zorder=2)
        for x0,x1,sd,lab,c in gl: arrow(ax,x0,x1,y,sd,c)
        for x0,x1,sd,lab,c in gl:
            if c=="other": continue
            txt=lab if c in("ARG","Metal") else shortlab(lab,c)
            if not txt: continue
            ax.text((x0+x1)/2,y+(H*1.6 if y>0 else -H*1.6),txt,ha="left" if y>0 else "right",
                    va="bottom" if y>0 else "top",fontsize=7.3 if c in("ARG","Metal") else 6.2,
                    fontweight="bold",color=COL[c],rotation=36,rotation_mode="anchor",zorder=6)
    # attC sites (re-anchored)
    for (contig,ac,asd,y) in [(tanc[0],tanc[1],tanc[2],yT),(banc[0],banc[1],banc[2],yB)]:
        for (a,b) in ATTC.get(contig,[]):
            mid=(a+b)/2; x=(mid-ac) if asd>0 else (ac-mid)
            if -WIN<=x<=WIN: ax.scatter(x,y,marker="v",s=46,facecolor=COL["Integron"],edgecolor="black",lw=0.6,zorder=8)
    # integron / ICE bracket above the genes
    if cs["note"]:
        ntyp,nlbl=cs["note"]; yb=yT+H*5.3
        if ntyp=="ICE": x0n,x1n=-WIN,WIN
        else:
            xs=[g[0] for g in tg if g[4]in("ARG","Integron")]+[g[1] for g in tg if g[4]in("ARG","Integron")]
            x0n,x1n=(min(xs),max(xs)) if xs else (-1000,1000)
        ax.plot([x0n,x1n],[yb,yb],color=COL[ntyp],lw=2.4,zorder=5)
        ax.plot([x0n,x0n],[yb,yb-H*0.5],color=COL[ntyp],lw=2.4,zorder=5); ax.plot([x1n,x1n],[yb,yb-H*0.5],color=COL[ntyp],lw=2.4,zorder=5)
        ax.text((x0n+x1n)/2,yb+H*0.25,nlbl,ha="center",va="bottom",fontsize=8.5,fontweight="bold",color=COL[ntyp])
    # left contig labels (Plasmid(tier): / Chromosome:)
    chrom=bool(cs["bot"][2])
    ax.text(-WIN*1.05,yT,plabel(cs['top'][0],False),ha="right",va="center",fontsize=7.5,fontweight="bold")
    ax.text(-WIN*1.05,yB,plabel(cs['bot'][0],chrom),ha="right",va="center",fontsize=7.5,fontweight="bold",color="0.35")
    ax.set_xlim(-WIN*1.07,WIN*1.10); ax.set_ylim(-1.7,3.2); ax.axis("off")
axes[-1].set_ylim(-2.6,3.2)   # room for scale bar (bottom-right)
_x1=WIN*1.0; _x0=_x1-2000
axes[-1].plot([_x0,_x1],[-2.2,-2.2],color="black",lw=2.5); axes[-1].text((_x0+_x1)/2,-2.4,"2 kb",ha="center",va="top",fontsize=8.5,fontweight="bold")
_hand=[Patch(facecolor=COL[k],edgecolor="black",label=k) for k in ["ARG","IS","ICE","Integron","other"]]
_hand.append(Line2D([0],[0],marker="v",color="w",markerfacecolor=COL["Integron"],markeredgecolor="black",markersize=8,label="attC"))
lg=axes[0].legend(handles=_hand,loc="lower center",bbox_to_anchor=(0.5,1.15),ncol=6,frameon=False,fontsize=9.5,handletextpad=0.4,columnspacing=1.1)
plt.setp(lg.get_texts(),fontweight="bold")
# save contig/label info to TSV (user types these into the PPT manually)
with open(f"{OUT}/Fig4b_cassette_labels.tsv","w") as fo:
    fo.write("description\ttier\ttop_label\ttop_contig\tbot_label\tbot_contig\tanchor_ARG\tgenome_wide_id\twindow\tARGs\n")
    for cs,tg,bg,tanc,banc in rows:
        args=",".join(sorted({g[3] for g in tg if g[4] in("ARG","Metal")}))
        chrom=bool(cs["bot"][2])
        fo.write(f"{cs['title']}\t{cs['tier']}\tPlasmid ({cs['tier']})\t{cs['top'][0]}\t{'Chromosome' if chrom else 'Plasmid ('+cs['tier']+')'}\t{cs['bot'][0]}\t{cs['anchor']}\t{cs['rid']}\t±10 kb (20 kb)\t{args}\n")
fig.savefig(f"{OUT}/Fig4b_synteny.png",dpi=600,bbox_inches="tight"); plt.close(fig)
from PIL import Image; im=Image.open(f"{OUT}/Fig4b_synteny.png"); w,h=im.size; im.resize((1100,int(1100*h/w))).save("/tmp/fig4b_main.png")
print("done")
