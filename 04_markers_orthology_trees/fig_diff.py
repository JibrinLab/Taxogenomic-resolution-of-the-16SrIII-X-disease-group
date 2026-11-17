# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
import os,re,sys,numpy as np,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
matplotlib.rcParams.update({"font.family":"DejaVu Sans","pdf.fonttype":42,"svg.fonttype":"none"})
B=os.environ.get("ANI_BASE","/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy")
OUT=B+"/15_figures/diff"; os.makedirs(OUT,exist_ok=True)
PPS=["PPS1A","PPS2B","PPS6A","PPS6B","PPS9B","PPS10A"]
G2=["CanS4","MW1","VAC","MA1"]; G1=["CX","6A1","ZT3_1","PR2021","JR1","PDA15"]
OT=["LWB","CicWB_2022","ChTDIII","OVSP_L01","OVSP_S01"]; UN=["ChTDIII","OVSP_L01","OVSP_S01"]
ALL=PPS+G2+G1+OT
GC={"G2":"#C2410C","G1":"#1D4E89","UN":"#2A8C6E"}
def fa(p):
    d={};k=None
    for l in open(p):
        l=l.rstrip()
        if l.startswith(">"): k=l[1:].split()[0]; d[k]=[]
        elif k: d[k].append(l)
    return {a:"".join(b) for a,b in d.items()}
def lab(name):
    n=name.replace("-","_")
    if "PBT" in n: return "PBT"
    for l in ALL:
        if re.search(r"(?<![A-Za-z0-9])"+re.escape(l)+r"(?![A-Za-z0-9])",n): return l
    return None
def col(l): return GC["G2"] if l in PPS+G2 else "#2A8C6E" if l in UN else GC["G1"] if l in G1+["LWB","CicWB_2022"] else "#555555"
def disp(l): return l.replace("_","-")
CM=ListedColormap(["#F4F4F4","#E8803B","#1D4E89","#B8B8B8","#6B3FA0"])
# ---------- 16S ----------
def fig16():
    d=fa(B+"/04_16S/16S_pruni_set.aln.fasta"); rows=[(n,lab(n)) for n in d if lab(n)]
    ref=[n for n,l in rows if l=="PPS1A"] or [n for n,l in rows if l in PPS]
    ref=ref[0]; R=d[ref].upper(); L=len(R)
    def code(s):
        s=s.upper(); a=np.zeros(L,int)
        for i,(x,y) in enumerate(zip(R,s)):
            if x=="-" and y=="-": a[i]=0
            elif x=="-" or y=="-": a[i]=2
            elif x not in "ACGT" or y not in "ACGT": a[i]=3
            elif x!=y: a[i]=1
        return a
    M=[];T=[]
    for n,l in rows:
        c=code(d[n]); M.append(c); s=d[n].upper()
        cmp_=sum(1 for x,y in zip(R,s) if x in "ACGT" and y in "ACGT")
        T.append((n,l,int((c==1).sum()),int((c==2).sum()),int((c==3).sum()),cmp_,100*(1-(c==1).sum()/max(cmp_,1))))
    order=sorted(range(len(rows)),key=lambda i:(T[i][2]+T[i][3],T[i][1]))
    M=np.array([M[i] for i in order]); T=[T[i] for i in order]
    with open(OUT+"/16S_differences_vs_PPS1A.tsv","w") as f:
        f.write("record\tlabel\tmismatches\tindel_columns\tambiguous\tcompared_sites\tidentity_pct\n")
        for t in T: f.write("\t".join(str(x) if not isinstance(x,float) else "%.2f"%x for x in t)+"\n")
    n=len(T); fig=plt.figure(figsize=(7.2,1.2+0.17*n))
    ax=fig.add_axes([0.22,0.17,0.60,0.78]); ax.imshow(M,cmap=CM,vmin=0,vmax=4,aspect="auto",interpolation="nearest")
    ax.set_yticks(range(n)); cnt={}
    for t in T: cnt[t[1]]=cnt.get(t[1],0)+1
    seen={}
    def nm_(t):
        seen[t[1]]=seen.get(t[1],0)+1
        return disp(t[1])+("" if cnt[t[1]]==1 else " (copy %d)"%seen[t[1]])
    ax.set_yticklabels([nm_(t) for t in T],fontsize=5.5)
    for tl,t in zip(ax.get_yticklabels(),T): tl.set_color(col(t[1]))
    ax.set_xlabel("Alignment position (reference: PPS1A, %d columns)"%L,fontsize=7,labelpad=2); ax.tick_params(axis="x",labelsize=6)
    ax2=ax.twinx(); ax2.set_ylim(ax.get_ylim()); ax2.set_yticks(range(n)); ax2.set_yticklabels(["%d / %d / %.2f%%"%(t[2],t[3],t[6]) for t in T],fontsize=5.5); ax2.tick_params(length=0)
    ax2.set_ylabel("mismatch / indel cols / identity",fontsize=6)
    fig.legend(handles=[Patch(fc=CM(1/4*i),ec="#999",lw=0.4,label=nm) for i,nm in [(0,"identical"),(1,"mismatch"),(2,"indel"),(3,"ambiguous")]],loc="lower center",bbox_to_anchor=(0.5,0.0),ncol=4,fontsize=6.5,frameon=False)
    for e in ("pdf","svg"): fig.savefig(OUT+"/Fig_16S_differences."+e)
    fig.savefig(OUT+"/Fig_16S_differences.png",dpi=220); plt.close(fig); print("16S rows",n,"ref",ref)
    print(open(OUT+"/16S_differences_vs_PPS1A.tsv").read()[:2500])
# ---------- genes ----------
def codon_aln(prot,nt):
    res=sum(1 for c in prot if c!="-")
    if len(nt) not in (3*res,3*res+3): return None
    out=[];k=0
    for c in prot:
        if c=="-": out.append("---")
        else: out.append(nt[3*k:3*k+3].upper()); k+=1
    return out
def figgenes(genes=("rpoB","secY","secA","rps3","rpl22")):
    M=B+"/06_core_genome_tree/markers/"; panels=[];summ=[]
    for g in genes:
        pa=M+g+".aln.faa"
        if os.path.exists(M+g+".pruni.aln.faa"): pa=M+g+".pruni.aln.faa"
        P=fa(pa); N=fa(M+g+".nt.fna"); C={}
        for n in P:
            l=lab(n)
            if l and l!="PBT" and n in N:
                ca=codon_aln(P[n],N[n])
                if ca: C[l]=(P[n],ca)
                else: print("skip (length)",g,n)
        if "PPS1A" not in C: print("no ref",g); continue
        rp,rc=C["PPS1A"]; Lc=len(rp); rows=[]
        for l,(p,c) in C.items():
            a=np.zeros(Lc,int); nd=0; cmp_=0
            for i in range(Lc):
                x,y=rc[i],c[i]
                if x=="---" and y=="---": continue
                if x=="---" or y=="---" or rp[i]=="X" or p[i]=="X": a[i]=3; continue
                if "N" in x or "N" in y: a[i]=3; continue
                d_=sum(1 for u,v in zip(x,y) if u!=v); nd+=d_; cmp_+=3
                if d_: a[i]=1 if rp[i]==p[i] else 2
            rows.append((l,a,nd,cmp_,100*(1-nd/max(cmp_,1)),int((a==2).sum()),int((a==1).sum())))
        pp=[r for r in rows if r[0] in PPS and r[0]!="PPS1A" and r[2]==0]
        drop={r[0] for r in pp}; rows=[r for r in rows if r[0] not in drop]
        nm=len(pp)+1
        rows.sort(key=lambda r:(r[2],r[0])); panels.append((g,Lc,rows,nm))
        for r in rows: summ.append((g,r[0],r[2],r[3],"%.2f"%r[4],r[6],r[5]))
    with open(OUT+"/gene_differences_vs_PPS1A.tsv","w") as f:
        f.write("gene\tgenome\tnt_diffs\tnt_compared\tnt_identity_pct\tsyn_codons\tnonsyn_codons\n")
        for s in summ: f.write("\t".join(map(str,s))+"\n")
    h=[len(p[2]) for p in panels]; fig=plt.figure(figsize=(7.2,1.0+0.15*sum(h)+0.35*len(h)))
    tot=sum(h)+2.5*len(h); y0=0.95
    for g,Lc,rows,nm in panels:
        hh=(len(rows))/tot*0.9; ax=fig.add_axes([0.20,y0-hh,0.58,hh])
        ax.imshow(np.array([r[1] for r in rows]),cmap=CM,vmin=0,vmax=4,aspect="auto",interpolation="nearest")
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels([("PPS strains (n=%d)"%nm if (r[0]=="PPS1A" and nm>1) else disp(r[0])) for r in rows],fontsize=5.5)
        for tl,r in zip(ax.get_yticklabels(),rows): tl.set_color(col(r[0]))
        ax.set_xticks([]); ax.set_title("%s  (%d codons; reference PPS1A)"%(g,Lc),fontsize=7,loc="left")
        ax2=ax.twinx(); ax2.set_ylim(ax.get_ylim()); ax2.set_yticks(range(len(rows)))
        ax2.set_yticklabels(["%d nt, %.2f%%"%(r[2],r[4]) for r in rows],fontsize=5.5); ax2.tick_params(length=0)
        y0-=hh+2.5/tot*0.9
    fig.legend(handles=[Patch(fc="#F4F4F4",ec="#999",lw=0.4,label="identical codon"),Patch(fc="#E8803B",label="synonymous"),Patch(fc="#1D4E89",label="non-synonymous"),Patch(fc="#B8B8B8",label="gap / unresolved")],loc="lower center",ncol=4,fontsize=6.5,frameon=False)
    for e in ("pdf","svg"): fig.savefig(OUT+"/Fig_gene_differences."+e)
    fig.savefig(OUT+"/Fig_gene_differences.png",dpi=220); plt.close(fig)
    print("genes",[p[0] for p in panels]); print(open(OUT+"/gene_differences_vs_PPS1A.tsv").read()[:3500])
if __name__=="__main__":
    try: fig16()
    except Exception as e: print("16S FAILED",repr(e))
    figgenes(); print("done",OUT)
