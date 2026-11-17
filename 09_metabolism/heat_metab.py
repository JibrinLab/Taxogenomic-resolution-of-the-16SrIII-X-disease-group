# ANALYSIS: Metabolic gene panel (66 genes, 7 pathways) in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 8, Tables S7.1 to S7.10
# FOLDER: 09_metabolism (see README.md in the parent folder)
import sys,pandas as pd,numpy as np,matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch,Rectangle
if len(sys.argv)<3: sys.exit("usage: heat_metab.py presence_absence_matrix.tsv OUTPREFIX")
df=pd.read_csv(sys.argv[1],sep="\t",dtype=str).fillna("")
GR=[("Group 1","#C2410C",["PPS1A","PPS2B","PPS6A","PPS6B","PPS9B","PPS10A","CanS4","MW1","VAC","MA1"]),
    ("Group 2","#2F6DB5",["CX","6A1","ZT3-1","PR2021","JR1","PDA15"]),
    ("Group 3","#2A8C6E",["CicWB-2022","ChTDIII","OVSP-L01","OVSP-S01","Vc33"]),("Group 4","#7E3F98",["LWB"])]
order=[s for _,_,m in GR for s in m]
miss=[s for s in order if s not in df.columns]
if miss: sys.exit("columns missing: %s ; present: %s"%(miss,list(df.columns)))
raw=df[order].apply(lambda c:c.str.strip())
print("distinct cell values:",sorted(set(raw.values.ravel())))
det=raw.apply(lambda c:c.map(lambda v: 0 if v.lower() in ("0","","absent","not detected","no") else 1)).values
path=df["Pathway"].tolist(); gene=df["Gene"].tolist()
n=len(gene); m=len(order); gap=1.1
y=[];cur=1.0;prev=None;starts={}
for i,p in enumerate(path):
    if p!=prev:
        if prev is not None: cur+=gap
        starts[p]=cur-1.0
    y.append(cur); cur+=1; prev=p
H=0.17*cur+1.9; W=7.6
fig=plt.figure(figsize=(W,H))
L,R=0.13,0.97; top=1-1.05/H; bot=0.9/H
ax=fig.add_axes([L,bot,R-L,top-bot-0.012])
for i in range(n):
    for j in range(m):
        ax.add_patch(Rectangle((j,y[i]),1,1,fc="#1D4E89" if det[i,j] else "#F2F2F2",ec="white",lw=0.6))
ax.set_xlim(0,m); ax.set_ylim(cur,0)
ax.set_yticks([v+0.5 for v in y]); ax.set_yticklabels(gene,fontsize=6.5,color="black",style="italic")
ax.set_xticks([j+0.5 for j in range(m)]); ax.set_xticklabels(order,rotation=90,fontsize=7,color="black")
ax.tick_params(length=0)
for sp in ax.spines.values(): sp.set_visible(False)
sa=fig.add_axes([L,top-0.004,R-L,0.014]); sa.set_xlim(0,m); sa.set_ylim(0,1); sa.axis("off")
x=0
for g,c,mem in GR:
    sa.add_patch(Rectangle((x,0),len(mem),1,fc=c,ec="none"))
    sa.text(x+len(mem)/2,1.6,g,ha="center",va="bottom",fontsize=7.5,color="black")
    x+=len(mem)
    if x<m: ax.axvline(x,color="black",lw=0.6)
for p,s0 in starts.items(): ax.text(-0.1,s0+0.9,p,ha="left",va="bottom",fontsize=7.5,fontweight="bold",color="black",transform=ax.transData,clip_on=False)
fig.text(0.02,1-0.18/H,"Metabolic Pathway gene content of 22 pruni and pruni-related genomes",fontsize=9.5,fontweight="bold",va="top")
fig.legend(handles=[Patch(fc="#1D4E89",label="Detected"),Patch(fc="#F2F2F2",ec="#888",label="Not detected")],loc="lower center",bbox_to_anchor=(0.5,0.0),ncol=2,frameon=False,fontsize=8)
for e in ("pdf","svg"): fig.savefig(sys.argv[2]+"."+e)
fig.savefig(sys.argv[2]+".png",dpi=300); print("wrote",sys.argv[2],"genes",n)
