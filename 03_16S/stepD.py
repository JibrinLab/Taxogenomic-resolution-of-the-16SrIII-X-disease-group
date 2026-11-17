# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
from lib16s import *
import collections
def colid(a,b):
    n=m=0
    for x,y in zip(a,b):
        if x in "ACGT" and y in "ACGT": n+=1; m+=x!=y
    return (100*(1-m/n),n) if n else (0.0,0)
cop=collections.defaultdict(list)
for h in aln: cop[lab(h)].append(h)
G3=("ChTDIII","CicWB_2022","Vc33","OVSP_L01","OVSP_S01")
print("== 16S copies of group 3 vs each genome (min-max over its copies)")
for g in G3:
    for h in sorted(cop[g]):
        print(" ",h,"|"," | ".join("%s %.1f-%.1f"%(o,min(v),max(v)) for o in ("PPS1A","CX","PR2021","LWB")+G3 if o!=g and cop.get(o) for v in [[colid(aln[h],aln[p])[0] for p in cop[o]]]))
print("\n== mismatch windows (alignment columns, per 100)")
for a,b in (("OVSP_L01_c2","OVSP_L01_c1"),("ChTDIII_c1","CicWB_2022_c1")):
    w=collections.Counter(j//100 for j,(x,y) in enumerate(zip(aln[a],aln[b])) if x in "ACGT" and y in "ACGT" and x!=y)
    print(" ",a,"vs",b,":"," ".join("%d:%d"%(k*100,v) for k,v in sorted(w.items())))
ref=aln["LWB_c1"]; lwpos={}; p=0
for j,c in enumerate(ref):
    if c!="-": p+=1; lwpos[j]=p
oth=[h for h in aln if lab(h)!="LWB"]; diag=[]
for j,p in lwpos.items():
    b={aln[h][j] for h in cop["LWB"]}
    if len(b)!=1 or list(b)[0] not in "ACGT": continue
    b=list(b)[0]; ob=collections.Counter(aln[h][j] for h in oth if aln[h][j] in "ACGT")
    if ob and b not in ob and ob.most_common(1)[0][1]/sum(ob.values())>=0.9: diag.append((p,b,ob.most_common(1)[0][0],j))
print("\n== LWB-specific 16S positions (LWB position, LWB base, base in all other genomes):",[(p,b,o) for p,b,o,_ in diag])
fa={};h=None
for l in open("/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/26_jibrin_update_2026-11-12/data/16SrIII-X_candidates.fa"):
    l=l.strip()
    if l.startswith(">"): h=l[1:].split()[0]; fa[h]=[]
    elif h: fa[h].append(l.upper())
print("\n== published 16SrIII-X candidates: identity to the five queries, and base at each LWB-specific position")
for h,s in fa.items():
    s="".join(s); best=None
    for st,bb in (("+",s),("-",rc(s))):
        x=al.align(raw["LWB_c1"],bb)[0]
        if best is None or x.score>best[0]: best=(x.score,x,bb)
    _,x,bb=best; m={}
    for (ts,te),(qs,qe) in zip(*x.aligned):
        for k in range(te-ts): m[ts+k+1]=bb[qs+k]
    ids=" ".join("%s %.2f"%(q.split("_c")[0],pid(raw[q],s)[0]) for q in ("PPS1A_c1","CX_c1","OVSP_L01_c1","CicWB_2022_c1","LWB_c1"))
    print(" ",h,len(s),"bp |",ids,"| bases:"," ".join("%d:%s"%(p,m.get(p,"-")) for p,_,_,_ in diag))
