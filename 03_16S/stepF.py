# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
from lib16s import *
import collections
G={}
for x in "PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1".split(): G[x]="group1"
for x in "CX 6A1 ZT3_1 PR2021 JR1 PDA15".split(): G[x]="group2"
for x in "CicWB_2022 Vc33 ChTDIII OVSP_L01 OVSP_S01".split(): G[x]="group3"
G["LWB"]="LWB"
cop=collections.defaultdict(list)
for h in aln: cop[lab(h)].append(h)
sel=[h for h in aln if lab(h) in G]
pp={};p=0
for j,c in enumerate(aln["PPS1A_c1"]):
    if c!="-": p+=1; pp[p]=j
refs=readfa("/projects/dsn011/pps/26_final_markers/16S/phytoplasma_16S_refs.fasta")
print("== identity of the five queries to the Ca. P. pruni reference sequences")
for h,s in refs.items():
    if "JQ044393" in h or "JQ044385" in h: print(" ",h[:60]," | ".join("%s %.2f"%(q.split("_c")[0],pid(raw[q],s)[0]) for q in ("PPS1A_c1","OVSP_L01_c1","CicWB_2022_c1","LWB_c1","CX_c1")))
print("\n== bases at the four positions, PPS1A numbering, copies per group (22 genomes)")
for p in (182,412,413,617):
    j=pp[p]; print("  %d (PPS1A %s):"%(p,aln["PPS1A_c1"][j])," | ".join("%s %s"%(g,dict(collections.Counter(aln[h][j] for h in sel if G[lab(h)]==g))) for g in ("group1","group2","group3","LWB")))
print("\n== PPS1A-numbered positions where both LWB copies carry a base absent from all other genomes of the 22")
odd=("ChTDIII_c1","OVSP_L01_c2"); res={}
for name,ex in (("all copies",()),("without ChTDIII_c1 and OVSP_L01_c2",odd)):
    out=[]
    for p,j in pp.items():
        b={aln[h][j] for h in cop["LWB"]}
        if len(b)!=1 or list(b)[0] not in "ACGT": continue
        b=list(b)[0]; o=[h for h in sel if lab(h)!="LWB" and h not in ex]
        if all(aln[h][j] in "ACGT" and aln[h][j]!=b for h in o): out.append((p,b))
    res[name]=out; print(" ",name,":"," ".join("%d:%s>%s"%(p,aln["PPS1A_c1"][pp[p]],b) for p,b in out))
cand={};h=None
for l in open("/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/26_jibrin_update_2026-11-12/data/16SrIII-X_candidates.fa"):
    l=l.strip()
    if l.startswith(">"): h=l[1:].split()[0]; cand[h]=[]
    elif h: cand[h].append(l.upper())
print("\n== published 16SrIII-X-like sequences at the four positions and at the LWB-specific positions")
want=sorted({p for p,_ in res["without ChTDIII_c1 and OVSP_L01_c2"]}|{182,412,413,617})
for h in ("KC412026.1","MH497013.1"):
    s="".join(cand[h]); best=None
    for st,bb in (("+",s),("-",rc(s))):
        x=al.align(raw["PPS1A_c1"],bb)[0]
        if best is None or x.score>best[0]: best=(x.score,x,bb)
    _,x,bb=best; m={}
    for (ts,te),(qs,qe) in zip(*x.aligned):
        for k in range(te-ts): m[ts+k+1]=bb[qs+k]
    print(" ",h,":"," ".join("%d:%s"%(p,m.get(p,"-")) for p in want))
