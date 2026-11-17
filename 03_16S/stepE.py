# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
from lib16s import *
import collections
cop=collections.defaultdict(list)
for h in aln: cop[lab(h)].append(h)
refs=readfa("/projects/dsn011/pps/26_final_markers/16S/phytoplasma_16S_refs.fasta")
hits=sorted(((pid(raw["ChTDIII_c1"],s)[0],h) for h,s in refs.items()),reverse=True)[:6]
print("== ChTDIII_c1 nearest of the %d reference 16S sequences:"%len(refs)); [print("  %.2f  %s"%(i,h[:70])) for i,h in hits]
def posmap(h):
    m={};p=0
    for j,c in enumerate(aln[h]):
        if c!="-": p+=1; m[p]=j
    return m
lw=posmap("LWB_c1"); pp=posmap("PPS1A_c1"); oth=[h for h in aln if lab(h)!="LWB"]
print("\n== LWB positions where both LWB copies share one base carried by at most 3 other copies")
for p,j in lw.items():
    b={aln[h][j] for h in cop["LWB"]}
    if len(b)!=1 or list(b)[0] not in "ACGT": continue
    b=list(b)[0]; sh=[h for h in oth if aln[h][j]==b]; ob=collections.Counter(aln[h][j] for h in oth if aln[h][j] in "ACGT" and aln[h][j]!=b)
    if len(sh)<=3 and ob and ob.most_common(1)[0][1]>=0.9*(len(oth)-len(sh)): print("  LWB pos %d: LWB %s, others %s, shared by %s"%(p,b,ob.most_common(1)[0][0],sh))
for name,mp in (("LWB numbering",lw),("PPS1A numbering",pp)):
    print("\n== positions 182, 412, 413, 617 by",name)
    for p in (182,412,413,617):
        j=mp[p]; c=collections.Counter(aln[h][j] for h in aln)
        print("  %d (column %d): LWB %s | PPS1A %s | all copies %s | LWB base also in %s"%(p,j+1,aln["LWB_c1"][j],aln["PPS1A_c1"][j],dict(c),sorted({lab(h) for h in oth if aln[h][j]==aln["LWB_c1"][j]})[:8]))
