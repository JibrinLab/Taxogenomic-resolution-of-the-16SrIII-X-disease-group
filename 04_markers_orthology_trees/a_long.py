# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
import glob,os,csv
M="/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/06_core_genome_tree/"
R=["PPS1A","OVSP_L01","CicWB_2022","LWB","CX"]
G={"caryae":"PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1","pruni":"CX 6A1 ZT3_1 PR2021 JR1 PDA15","manihotii":"CicWB_2022 Vc33 ChTDIII OVSP_L01 OVSP_S01","lactucae":"LWB"}
grp={g:k for k,v in G.items() for g in v.split()}
def fa(f):
    d={};k=None
    for l in open(f):
        l=l.strip()
        if l[:1]==">": k=l[1:].split()[0]; d[k]=[]
        elif k: d[k].append(l)
    return {k:"".join(v) for k,v in d.items()}
def idn(a,b):
    p=[(x,y) for x,y in zip(a,b) if x not in "-X*" and y not in "-X*"]
    return (100*sum(x==y for x,y in p)/len(p),len(p)) if p else (None,0)
files=[f for f in sorted(glob.glob(M+"markers/*.aln.faa")) if ".pruni." not in f and "_plus" not in f]+[M+"markers/rpoB_plusCX.aln.faa",M+"concat95.aln.faa"]
w=csv.writer(open("marker_identity_long.tsv","w"),delimiter="\t"); w.writerow(["marker","reference","genome","group","identity","aligned_sites"])
for f in files:
    m=os.path.basename(f).replace(".aln.faa",""); A=fa(f)
    for r in R:
        if r not in A: print("missing",r,"in",m); continue
        for g,s in A.items():
            if g==r: continue
            i,n=idn(A[r],s); w.writerow([m,r,g,grp.get(g,"other"),"" if i is None else round(i,2),n])
print(len(files),"alignments:",[os.path.basename(f).replace(".aln.faa","") for f in files])
