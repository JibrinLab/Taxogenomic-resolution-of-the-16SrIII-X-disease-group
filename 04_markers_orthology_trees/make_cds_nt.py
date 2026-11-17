# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
import glob, os
from Bio.Seq import Seq
out=open("cds_nt.fasta","w"); tot=0
for f in sorted(glob.glob("inputs/*.gff")):
    lab=os.path.basename(f)[:-4]
    gff,fa=open(f).read().split("##FASTA\n",1)
    seq={};k=None
    for l in fa.splitlines():
        if l.startswith(">"): k=l[1:].split()[0]; seq[k]=[]
        else: seq[k].append(l.strip())
    seq={a:"".join(b) for a,b in seq.items()}
    n=0
    for l in gff.splitlines():
        c=l.split("\t")
        if len(c)<9 or c[2]!="CDS" or c[0] not in seq: continue
        s,e=int(c[3]),int(c[4]); ph=int(c[7]) if c[7] in ("0","1","2") else 0
        nt=Seq(seq[c[0]][s-1:e])
        if c[6]=="-": nt=nt.reverse_complement()
        nt=nt[ph:]; nt=nt[:len(nt)//3*3]
        p=str(nt.translate(table=11))
        if p.endswith("*"): nt=nt[:-3]; p=p[:-1]
        if "*" in p or len(p)<30: continue
        n+=1; tot+=1; out.write(f">{lab}_{n:05d}\n{str(nt)}\n")
out.close(); print("CDS written:",tot)
