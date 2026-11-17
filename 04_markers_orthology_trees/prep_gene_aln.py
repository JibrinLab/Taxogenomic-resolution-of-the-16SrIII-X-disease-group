# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
import subprocess,sys,os
sys.path.insert(0,os.getcwd())
from fig_diff import fa,lab,B
M=B+"/06_core_genome_tree/markers/"
for g in (os.environ.get("GENES","rpoB,secY,secA,rps3,rpl22")).split(","):
    P=fa(M+g+".faa"); keep={n:s for n,s in P.items() if lab(n) and lab(n)!="PBT"}
    open(M+g+".pruni.faa","w").write("".join(">%s\n%s\n"%(n,s) for n,s in keep.items()))
    r=subprocess.run(["mafft","--auto","--quiet",M+g+".pruni.faa"],capture_output=True,text=True)
    open(M+g+".pruni.aln.faa","w").write(r.stdout); print(g,len(keep),r.stdout.count(">"))
