# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
from Bio import SeqIO
B="/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/"
for r in SeqIO.parse(B+"14_pbt_genbank_refs/curated84_contigs.fna","fasta"):
    if r.id=="CX|NZ_LHCF01000002.1": break
nt=r.seq[16823:20547].reverse_complement(); print("region",len(nt))
p1=str(nt[0:1224].translate(table=11)); p2=str(nt[1270:3724].translate(table=11))
print(len(p1),len(p2),"stops",p1.count("*"),p2.count("*"))
prot=p1+"X"*16+p2; print(len(prot),prot[:20],prot[-10:])
ref=[str(x.seq) for x in SeqIO.parse(B+"06_core_genome_tree/markers/rpoB.faa","fasta") if x.id=="6A1"][0]
print("6A1 len",len(ref),"identical",sum(a==b for a,b in zip(prot,ref) if a!="X"),"of",len(ref)-16)
open(B+"06_core_genome_tree/markers/CX_rpoB_masked.faa","w").write(">CX\n"+prot+"\n")
