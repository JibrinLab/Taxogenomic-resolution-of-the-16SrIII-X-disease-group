# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
from lib16s import *
refs=readfa("/projects/dsn011/pps/26_final_markers/16S/phytoplasma_16S_refs.fasta")
s=refs[[k for k in refs if "JQ044393" in k][0]]
out=open(D+"/04_16S/all_described_species_identity_combined.tsv","a")
for q in ("PPS1A_c1","CX_c1","OVSP_L01_c1","CicWB_2022_c1","LWB_c1"):
    i,n,st=pid(raw[q],s); out.write("%s\tpruni\tJQ044393.1\t%.2f\t%d\tPX11Ct1\n"%(lab(q),i,n)); print(lab(q),"%.2f%% (%d sites)"%(i,n))
out.close()
