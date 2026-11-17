# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
import collections
import pandas as pd
P=pd.read_csv("protein_products.tsv",sep="\t").fillna("")
prod=dict(zip(P.gene,P["product"]))
t=pd.read_csv("proteins/OrthoFinder/Results_Oct20/Orthogroups/Orthogroups.tsv",sep="\t",index_col=0).fillna("")
out=[]
for og,row in t.iterrows():
    genes=[g for c in t.columns for g in row[c].split(", ") if g]
    top,n=collections.Counter(prod.get(g,"") for g in genes).most_common(1)[0]
    out.append((og,top,n,len(genes),sum(1 for c in t.columns if row[c])))
O=pd.DataFrame(out,columns=["OG","top_product","n_top","n_genes","n_genomes"])
pat=r"RNA polymerase subunit beta$|SecY|SecA|ribosomal protein"
M=O[O.top_product.str.contains(pat,case=False)&~O.top_product.str.contains("beta'|prime",case=False)]
M.to_csv("marker_orthogroups.tsv",sep="\t",index=False)
print("marker orthogroups:",len(M))
print(M.sort_values("n_genomes",ascending=False).head(50).to_string(index=False))
