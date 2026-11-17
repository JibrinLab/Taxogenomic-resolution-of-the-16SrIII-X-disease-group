# ANALYSIS: Average nucleotide identity (FastANI, 83 curated and 299 genomes; linkage sensitivity)
# MANUSCRIPT SECTION: Average nucleotide identity
# USED FOR: ANI levels n=83 and n=299; linkage sensitivity
# FOLDER: 02_ANI (see README.md in the parent folder)
import pandas as pd, numpy as np
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform
B="/projects/dsn011/ppss_publication/ani_all299/"
ids=[l.strip() for l in open("curated_genome_ids.txt") if l.strip()]
def load(m,k):
    d=pd.read_csv(f"{B}pyani_all299_{m}/{'ANIm' if m=='anim' else 'ANIb'}_{k}.tab",sep="\t",index_col=0)
    d.index=d.index.astype(str); d.columns=d.columns.astype(str)
    return d*100 if d.values.max()<=1.5 else d
mi=load("anim","percentage_identity"); bi=load("anib","percentage_identity"); mc=load("anim","alignment_coverage")
miss=[i for i in ids if i not in mi.index]; print("missing from matrices:",miss)
ids=[i for i in ids if i in mi.index]
sym=lambda d:((d.loc[ids,ids]+d.loc[ids,ids].T)/2).round(2)
ANIm=sym(mi); ANIb=sym(bi); COV=mc.loc[ids,ids].round(1)
rows=[]
for nm,M in (("ANIm",ANIm),("ANIb",ANIb)):
    A=M.values.copy(); np.fill_diagonal(A,100)
    Z=linkage(squareform(100-A,checks=False),"average")
    for t in (95,96):
        cl=fcluster(Z,100-t,"distance")
        print(f"{nm} @{t}%: {len(set(cl))} clusters from {len(ids)} genomes")
        M.attrs[t]=cl
nn=[]
for i,g in enumerate(ids):
    r=ANIm.loc[g].drop(g); nn.append((g,r.idxmax(),r.max(),ANIm.attrs[95][i],ANIm.attrs[96][i],ANIb.attrs[95][i],ANIb.attrs[96][i]))
NN=pd.DataFrame(nn,columns=["genome","nearest_neighbour","ANIm_to_nearest","ANIm_cluster95","ANIm_cluster96","ANIb_cluster95","ANIb_cluster96"])
with pd.ExcelWriter("curated83_ANI.xlsx") as w:
    NN.to_excel(w,"nearest_and_clusters",index=False)
    ANIm.to_excel(w,"ANIm_matrix"); ANIb.to_excel(w,"ANIb_matrix"); COV.to_excel(w,"ANIm_aligned_pct")
cl=NN.groupby("ANIm_cluster95").genome.apply(lambda s:[x.split('_',2)[-1] if x.startswith('GC') else x for x in s])
print("\nANIm 95% clusters with 1 member (singletons):")
print(cl[cl.apply(len)==1].to_string())
print("\nlargest clusters:"); print(cl.apply(len).sort_values(ascending=False).head(8).to_string())
