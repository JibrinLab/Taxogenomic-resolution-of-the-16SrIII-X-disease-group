# ANALYSIS: Average nucleotide identity (FastANI, 83 curated and 299 genomes; linkage sensitivity)
# MANUSCRIPT SECTION: Average nucleotide identity
# USED FOR: ANI levels n=83 and n=299; linkage sensitivity
# FOLDER: 02_ANI (see README.md in the parent folder)
import pandas as pd, numpy as np, re, os
from openpyxl.formatting.rule import ColorScaleRule
from openpyxl.styles import Font, PatternFill, Alignment
T="/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/"
ids=[l.strip() for l in open("curated_genome_ids.txt") if l.strip()]
cm=pd.read_csv(T+"01_genome_quality/checkm_with_organism.tsv",sep="\t")
fl=pd.read_csv(T+"01_genome_quality/curated_final_list.tsv",sep="\t")
cand=pd.read_csv(T+"01_genome_quality/curated_candidates.tsv",sep="\t")
x=pd.read_excel("curated83_ANI.xlsx",sheet_name=None,index_col=0)
nn=pd.read_excel("curated83_ANI.xlsx",sheet_name="nearest_and_clusters")
s16=pd.read_csv(T+"04_16S/16S_summary.tsv",sep="\t")
pel=set(open("pel_ids.txt").read().split())
rows=[]
FA="/projects/dsn011/ppss_publication/ani_all299/fastani_all299.tsv"
pat=r"(.+?\.fna)\s*(/.+?\.fna)\s+([\d.]+)\s+(\d+)\s+(\d+)"
for l in open(FA):
    m=re.match(pat,l.strip())
    if m:
        q=os.path.basename(m[1])[:-4]
        r=os.path.basename(m[2])[:-4]
        rows.append((q,r,float(m[3]),int(m[4])/int(m[5])*100))
fa=pd.DataFrame(rows,columns=["q","r","ani","frac"])
fa=fa[fa.q.isin(ids)&fa.r.isin(ids)]
A=fa.pivot_table(index="q",columns="r",values="ani").reindex(index=ids,columns=ids)
F=fa.pivot_table(index="q",columns="r",values="frac").reindex(index=ids,columns=ids)
A=((A+A.T)/2).fillna(A).fillna(A.T).round(2)
F=F.round(1)
keep=["genome","acc","asm","organism","strain","level","date","scaffolds"]
t=fl.merge(cm[keep],on="genome",how="left")
t=t.merge(nn,on="genome",how="left")
s=s16.rename(columns={"n_16S":"n_16S_copies","lengths":"16S_lengths"})
t=t.merge(s[["genome","n_16S_copies","16S_lengths"]],on="genome",how="left")
t["organism"]=t.organism.fillna("Candidatus Phytoplasma sp. (PPS, this study)")
t["strain"]=np.where(t.genome.str.startswith("PPS"),t.genome,t.strain)
t["in_Pellegrinetti_2026"]=t.genome.str.extract(r"(\d{9})")[0].isin(pel)
t["source"]=t.reason.str.rstrip(";")
t=t.drop(columns="reason")
ren={"genome":"Genome","acc":"Assembly accession","asm":"Assembly name"}
ren.update({"organism":"Organism name","strain":"Strain","level":"Assembly level"})
ren.update({"date":"Release date","scaffolds":"Scaffolds"})
ren.update({"Completeness":"CheckM completeness (%)","Contamination":"CheckM contamination (%)"})
ren.update({"ANIm_to_nearest":"ANIm to nearest (%)","nearest_neighbour":"Nearest neighbour"})
t=t.rename(columns=ren)
bad="; ".join(cand[~cand["pass"]].genome)
items=["Quality tool","Inclusion rule","Genomes scored","Pass","Curated set","Excluded","ANI","Clustering","Species boundary","Date"]
vals=["CheckM v1 lineage_wf","completeness >= 90% and contamination <= 5%",len(cm),int(cm["pass"].sum()),len(t),bad]
vals+=["pyANI ANIm and ANIb (both directions averaged); FastANI fragLen 3000 minFraction 0.2","average linkage of ANI matrix","95 to 96% ANI (IRPCM rule c)","2026-10-20"]
crit=pd.DataFrame({"Item":items,"Value":vals})
ex=cand[~cand["pass"]][["genome","Completeness","Contamination","reason"]]
w=pd.ExcelWriter("Curated83_Supplementary.xlsx",engine="openpyxl")
crit.to_excel(w,sheet_name="Criteria",index=False)
t.to_excel(w,sheet_name="Curated83",index=False)
ex.to_excel(w,sheet_name="Excluded",index=False)
cm.drop(columns=["in_cand"],errors="ignore").to_excel(w,sheet_name="CheckM_all299",index=False)
x["ANIm_matrix"].to_excel(w,sheet_name="ANIm_ANI")
x["ANIb_matrix"].to_excel(w,sheet_name="ANIb_ANI")
x["ANIm_aligned_pct"].to_excel(w,sheet_name="ANIm_aligned_pct")
A.to_excel(w,sheet_name="FastANI_ANI")
F.to_excel(w,sheet_name="FastANI_fragment_pct")
hf=PatternFill("solid",fgColor="1F3864")
for ws in w.book.worksheets:
    for c in ws[1]:
        c.font=Font(bold=True,color="FFFFFF")
        c.fill=hf
        c.alignment=Alignment(wrap_text=True,vertical="top")
    ws.freeze_panes="B2"
    for col in ws.columns:
        n=max(len(str(c.value or "")) for c in col[:60])
        ws.column_dimensions[col[0].column_letter].width=min(40,max(10,n+2))
    if ws.title in ("ANIm_ANI","ANIb_ANI","FastANI_ANI"):
        rule=ColorScaleRule(start_type="num",start_value=90,start_color="FFFFFF",end_type="num",end_value=100,end_color="2F6FB0")
        ws.conditional_formatting.add(ws.dimensions,rule)
w.close()
print("rows:",len(t),"| FastANI pairs kept:",len(fa),"| in Pellegrinetti:",int(t.in_Pellegrinetti_2026.sum()))
