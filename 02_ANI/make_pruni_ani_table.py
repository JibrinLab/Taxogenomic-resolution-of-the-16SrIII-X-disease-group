# ANALYSIS: Average nucleotide identity (FastANI, 83 curated and 299 genomes; linkage sensitivity)
# MANUSCRIPT SECTION: Average nucleotide identity
# USED FOR: ANI levels n=83 and n=299; linkage sensitivity
# FOLDER: 02_ANI (see README.md in the parent folder)
import pandas as pd, glob, itertools, os, sys
B="/projects/dsn011/ppss_publication/ani_all299/"
def load(m, kind):
    f=glob.glob(f"{B}pyani_all299_{m}/{m.upper().replace('ANIM','ANIm').replace('ANIB','ANIb')}_{kind}.tab")
    if not f: sys.exit(f"missing {m} {kind}")
    d=pd.read_csv(f[0],sep="\t",index_col=0); d.columns=d.columns.astype(str); d.index=d.index.astype(str)
    return d
m=load("anim","percentage_identity"); b=load("anib","percentage_identity")
mc=load("anim","alignment_coverage")
G1={"CX":"001277135","Vc33":"001623385","LWB":"051397825","JR1":"000309465","CicWB2022":"035853675",
    "6A1":"033391615","ZT3-1":"024586395","PR2021":"029746895","PDA15":"041984945"}
G2={"CanS4":"050286905","MW1":"000309485","VAC":"000309405","MA1":"000300695"}
OTH={"ChTDIII":"013391955","OVSP_L01":"059597165","OVSP_S01":"059597155"}
CODE={"CX":"CX group","6A1":"CX group","ZT3-1":"CX group","PR2021":"CX group","JR1":"CX group","PDA15":"CX group","CicWB2022":"CicWB-2022 group","Vc33":"CicWB-2022 group","ChTDIII":"CicWB-2022 group","OVSP_L01":"CicWB-2022 group","OVSP_S01":"CicWB-2022 group","LWB":"LWB","CanS4":"CanS4/MW1/VAC/MA1","MW1":"CanS4/MW1/VAC/MA1","VAC":"CanS4/MW1/VAC/MA1","MA1":"CanS4/MW1/VAC/MA1"}
def find(k):
    h=[c for c in m.index if k in c]; return h[0] if h else None
grp={}
for n,k in G1.items(): grp[n]=(CODE[n],find(k))
for n,k in G2.items(): grp[n]=(CODE[n],find(k))
for n,k in OTH.items(): grp[n]=(CODE[n],find(k))
for c in m.index:
    if c.startswith("PPS"): grp[c]=("PPS strains",c)
miss=[n for n,(g,c) in grp.items() if c is None]
print("not found:",miss)
grp={n:v for n,v in grp.items() if v[1]}
rows=[]
for a,b_ in itertools.combinations(grp,2):
    ga,ca=grp[a]; gb,cb=grp[b_]
    anim=(m.loc[ca,cb]+m.loc[cb,ca])/2*100 if m.loc[ca,cb]<=1.0 else (m.loc[ca,cb]+m.loc[cb,ca])/2
    anib=(b.loc[ca,cb]+b.loc[cb,ca])/2*100 if b.loc[ca,cb]<=1.0 else (b.loc[ca,cb]+b.loc[cb,ca])/2
    cov=min(mc.loc[ca,cb],mc.loc[cb,ca]); cov=cov*100
    cat=" vs ".join(sorted([ga,gb]))
    rows.append((a,b_,cat,round(anim,2),round(anib,2),round(cov,1)))
df=pd.DataFrame(rows,columns=["A","B","Category","ANIm","ANIb","min_ANIm_aligned_%"])
print(df.groupby("Category")[["ANIm","ANIb"]].agg(["count","min","median","max"]).round(2).to_string())
print("\nPPS vs each reference, ANIm (rows sorted high to low):")
p=df[df.Category.str.contains("PPS")].sort_values("ANIm",ascending=False)
print(p.head(40).to_string(index=False))
with pd.ExcelWriter("pruni_ANI_tables.xlsx") as w:
    df.to_excel(w,"pairs",index=False)
    for v,nm in ((m,"ANIm"),(b,"ANIb")):
        ids=[c for _,c in grp.values()]; v.loc[ids,ids].to_excel(w,nm+"_matrix")
