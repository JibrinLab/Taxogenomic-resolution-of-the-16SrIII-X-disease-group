# ANALYSIS: Average nucleotide identity of the 22 genomes plus PG-8A (FastANI, ANIm, ANIb)
# MANUSCRIPT SECTION: Average nucleotide identity
# USED FOR: Fig. 2, Fig. S3.1, Tables S3.1 to S3.5
# FOLDER: 02_ANI/22genomes_PG8A (see README.md in the parent folder)
import pandas as pd, glob, re, itertools
G={'1':'PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1','2':'CX 6A1 ZT3-1 PR2021 JR1 PDA15','3':'CicWB-2022 ChTDIII OVSP-L01 OVSP-S01 Vc33','4':'LWB'}
N=lambda x: re.sub(r'[^A-Z0-9]','',x.upper())
K={N(g):k for k,v in G.items() for g in v.split()}
def grp(l):
    n=N(l)
    if n in K: return K[n]
    c=[K[k] for k in K if n.startswith(k)]
    return c[0] if len(set(c))==1 else None
for f in sorted(glob.glob('/projects/dsn011/ppss_publication/ani_23_with_PG8A/tables/*_22genomes_symmetrized_mean.tsv')):
    m=pd.read_csv(f,sep='\t',index_col=0)
    if m.values.max()<=1.0001: m=m*100
    g={l:grp(l) for l in m.index}
    bad=[l for l in g if g[l] is None]
    print('\n=====',f.split('/')[-1],'genomes',len(m),'unmatched',bad)
    for a,b in itertools.combinations_with_replacement('1234',2):
        A=[l for l in m.index if g[l]==a]; B=[l for l in m.columns if g[l]==b]
        pr=[(m.loc[x,y],x,y) for x in A for y in B if x!=y and (a!=b or x<y)]
        if not pr: print(f'cluster {a} vs {b}: not applicable (one genome)'); continue
        lo=min(pr); hi=max(pr)
        print(f'cluster {a} vs {b}: {lo[0]:.2f}-{hi[0]:.2f}  n={len(pr)}  min {lo[1]}/{lo[2]}  max {hi[1]}/{hi[2]}')
