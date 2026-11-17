# ANALYSIS: Average nucleotide identity of the 22 genomes plus PG-8A (FastANI, ANIm, ANIb)
# MANUSCRIPT SECTION: Average nucleotide identity
# USED FOR: Fig. 2, Fig. S3.1, Tables S3.1 to S3.5
# FOLDER: 02_ANI/22genomes_PG8A (see README.md in the parent folder)
import csv, itertools, os, re
G={'PPS':['PPS1A','PPS2B','PPS6A','PPS6B','PPS9B','PPS10A','CanS4','MW1','VAC','MA1'],
   'CX':['CX','6A1','ZT3-1','PR2021','JR1','PDA15'],
   'CicWB':['CicWB-2022','Vc33','ChTDIII','OVSP_L01','OVSP_S01'],
   'LWB':['LWB']}
grp={g:k for k,v in G.items() for g in v}
def tab(p):
    r=list(csv.reader(open(p),delimiter='\t')); h=r[0][1:]
    return {(row[0],h[j]):float(x)*100 for row in r[1:] for j,x in enumerate(row[1:])}
def fani(p):
    d={}
    for l in open(p):
        a,b,v=l.split()[:3]
        d[(os.path.basename(a).rsplit('.',1)[0],os.path.basename(b).rsplit('.',1)[0])]=float(v)
    return d
def stats(name,m):
    ks=list(grp); w=[];b=[]
    for x,y in itertools.combinations(ks,2):
        if (x,y) in m and (y,x) in m:
            v=(m[(x,y)]+m[(y,x)])/2
            (w if grp[x]==grp[y] else b).append((v,x,y))
    print(f'{name}: pairs={len(w)+len(b)} within_min={min(w)} between_max={max(b)}')
A=tab('pyani_anim/ANIm_percentage_identity.tab'); stats('ANIm',A)
B=tab('pyani_anib/ANIb_percentage_identity.tab'); stats('ANIb',B)
F=fani('fastani_23.tsv'); stats('FastANI',F)
print('FastANI rows with A_laidlawii:',[k for k in F if 'A_laidlawii' in k and k[0]!=k[1]])
C=tab('pyani_anib/ANIb_alignment_coverage.tab')
print('ANIb PG-8A coverage range (%):',min(v for (x,y),v in C.items() if 'A_laidlawii' in (x,y) and x!=y),max(v for (x,y),v in C.items() if 'A_laidlawii' in (x,y) and x!=y))
print('ANIb PG-8A identity range (%):',min(v for (x,y),v in B.items() if 'A_laidlawii' in (x,y) and x!=y),max(v for (x,y),v in B.items() if 'A_laidlawii' in (x,y) and x!=y))
print('MA1-PDA15 new ANIm:',A[('MA1','PDA15')],A[('PDA15','MA1')])
r=list(csv.reader(open('../ani_all299/pyani_all299_anim/ANIm_percentage_identity.tab'),delimiter='\t')); h=r[0][1:]
for row in r[1:]:
    if row[0].startswith('MA1') or row[0]=='MA1':
        for j,c in enumerate(h):
            if c.startswith('PDA15'): print('299 run',row[0],c,float(row[j+1])*100)
