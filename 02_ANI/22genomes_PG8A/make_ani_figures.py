# ANALYSIS: Average nucleotide identity of the 22 genomes plus PG-8A (FastANI, ANIm, ANIb)
# MANUSCRIPT SECTION: Average nucleotide identity
# USED FOR: Fig. 2, Fig. S3.1, Tables S3.1 to S3.5
# FOLDER: 02_ANI/22genomes_PG8A (see README.md in the parent folder)
import csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
G=[['PPS1A','PPS2B','PPS6A','PPS6B','PPS9B','PPS10A','CanS4','MW1','VAC','MA1'],
   ['CX','6A1','ZT3-1','PR2021','JR1','PDA15'],
   ['CicWB-2022','Vc33','ChTDIII','OVSP_L01','OVSP_S01'],
   ['LWB']]
PG='A_laidlawii'
order=[g for grp in G for g in grp]
allg=order+[PG]
def disp(n):
    return 'A. laidlawii PG-8A' if n==PG else n.replace('_','-')
def tab(p):
    r=list(csv.reader(open(p),delimiter='\t'))
    h=r[0][1:]
    return {(row[0],h[j]):float(x)*100 for row in r[1:] for j,x in enumerate(row[1:])}
def fani(p):
    d={}
    for l in open(p):
        a,b,v=l.split()[:3]
        ka=os.path.basename(a).rsplit('.',1)[0]
        kb=os.path.basename(b).rsplit('.',1)[0]
        d[(ka,kb)]=float(v)
    return d
src={'ANIm':tab('pyani_anim/ANIm_percentage_identity.tab'),
     'ANIb':tab('pyani_anib/ANIb_percentage_identity.tab'),
     'FastANI':fani('fastani_23.tsv')}
n=len(allg)
cmap=LinearSegmentedColormap.from_list('w2b',['#f7fbff','#9ecae1','#08306b'])
cmap.set_bad('#d9d9d9')
bounds=np.cumsum([len(g) for g in G])
for name,m in src.items():
    D=np.full((n,n),np.nan)
    for i,a in enumerate(order):
        for j,b in enumerate(order):
            D[i,j]=100.0 if a==b else (m[(a,b)]+m[(b,a)])/2
    f=open(f'tables/{name}_22genomes_symmetrized_mean.tsv','w')
    f.write('\t'+'\t'.join(allg)+'\n')
    for i,a in enumerate(allg):
        f.write(a+'\t'+'\t'.join('NA' if np.isnan(x) else f'{x:.4f}' for x in D[i])+'\n')
    f.close()
    off=D[:22,:22][~np.eye(22,dtype=bool)]
    vmin=np.floor(off.min())
    fig,ax=plt.subplots(figsize=(10,9))
    im=ax.imshow(np.ma.masked_invalid(D),cmap=cmap,vmin=vmin,vmax=100)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels([disp(x) for x in allg],rotation=90,fontsize=8,color='black')
    ax.set_yticklabels([disp(x) for x in allg],fontsize=8,color='black')
    for i in range(22):
        for j in range(22):
            v=D[i,j]
            dark=(v-vmin)/(100-vmin)>0.6
            ax.text(j,i,f'{v:.1f}',ha='center',va='center',fontsize=4.5,color='white' if dark else 'black')
    for b in list(bounds[:-1])+[22]:
        ax.axhline(b-0.5,color='white',lw=1.8)
        ax.axvline(b-0.5,color='white',lw=1.8)
    cb=fig.colorbar(im,ax=ax,fraction=0.035,pad=0.02)
    cb.set_label(name+' identity (%), mean of both directions',color='black')
    ax.set_title(name+', 22 pruni and pruni-related genomes + A. laidlawii PG-8A (outgroup)',fontsize=10,color='black')
    fig.text(0.02,0.01,'Grey: ANI not reportable for PG-8A (see tables/PG8A_pairs_summary.tsv).',fontsize=7,color='black')
    fig.tight_layout(rect=(0,0.02,1,1))
    for ext in ('png','pdf'):
        fig.savefig('figures/ANI_'+name+'_heatmap.'+ext,dpi=300)
    plt.close(fig)
Am=src['ANIm']
Ab=src['ANIb']
F=src['FastANI']
Ac=tab('pyani_anim/ANIm_alignment_coverage.tab')
Al=tab('pyani_anim/ANIm_alignment_lengths.tab')
Bc=tab('pyani_anib/ANIb_alignment_coverage.tab')
hdr='genome ANIm_id_PG8A_query ANIm_id_PG8A_subject ANIm_cov_PG8A_query ANIm_cov_PG8A_subject ANIm_aligned_bp ANIb_id_PG8A_query ANIb_id_PG8A_subject ANIb_cov_PG8A_query ANIb_cov_PG8A_subject FastANI'
f=open('tables/PG8A_pairs_summary.tsv','w')
f.write(hdr.replace(' ','\t')+'\n')
for g in order:
    row=[Am[(PG,g)],Am[(g,PG)],Ac[(PG,g)],Ac[(g,PG)],Al[(PG,g)]/100,Ab[(PG,g)],Ab[(g,PG)],Bc[(PG,g)],Bc[(g,PG)]]
    fv=(F.get((PG,g)),F.get((g,PG)))
    f.write(g+'\t'+'\t'.join('%.4f'%x for x in row)+'\t'+('no_value' if fv==(None,None) else str(fv))+'\n')
f.close()
