# ANALYSIS: Average nucleotide identity of the 22 genomes plus PG-8A (FastANI, ANIm, ANIb)
# MANUSCRIPT SECTION: Average nucleotide identity
# USED FOR: Fig. 2, Fig. S3.1, Tables S3.1 to S3.5
# FOLDER: 02_ANI/22genomes_PG8A (see README.md in the parent folder)
import csv, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans']
plt.rcParams['pdf.fonttype']=42
plt.rcParams['svg.fonttype']='none'
PG='A_laidlawii'
GROUPS=[('Group 1',0,9),('Group 2',10,15),('Group 3',16,20),('Group 4',21,21)]
def disp(n):
    return 'A. laidlawii PG-8A' if n==PG else n.replace('_','-')
def load(name):
    r=list(csv.reader(open('02_tables/'+name+'_22genomes_symmetrized_mean.tsv'),delimiter='\t'))
    labs=r[0][1:]
    D=np.array([[np.nan if x=='NA' else float(x) for x in row[1:]] for row in r[1:]])
    return labs,D
cmap=__import__('copy').copy(plt.get_cmap('YlGnBu'))
cmap.set_bad('#d9d9d9')
def plot(name):
    labs,D=load(name)
    n=len(labs)
    off=D[:22,:22][~np.eye(22,dtype=bool)]
    vmin=np.floor(off.min())
    vmax=100.0
    fig,ax=plt.subplots(figsize=(15,14))
    im=ax.imshow(np.ma.masked_invalid(D),cmap=cmap,vmin=vmin,vmax=vmax)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels([disp(x) for x in labs],rotation=90,fontsize=12,fontweight='bold',color='black')
    ax.set_yticklabels([disp(x) for x in labs],fontsize=12,fontweight='bold',color='black')
    for t in (ax.get_xticklabels()[-1],ax.get_yticklabels()[-1]):
        t.set_style('italic')
    for i in range(22):
        for j in range(22):
            v=D[i,j]
            c='white' if (v-vmin)/(vmax-vmin)>0.5 else 'black'
            ax.text(j,i,'%.1f'%v,ha='center',va='center',fontsize=7.5,fontweight='bold',color=c)
    ax.set_xticks(np.arange(-0.5,n,1),minor=True)
    ax.set_yticks(np.arange(-0.5,n,1),minor=True)
    ax.grid(which='minor',color='white',lw=0.6)
    ax.tick_params(which='minor',length=0)
    for b in (9.5,15.5,20.5,21.5):
        ax.axhline(b,color='white',lw=3)
        ax.axvline(b,color='white',lw=3)
    for g,a,z in GROUPS:
        ax.plot([a-0.4,z+0.4],[-0.9,-0.9],color='black',lw=2.5,clip_on=False)
    ax.set_xlim(-0.5,n-0.5)
    ax.set_ylim(n-0.5,-0.5)
    cb=fig.colorbar(im,ax=ax,fraction=0.03,pad=0.02)
    cb.set_label(name+' identity (%)',fontsize=14,fontweight='bold')
    cb.ax.tick_params(labelsize=12)
    for l in cb.ax.get_yticklabels():
        l.set_fontweight('bold')
    os.makedirs('01_figures/publication_style',exist_ok=True)
    for ext in ('png','pdf','svg'):
        fig.savefig('01_figures/publication_style/ANI_'+name+'_heatmap_pub.'+ext,dpi=600,bbox_inches='tight')
    plt.close(fig)
for nm in ('ANIm','ANIb','FastANI'):
    plot(nm)
