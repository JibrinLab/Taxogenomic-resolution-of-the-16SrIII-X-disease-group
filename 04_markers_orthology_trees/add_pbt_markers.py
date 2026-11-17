# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
import subprocess,os
T="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
B="TCAG"
cod={a+b+c:T[i] for i,(a,b,c) in enumerate([(a,b,c) for a in B for b in B for c in B])}
def rd(f): return "".join(l.strip() for l in open(f) if l[0]!=">").upper()
def tr(s): return "".join(cod.get(s[i:i+3],"X") for i in range(0,len(s)-2,3))
def fa(f):
    d={};k=None
    for l in open(f):
        l=l.rstrip()
        if l.startswith(">"): k=l[1:].split()[0]; d[k]=""
        elif k: d[k]+=l
    return d
src={"GU004361":["secY"],"EF186807":["rps3","rpl22"]}
with open("pbt_orfs.faa","w") as o:
    for acc in src:
        s=rd(acc+".fa")
        for fr in range(3):
            pos=0
            for seg in tr(s[fr:]).split("*"):
                if len(seg)>=40: o.write(f">{acc}|{fr}|{pos}\n{seg}\n")
                pos+=len(seg)+1
orf=fa("pbt_orfs.faa")
env=dict(os.environ,PATH="/projects/dsn011/urbello_conda/envs/mafft_env/bin:"+os.environ["PATH"])
for acc,ms in src.items():
    for m in ms:
        subprocess.run(["hmmsearch","--noali","--tblout",f"{m}.tbl",f"../11_tree_all_genomes/hmm/{m}.hmm","pbt_orfs.faa"],capture_output=True)
        hits=[l.split() for l in open(f"{m}.tbl") if l[0]!="#"]
        hits=[h for h in hits if h[0].startswith(acc)]
        if not hits: print(m,"no ORF hit"); continue
        h=max(hits,key=lambda x:float(x[5]))
        nm=f"PBT_{acc}_{m}"
        open(f"{m}_pbt.faa","w").write(f">{nm}\n{orf[h[0]]}\n")
        base=f"../06_core_genome_tree/markers/{m}.aln.faa"
        out=subprocess.run(["mafft","--add",f"{m}_pbt.faa","--keeplength","--thread","8",base],capture_output=True,text=True,env=env).stdout
        open(f"{m}_plusPBT.aln.faa","w").write(out)
        A=fa(f"{m}_plusPBT.aln.faa"); q=A[nm]; r=[]
        for g,s in A.items():
            if g==nm: continue
            p=[(x,y) for x,y in zip(q,s) if x not in "-X" and y not in "-X"]
            if p: r.append((sum(x==y for x,y in p)/len(p),len(p),g))
        r.sort(reverse=True)
        print(m,"ORF",h[0],"hmm score",h[5],"aligned aa",len(q.replace("-","")))
        for i,n,g in r[:8]: print(f"   {g:10s} {100*i:6.2f}% over {n} aa")
