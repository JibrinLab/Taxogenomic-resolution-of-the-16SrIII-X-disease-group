# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
import collections, glob, os, subprocess
R="proteins/OrthoFinder/Results_Oct20/"
labs=[os.path.basename(f)[:-4] for f in glob.glob("inputs/*.gff")]
ing=[l for l in labs if l!="A_laidlawii"]; N=len(ing)
env=dict(os.environ,PATH="/projects/dsn011/urbello_conda/envs/mafft_env/bin:"+os.environ["PATH"])
sel=[]
for l in open(R+"Orthogroups/Orthogroups.txt"):
    k,v=l.rstrip("\n").split(": ",1)
    cn=collections.Counter(x.rsplit("_",1)[0] for x in v.split())
    pres=sum(1 for s in ing if cn.get(s,0)>0)
    if pres>=0.95*N and all(cn.get(s,0)<=1 for s in ing) and cn.get("A_laidlawii",0)==1: sel.append(k)
print("selected OGs:",len(sel))
os.makedirs("concat",exist_ok=True)
def readfa(t):
    d={};k=None
    for l in t.splitlines():
        if l.startswith(">"): k=l[1:].split()[0]; d[k]=[]
        elif k: d[k].append(l.strip())
    return {a:"".join(b) for a,b in d.items()}
aln={};info=[]
for og in sel:
    seqs=readfa(open(R+f"Orthogroup_Sequences/{og}.fa").read())
    seqs={g.rsplit("_",1)[0]:s for g,s in seqs.items()}
    tmp=f"concat/{og}.fa"
    open(tmp,"w").write("".join(f">{g}\n{s}\n" for g,s in seqs.items()))
    out=subprocess.run(["mafft","--auto","--thread","4",tmp],capture_output=True,text=True,env=env).stdout
    a={g:s.upper() for g,s in readfa(out).items()}
    L=len(next(iter(a.values())))
    keep=[i for i in range(L) if sum(1 for s in a.values() if s[i]=="-")<=0.5*len(a)]
    if len(keep)<80: continue
    aln[og]={g:"".join(s[i] for i in keep) for g,s in a.items()}
    info.append((og,len(a),len(keep)))
with open("concat95.aln.faa","w") as f:
    for g in labs:
        f.write(">"+g+"\n"+"".join(aln[o].get(g,"-"*l) for o,n,l in info)+"\n")
with open("concat95_partitions.tsv","w") as f:
    p=1
    for o,n,l in info:
        f.write(f"{o}\t{n}\t{l}\t{p}-{p+l-1}\n"); p+=l
print("OGs used:",len(info),"| total sites:",sum(l for o,n,l in info))
