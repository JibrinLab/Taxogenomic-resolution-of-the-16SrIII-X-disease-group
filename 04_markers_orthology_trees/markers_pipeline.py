# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
import glob, os, subprocess
R="proteins/OrthoFinder/Results_Oct20/Orthogroup_Sequences/"
env=dict(os.environ,PATH="/projects/dsn011/urbello_conda/envs/mafft_env/bin:"+os.environ["PATH"])
SINGLE={"rpoB":"OG0000200","secA":"OG0000251","secY":"OG0000279","rps3":"OG0000193","rpl22":"OG0000115"}
RP={"rps10":"OG0000269","rpl3":"OG0000243","rpl4":"OG0000242","rpl2":"OG0000280","rps19":"OG0000116","rpl22":"OG0000115","rps3":"OG0000193","rpl16":"OG0000192","rpl29":"OG0000114","rps17":"OG0000191","rpl14":"OG0000190","rpl24":"OG0000240","rpl5":"OG0000189","rps14":"OG0000113","rps8":"OG0000188","rpl6":"OG0000187","rpl18":"OG0000186","rps5":"OG0000185","rpl30":"OG0000184","rpl15":"OG0000239"}
def parse(t):
    d={};k=None
    for l in t.splitlines():
        if l.startswith(">"): k=l[1:].split()[0]; d[k]=[]
        elif k: d[k].append(l.strip())
    return {a:"".join(b) for a,b in d.items()}
readfa=lambda p:parse(open(p).read())
NT=readfa("cds_nt.fasta")
labs=sorted(os.path.basename(f)[:-4] for f in glob.glob("inputs/*.gff"))
os.makedirs("markers",exist_ok=True)
def prep(name,og):
    best={}
    for g,s in readfa(R+og+".fa").items():
        lb=g.rsplit("_",1)[0]
        if lb not in best or len(s)>len(best[lb][1]): best[lb]=(g,s)
    open(f"markers/{name}.faa","w").write("".join(f">{lb}\n{s}\n" for lb,(g,s) in best.items()))
    open(f"markers/{name}.nt.fna","w").write("".join(f">{lb}\n{NT[g]}\n" for lb,(g,s) in best.items() if g in NT))
    o=subprocess.run(["mafft","--auto","--thread","4",f"markers/{name}.faa"],capture_output=True,text=True,env=env).stdout
    a={k:v.upper() for k,v in parse(o).items()}
    L=len(next(iter(a.values())))
    keep=[i for i in range(L) if sum(1 for s in a.values() if s[i]=="-")<=0.5*len(a)]
    a={k:"".join(s[i] for i in keep) for k,s in a.items()}
    with open(f"markers/{name}.aln.faa","w") as f:
        for k,s in a.items(): f.write(f">{k}\n{s}\n")
    print(f"{name:7s} genomes {len(a):3d} | sites {len(keep):4d}")
    return a
A={}
for n,og in {**SINGLE,**RP}.items(): A[n]=prep(n,og)
def concat(names,out):
    with open(out,"w") as f, open(out+".parts.tsv","w") as p:
        pos=1
        for n in names: l=len(next(iter(A[n].values()))); p.write(f"{n}\t{pos}-{pos+l-1}\n"); pos+=l
        for lb in labs:
            f.write(">"+lb+"\n"+"".join(A[n].get(lb,"-"*len(next(iter(A[n].values())))) for n in names)+"\n")
    print(out,"sites",pos-1)
concat(list(RP),"markers/rp_operon.concat.faa")
concat(["rpoB","secA","secY"]+list(RP),"markers/markers_all.concat.faa")
