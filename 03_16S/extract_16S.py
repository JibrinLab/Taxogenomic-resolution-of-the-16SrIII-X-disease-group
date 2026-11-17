# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
import glob, re, os, sys
R="/projects/dsn011/ppss_publication/"
T=R+"2026-10-20_new_species_taxonomy/"
ids=[l.strip() for l in open(T+"03_ani_curated/curated_genome_ids.txt") if l.strip()]
comp=str.maketrans("ACGTacgtNn","TGCAtgcaNn")
def fasta(p):
    d={};k=None
    for l in open(p):
        if l.startswith(">"): k=l[1:].split()[0]; d[k]=[]
        else: d[k].append(l.strip())
    return {a:"".join(b) for a,b in d.items()}
def gff_for(g):
    n=re.search(r"\d{9}",g)
    if g.startswith("PPS"):
        c=glob.glob(R+f"gff_clean/{g}.gff")+glob.glob(R+f"pps/22_prokka_final/{g}_*/*.gff")
    elif "OVSP" in g:
        s="OVSPS01" if "S01" in g else "OVSPL01"
        c=glob.glob(T+f"10_annotation_missing/prokka_t11/{s}/{s}.gff")
    else:
        c=glob.glob(R+f"annotations_gff/{g}.gff")+glob.glob(R+f"annotations_gff/*{n.group(0)}*.gff")
    return c[0] if c else None
out=open("16S_all_operons.fasta","w"); rows=[]
for g in ids:
    f=gff_for(g); gen=R+f"ani_all299/genomes/{g}.fna"
    if not f or not os.path.exists(gen): rows.append((g,"NO_GFF" if not f else "NO_FNA",0,"")); continue
    seq=fasta(gen); L=[]
    for l in open(f):
        if l.startswith("#") or l.count("\t")<8: continue
        c=l.rstrip("\n").split("\t")
        if c[2]!="rRNA" or not re.search(r"16S",c[8],re.I): continue
        if c[0] not in seq: continue
        s,e=int(c[3]),int(c[4]); x=seq[c[0]][s-1:e]
        if c[6]=="-": x=x.translate(comp)[::-1]
        L.append(len(x)); out.write(f">{g}|{c[0]}:{s}-{e}({c[6]})|{len(x)}\n{x}\n")
    rows.append((g,os.path.basename(f),len(L),",".join(map(str,L))))
out.close()
with open("16S_summary.tsv","w") as o:
    o.write("genome\tgff\tn_16S\tlengths\n")
    for r in rows: o.write("\t".join(map(str,r))+"\n")
n0=sum(1 for r in rows if r[2]==0); n12=sum(1 for r in rows if any(int(x)>=1200 for x in r[3].split(",") if x))
print("genomes:",len(rows),"| with no 16S found:",n0,"| with a copy >=1200 bp:",n12)
for r in rows:
    if r[2]==0: print("  none:",r[0],r[1])
