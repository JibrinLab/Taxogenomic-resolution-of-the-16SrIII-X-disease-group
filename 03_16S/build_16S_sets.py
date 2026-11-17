# ANALYSIS: 16S rRNA gene analysis
# MANUSCRIPT SECTION: Multilocus sequence comparison and 16S classification
# USED FOR: Figs. S2.1 to S2.7, Tables S2.1 to S2.6
# FOLDER: 03_16S (see README.md in the parent folder)
import re
T="/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/"
AC=T+"00_inputs/aclai_full/ncbi_dataset/data/GCF_000018785.1/"
lab=dict(l.rstrip("\n").split("\t") for l in open(T+"06_core_genome_tree/genome_labels.tsv"))
comp=str.maketrans("ACGTacgtNn","TGCAtgcaNn")
recs=[]
k=None
for l in open("16S_all_operons.fasta"):
    if l.startswith(">"): k=l[1:].strip().split("|"); recs.append([k[0],k[1],int(k[2]),[]])
    else: recs[-1][3].append(l.strip())
recs=[(g,c,n,"".join(s)) for g,c,n,s in recs]
d={};k=None
for l in open(AC+"GCF_000018785.1_ASM1878v1_genomic.fna"):
    if l.startswith(">"): k=l[1:].split()[0]; d[k]=[]
    else: d[k].append(l.strip())
d={a:"".join(b) for a,b in d.items()}
for l in open(AC+"genomic.gff"):
    c=l.rstrip("\n").split("\t")
    if len(c)>8 and c[2]=="rRNA" and re.search("16S",c[8],re.I):
        s=d[c[0]][int(c[3])-1:int(c[4])]
        if c[6]=="-": s=s.translate(comp)[::-1]
        recs.append(("A_laidlawii",f"{c[0]}:{c[3]}-{c[4]}",len(s),s))
by={}
for g,c,n,s in recs:
    if n>=1200: by.setdefault(g,[]).append((n,c,s))
fa=open("16S_allcopies.fasta","w"); fr=open("16S_rep.fasta","w"); tb=open("16S_set_table.tsv","w")
tb.write("genome\tlabel\tcopy\tlength\tcontig_coords\n")
for g,L in by.items():
    lb=lab.get(g,g)
    L.sort(reverse=True)
    for i,(n,c,s) in enumerate(L,1):
        fa.write(f">{lb}_c{i}\n{s}\n"); tb.write(f"{g}\t{lb}\t{i}\t{n}\t{c}\n")
    fr.write(f">{lb}\n{L[0][2]}\n")
print("genomes with >=1200 bp 16S (incl. A_laidlawii):",len(by))
print("total copies:",sum(len(v) for v in by.values()))
print("A_laidlawii copies >=1200:",len(by.get("A_laidlawii",[])))
