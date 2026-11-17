# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
import glob, re, os
R="/projects/dsn011/ppss_publication/"
T=R+"2026-10-20_new_species_taxonomy/"
ids=[l.strip() for l in open(T+"03_ani_curated/curated_genome_ids.txt") if l.strip()]
AC=T+"00_inputs/aclai_full/ncbi_dataset/data/GCF_000018785.1/"
os.makedirs("inputs",exist_ok=True)
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
def label(g):
    if g.startswith("GC"): return re.sub(r"[^A-Za-z0-9]+","_",g.split("_",2)[-1])
    return g
def fasta_ids(t): return set(l[1:].split()[0] for l in t.splitlines() if l.startswith(">"))
jobs=[(g,gff_for(g),R+f"ani_all299/genomes/{g}.fna") for g in ids]
jobs.append(("A_laidlawii",AC+"genomic.gff",AC+"GCF_000018785.1_ASM1878v1_genomic.fna"))
seen={};m=open("genome_labels.tsv","w");m.write("genome\tlabel\n");bad=[]
for g,f,fa in jobs:
    lab=label(g)
    if lab in seen: lab=lab+"_"+re.search(r"\d{9}",g).group(0)
    seen[lab]=g
    if not f or not os.path.exists(fa): bad.append((g,"missing file"));continue
    seq=open(fa).read(); gf=open(f).read().split("##FASTA")[0]
    ctg=set(l.split("\t")[0] for l in gf.splitlines() if l and not l.startswith("#") and l.count("\t")>=8)
    if not ctg<=fasta_ids(seq): bad.append((g,"contig names differ"));continue
    open(f"inputs/{lab}.gff","w").write(gf.rstrip("\n")+"\n##FASTA\n"+seq)
    m.write(f"{g}\t{lab}\n")
m.close()
print("written:",len(jobs)-len(bad),"of",len(jobs));print("problems:",bad)
