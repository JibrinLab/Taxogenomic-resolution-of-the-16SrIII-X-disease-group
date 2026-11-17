#!/bin/bash
# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
#SBATCH -J pbt16S
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 04:00:00
#SBATCH -o pbt16S_%j.out
export PATH=/projects/dsn011/urbello_conda/envs/mafft_env/bin:$PATH
module load iqtree/2.4.0
cd /projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/14_pbt_genbank_refs
mafft --add GU004371.fa --keeplength --thread 8 ../04_16S/16S_rep.aln.fasta > 16S_rep_plusPBT.aln.fasta
IQ=$(which iqtree2 2>/dev/null || which iqtree)
$IQ -s 16S_rep_plusPBT.aln.fasta -m MFP -B 1000 -T 8 --prefix 16S_rep_plusPBT_tree -o A_laidlawii
/opt/anaconda3/2021.05/bin/python - <<'PY'
d={};k=None
for l in open("16S_rep_plusPBT.aln.fasta"):
    l=l.strip()
    if l.startswith(">"): k=l[1:].split()[0]; d[k]=""
    elif k: d[k]+=l.upper()
q=d["PBT_GU004371"];r=[]
for g,s in d.items():
    if g==q: continue
    p=[(a,b) for a,b in zip(q,s) if a in "ACGT" and b in "ACGT"]
    if p and g!="PBT_GU004371": r.append((sum(a==b for a,b in p)/len(p),len(p),g))
r.sort(reverse=True)
for i,n,g in r[:25]: print(f"{g:14s} {100*i:7.3f}% over {n} bp")
PY
echo done
