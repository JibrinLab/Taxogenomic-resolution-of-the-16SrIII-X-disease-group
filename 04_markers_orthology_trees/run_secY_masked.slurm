#!/bin/bash
# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
#SBATCH -J secYmask
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -t 04:00:00
#SBATCH -o secYmask_%j.out
module load iqtree/2.4.0
cd /projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/14_pbt_genbank_refs
/opt/anaconda3/2021.05/bin/python - <<'PY'
d={};k=None
for l in open("secY_plusPBT.aln.faa"):
    l=l.rstrip()
    if l.startswith(">"): k=l[1:].split()[0]; d[k]=[]
    elif k: d[k].append(l)
d={a:"".join(b) for a,b in d.items()}
k="PBT_GU004361_secY"; s=list(d[k])
idx=[i for i,c in enumerate(s) if c!="-"][-10:]
for i in idx: s[i]="X"
d[k]="".join(s)
open("secY_plusPBT_masked.aln.faa","w").write("".join(f">{a}\n{b}\n" for a,b in d.items()))
print("masked",len(idx),"residues")
PY
IQ=$(which iqtree2 2>/dev/null || which iqtree)
$IQ -s secY_plusPBT_masked.aln.faa -m MFP -mset LG,WAG,JTT -B 1000 -T 4 --prefix tree_secY_plusPBT_masked -quiet
echo done
