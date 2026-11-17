#!/bin/bash
# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
#SBATCH -J pbtmark
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 08:00:00
#SBATCH -o pbtmark_%j.out
module load hmmer
module load iqtree/2.4.0
cd /projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/14_pbt_genbank_refs
/opt/anaconda3/2021.05/bin/python add_pbt_markers.py
IQ=$(which iqtree2 2>/dev/null || which iqtree)
for m in secY rps3 rpl22; do
  [ -f ${m}_plusPBT.aln.faa ] && $IQ -s ${m}_plusPBT.aln.faa -m MFP -mset LG,WAG,JTT -B 1000 -T 2 --prefix tree_${m}_plusPBT -quiet &
done
wait
echo done
