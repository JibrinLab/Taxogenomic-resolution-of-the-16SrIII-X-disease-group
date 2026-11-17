#!/bin/bash
# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
#SBATCH -J markers
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 12:00:00
#SBATCH -o markers_%j.out
module load iqtree/2.4.0
cd /projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/06_core_genome_tree
PY=/opt/anaconda3/2021.05/bin/python
$PY make_cds_nt.py && $PY markers_pipeline.py
IQ=$(which iqtree2 2>/dev/null || which iqtree)
for a in rpoB secA secY rps3 rpl22 rp_operon.concat markers_all.concat; do
  f=markers/$a.aln.faa; [ -f $f ] || f=markers/$a.faa
  $IQ -s $f -m MFP -mset LG,WAG,JTT -B 1000 -T 2 --prefix markers/tree_$a -quiet &
done
wait
echo done
