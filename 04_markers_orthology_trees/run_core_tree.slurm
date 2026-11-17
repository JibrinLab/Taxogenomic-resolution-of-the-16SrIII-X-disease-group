#!/bin/bash
# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
#SBATCH -J core_tree
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 24:00:00
#SBATCH -o core_tree_%j.out
module load iqtree/2.4.0
cd /projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/06_core_genome_tree
IQ=$(which iqtree2 2>/dev/null || which iqtree)
$IQ -s panaroo_out/core_gene_alignment.aln -m MFP -B 1000 -T 16 --prefix core84_tree -o A_laidlawii
echo done
