#!/bin/bash
# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
#SBATCH -J orthofinder84
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 24:00:00
#SBATCH -o orthofinder84_%j.out
export PATH=/projects/dsn011/urbello_conda/envs/orthofinder_env/bin:$PATH
cd /projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/06_core_genome_tree
orthofinder -f proteins -t 16 -a 8 -S diamond -M msa
echo done
