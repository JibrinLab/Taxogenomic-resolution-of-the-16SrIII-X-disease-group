#!/bin/bash
# ANALYSIS: Trees with the four species marked
# MANUSCRIPT SECTION: Figures
# USED FOR: Tree figures
# FOLDER: 13_trees_taxa (see README.md in the parent folder)
#SBATCH --job-name=taxa_trees_v3
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=draw_taxa_trees_v3_%j.out
cd /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis
R=/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy
M=$R/06_core_genome_tree/markers
rm -rf taxa_trees_v3 taxa_core_four_groups
# all trees, every taxon labelled (branches flipped so each taxon sits together)
python3 tree_taxa.py --map taxon_map.tsv taxa_trees_v3 \
  $R/06_core_genome_tree/core84_tree.contree \
  $R/06_core_genome_tree/concat95_tree.contree \
  $R/04_16S/16S_rep_tree.contree \
  $M/tree_rp_operon.concat.contree \
  $M/tree_secY.contree \
  $M/tree_rpoB.contree \
  $M/tree_markers_all.concat.contree \
  12_marker_gene_trees/tuf_tree.contree
# core gene tree with only the four groups (22 genomes plus outgroup)
python3 tree_taxa.py --groups-only taxa_core_four_groups $R/06_core_genome_tree/core84_tree.contree
