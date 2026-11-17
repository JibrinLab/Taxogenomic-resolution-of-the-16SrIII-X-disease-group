#!/bin/bash
# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
#SBATCH -J panaroo84
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -t 24:00:00
#SBATCH -o panaroo84_%j.out
export PATH=/projects/dsn011/urbello_conda/envs/panaroo_env/bin:/projects/dsn011/urbello_conda/envs/mafft_env/bin:$PATH
cd /projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/06_core_genome_tree
panaroo -i inputs/*.gff -o panaroo_out --clean-mode sensitive -c 0.7 --len_dif_percent 0.7 \
        --remove-invalid-genes -a core --aligner mafft --core_threshold 0.90 -t 16
echo done
