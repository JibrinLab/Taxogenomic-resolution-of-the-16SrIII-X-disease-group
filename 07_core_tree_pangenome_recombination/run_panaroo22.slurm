#!/bin/bash
# ANALYSIS: Core genome tree, pangenome and recombination
# MANUSCRIPT SECTION: Phylogenetic reconstructions; Comparative Genomics Analyses
# USED FOR: Fig. 5, Fig. 6, Figs. S5.1 to S5.6
# FOLDER: 07_core_tree_pangenome_recombination (see README.md in the parent folder)
#SBATCH -J pan22
#SBATCH -p batch
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o pan22_%j.out
source /opt/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate panaroo_env
cd /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/09_core_tree/pps_pruni22
panaroo --version; mafft --version
panaroo -i inputs/*.gff -o panaroo_out --clean-mode sensitive -c 0.7 --len_dif_percent 0.7 --remove-invalid-genes -a core --aligner mafft --core_threshold 0.90 -t 16
module load iqtree/2.4.0
iqtree2 -s panaroo_out/core_gene_alignment.aln -m TEST -B 1000 -T 16 --prefix panaroo_out/core_tree_panaroo_new --redo
echo done
