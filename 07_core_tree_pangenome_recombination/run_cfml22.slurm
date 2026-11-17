#!/bin/bash
# ANALYSIS: Core genome tree, pangenome and recombination
# MANUSCRIPT SECTION: Phylogenetic reconstructions; Comparative Genomics Analyses
# USED FOR: Fig. 5, Fig. 6, Figs. S5.1 to S5.6
# FOLDER: 07_core_tree_pangenome_recombination (see README.md in the parent folder)
#SBATCH -J cfml22
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH -o cfml22_%j.out
set -x
W=/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/09_core_tree/pps_pruni22
ALN=$W/panaroo_out/core_gene_alignment_filtered.aln
mkdir -p $W/iqtree_filtered $W/clonalframe
module load iqtree/2.4.0
iqtree2 -s $ALN -m MFP -B 1000 -T 8 --prefix $W/iqtree_filtered/iq22
source /opt/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate clonalframeml_env
ClonalFrameML $W/iqtree_filtered/iq22.treefile $ALN $W/clonalframe/cfml_iq22 -em true -show_progress true
echo CFML_EXIT=$?
