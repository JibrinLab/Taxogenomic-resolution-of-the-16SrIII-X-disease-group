#!/bin/bash
# ANALYSIS: Marker genes, orthology and phylogenetic trees
# MANUSCRIPT SECTION: Multilocus sequence comparison; Phylogenetic reconstructions
# USED FOR: Fig. S2.8, Figs. S4.1 to S4.6
# FOLDER: 04_markers_orthology_trees (see README.md in the parent folder)
#SBATCH -J cxrpoB
#SBATCH -p batch
#SBATCH -c 4
#SBATCH -t 02:00:00
#SBATCH -o cxrpoB_%j.out
module load iqtree/2.4.0
export PATH=$PATH:/projects/dsn011/urbello_conda/envs/panaroo_env/bin
cd /projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy
/projects/dsn011/urbello_conda/envs/panaroo_env/bin/python add_cx_rpoB.py
cd 06_core_genome_tree/markers
which mafft
mafft --add CX_rpoB_masked.faa --keeplength rpoB.aln.faa > rpoB_plusCX.aln.faa
grep -c ">" rpoB_plusCX.aln.faa
IQ=$(which iqtree2 2>/dev/null || which iqtree)
$IQ -s rpoB_plusCX.aln.faa -m MFP -mset LG,WAG,JTT -B 1000 -T 4 --prefix tree_rpoB_plusCX -quiet
echo done
