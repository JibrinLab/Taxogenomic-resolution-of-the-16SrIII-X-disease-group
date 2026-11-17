#!/bin/bash
# ANALYSIS: Core genome tree, pangenome and recombination
# MANUSCRIPT SECTION: Phylogenetic reconstructions; Comparative Genomics Analyses
# USED FOR: Fig. 5, Fig. 6, Figs. S5.1 to S5.6
# FOLDER: 07_core_tree_pangenome_recombination (see README.md in the parent folder)
#SBATCH -J grpPR
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 06:00:00
#SBATCH --array=1-3
#SBATCH -o grp%a_%A.out
# Per-group Panaroo (relaxed) + IQ-TREE + ClonalFrameML. Submit from the folder that holds this file.
set -x
P=/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis
W=$P/15_group_pangenome_recombination
GFF=$P/05_metabolism_22genomes/gff_tree
case $SLURM_ARRAY_TASK_ID in
 1) IDS="PPS1A PPS2B PPS6A PPS6B PPS9B PPS10A CanS4 MW1 VAC MA1";;
 2) IDS="CX 6A1 ZT3-1 PR2021 JR1 PDA15";;
 3) IDS="CicWB-2022 ChTDIII OVSP-L01 OVSP-S01 Vc33";;
esac
G=$W/Group$SLURM_ARRAY_TASK_ID
mkdir -p $G/inputs; cd $G
for i in $IDS; do ln -sf $GFF/$i/$i.gff inputs/$i.gff; done
ls inputs | wc -l

source /opt/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate panaroo_env
panaroo --version; mafft --version
panaroo -i inputs/*.gff -o panaroo_out --clean-mode sensitive -c 0.7 --len_dif_percent 0.7 \
  --remove-invalid-genes -a core --aligner mafft --core_threshold 0.90 -t 8
echo PANAROO_EXIT=$?

ALN=panaroo_out/core_gene_alignment_filtered.aln
[ -s $ALN ] || ALN=panaroo_out/core_gene_alignment.aln
mkdir -p clonalframe iqtree
# ClonalFrameML accepts A, C, G, T, N and gaps only
awk '/^>/{print;next}{gsub(/[^ACGTacgt-]/,"N");print}' $ALN > clonalframe/core.ACGTN.aln

module load iqtree/2.4.0
iqtree2 -s clonalframe/core.ACGTN.aln -m MFP -B 1000 -T 8 --prefix iqtree/tree --redo
echo IQTREE_EXIT=$?

conda activate clonalframeml_env
ClonalFrameML iqtree/tree.treefile clonalframe/core.ACGTN.aln clonalframe/cfml -em true -show_progress true
echo CFML_EXIT=$?
