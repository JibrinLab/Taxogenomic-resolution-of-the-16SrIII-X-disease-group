#!/bin/bash
# ANALYSIS: Metabolic gene panel (66 genes, 7 pathways) in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 8, Tables S7.1 to S7.10
# FOLDER: 09_metabolism (see README.md in the parent folder)
#SBATCH -J prokka7
#SBATCH -p batch
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH --array=1-7
#SBATCH -o prokka7_%A_%a.out
module load anaconda3/2020.07
export PATH=/opt/anaconda3/2020.07/bin:$PATH
W=/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/05_metabolism_22genomes
read ID FA < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $W/genomes_new7.tsv)
LT=${ID//[^A-Za-z0-9]/}
prokka --version 2>&1
prokka --kingdom Bacteria --genus Phytoplasma --locustag $LT --prefix $ID --outdir $W/prokka_new7/$ID --cpus 4 --rfam --force --compliant $FA
ls $W/prokka_new7/$ID/$ID.gff && echo OUTPUT_OK
