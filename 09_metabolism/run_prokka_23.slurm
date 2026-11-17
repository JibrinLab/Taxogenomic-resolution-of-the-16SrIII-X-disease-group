#!/bin/bash
# ANALYSIS: Metabolic gene panel (66 genes, 7 pathways) in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. 8, Tables S7.1 to S7.10
# FOLDER: 09_metabolism (see README.md in the parent folder)
#SBATCH -J reannot
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 04:00:00
#SBATCH --array=1-23%8
#SBATCH -o /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/08_reannotation_prokka_23genomes/logs/reannot_%A_%a.out
W=/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/08_reannotation_prokka_23genomes
export PATH=/opt/anaconda3/2020.07/bin:$PATH
read ID FA < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $W/genomes_for_annotation.tsv)
LT=${ID//[^A-Za-z0-9]/}
echo "task ${SLURM_ARRAY_TASK_ID}: $ID  tag $LT  $FA"
prokka --version 2>&1 | head -1
prokka --kingdom Bacteria --genus Phytoplasma --locustag $LT --prefix $ID --outdir $W/prokka/$ID --cpus ${SLURM_CPUS_PER_TASK} --rfam --force $FA || echo "prokka exit code $?"
test -s $W/prokka/$ID/$ID.gff && test -s $W/prokka/$ID/$ID.faa && echo "OUTPUT_OK $ID" || echo "OUTPUT_MISSING $ID"
