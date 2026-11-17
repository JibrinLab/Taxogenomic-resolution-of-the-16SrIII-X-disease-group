#!/bin/bash
# ANALYSIS: OrthoFinder orthogroups of the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Figs. S6.3 to S6.5
# FOLDER: 12_orthogroups_22genomes (see README.md in the parent folder)
#SBATCH -J of22
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/14_orthofinder_22genomes/of22_%j.out
export PATH=/projects/dsn011/urbello_conda/envs/orthofinder_env/bin:$PATH
cd /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/14_orthofinder_22genomes
orthofinder -h 2>&1 | head -n 3
orthofinder -f proteins22 -t 16 -a 8 -S diamond -M msa -n r22
orthofinder -f proteins20 -t 16 -a 8 -S diamond -M msa -n r4pps
echo done
