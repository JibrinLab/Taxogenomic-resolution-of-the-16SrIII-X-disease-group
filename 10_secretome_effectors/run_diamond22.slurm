#!/bin/bash
# ANALYSIS: Secretome (SignalP, DIAMOND) and effector search in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. S6.2, Tables S8.x, S9.x, S10.x
# FOLDER: 10_secretome_effectors (see README.md in the parent folder)
#SBATCH -J dmnd22
#SBATCH -p express
#SBATCH -N 1 -n 8
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -o /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/06_secretome_effector_22genomes/diamond/dmnd22_%j.out
module load diamond/2.1.13
cd /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/06_secretome_effector_22genomes
diamond version
diamond makedb --in blast/all_proteomes.faa -d diamond/all22 --threads 8
diamond blastp -q diamond/candidates_SPSec.faa -d diamond/all22 --evalue 1e-5 --id 30 -k 0 --threads 8 --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp evalue bitscore -o diamond/hits.tsv
echo done; wc -l diamond/hits.tsv
