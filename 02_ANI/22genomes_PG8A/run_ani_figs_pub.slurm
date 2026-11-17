#!/bin/bash
# ANALYSIS: Average nucleotide identity of the 22 genomes plus PG-8A (FastANI, ANIm, ANIb)
# MANUSCRIPT SECTION: Average nucleotide identity
# USED FOR: Fig. 2, Fig. S3.1, Tables S3.1 to S3.5
# FOLDER: 02_ANI/22genomes_PG8A (see README.md in the parent folder)
#SBATCH -J anifig
#SBATCH -p express
#SBATCH -N 1 -n 2
#SBATCH --mem=24G
#SBATCH -t 0:30:00
#SBATCH -o /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/01_ANI_23genomes_PG8A/05_logs/anifig_%j.out
cd /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/01_ANI_23genomes_PG8A
/projects/dsn011/urbello_conda/pyani_venv/bin/python 04_scripts/make_ani_figures_pub.py && echo done
ls -la 01_figures/publication_style | grep -E "ANIm|ANIb|FastANI"
