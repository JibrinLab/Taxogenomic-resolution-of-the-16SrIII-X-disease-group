#!/bin/bash
# ANALYSIS: CheckM genome quality
# MANUSCRIPT SECTION: Additional genomes used
# USED FOR: Table S1.6
# FOLDER: 01_genome_quality (see README.md in the parent folder)
#SBATCH -J checkm_aclai
#SBATCH -p batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH -t 04:00:00
#SBATCH -o checkm_aclai_%j.out
export PATH=/projects/dsn011/urbello_conda/envs/checkm_env/bin:$PATH
IN=/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/00_inputs/aclai_checkm
OUT=/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/01_genome_quality/checkm_aclai
checkm lineage_wf -t 16 -x fna --tab_table -f $OUT.summary.tsv $IN $OUT
