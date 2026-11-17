#!/bin/bash
# ANALYSIS: CheckM genome quality
# MANUSCRIPT SECTION: Additional genomes used
# USED FOR: Table S1.6
# FOLDER: 01_genome_quality (see README.md in the parent folder)
#SBATCH -J checkm299
#SBATCH -p batch
#SBATCH -t 48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -o checkm299_%j.out
export PATH=/projects/dsn011/urbello_conda/envs/checkm_env/bin:$PATH
IN=/projects/dsn011/ppss_publication/ani_all299/genomes
OUT=/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy/01_genome_quality/checkm_all299
echo "start $(date)"
checkm lineage_wf -t 32 -x fna --tab_table -f checkm_lineage_wf_summary.tsv $IN $OUT
checkm qa -o 2 --tab_table -f checkm_all299_extended.tsv -t 32 $OUT/lineage.ms $OUT
echo "end $(date) exit=$?"
