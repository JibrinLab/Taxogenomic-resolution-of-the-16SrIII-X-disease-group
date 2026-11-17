#!/bin/bash
# ANALYSIS: GTDB-Tk classification of 53 genomes
# MANUSCRIPT SECTION: ANI, dDDH and GTDB clustering
# USED FOR: Tables S5.1 to S5.3
# FOLDER: 06_GTDB-Tk (see README.md in the parent folder)
#SBATCH -J gtdbtk53
#SBATCH -p bigmem
#SBATCH -N 1 -n 32
#SBATCH --mem=700G
#SBATCH -t 24:00:00
#SBATCH -o gtdbtk53_%j.out
module load gtdbtk/2.6.1
module load gtdb/r226
cd /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/10_gtdbtk
gtdbtk --version
module list 2>&1
echo "GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH"
gtdbtk classify_wf --genome_dir genomes --out_dir gtdbtk_out --extension fna --cpus 32 --pplacer_cpus 1
echo "done: $(date)"
