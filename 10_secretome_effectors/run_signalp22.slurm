#!/bin/bash
# ANALYSIS: Secretome (SignalP, DIAMOND) and effector search in the 22 genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Fig. S6.2, Tables S8.x, S9.x, S10.x
# FOLDER: 10_secretome_effectors (see README.md in the parent folder)
#SBATCH --job-name=signalp22
#SBATCH --output=/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/06_secretome_effector_22genomes/logs/signalp22_%j.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=batch
W=/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/06_secretome_effector_22genomes
cd /scratch/urbello/envs/funannotate/bin
for F in $W/proteomes22/*.faa; do
  L=$(basename $F .faa); [ -s "$F" ] || { echo "EMPTY $L"; continue; }
  echo "=== $L ==="
  mkdir -p $W/signalp_out/tmp_$L; ./signalp -fasta $(readlink -f $F) -org gram+ -format short -prefix $W/signalp_out/${L}_signalp -tmp $W/signalp_out/tmp_$L
done
echo DONE
