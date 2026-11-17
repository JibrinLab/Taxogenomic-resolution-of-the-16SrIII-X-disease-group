#!/bin/bash
# ANALYSIS: Whole-genome synteny of the PPS strains against four reference genomes
# MANUSCRIPT SECTION: Comparative Genomics Analyses
# USED FOR: Figs. S1.4 to S1.7
# FOLDER: 11_synteny (see README.md in the parent folder)
#SBATCH -J synt4
#SBATCH -p batch
#SBATCH -N 1 -n 16
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -o /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/07_pangenome_synteny_SNP/synteny_4reps/synt4_%j.out
module load mummer/4.0.0rc1
export G=/projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/06_secretome_effector_22genomes/genomes22
cd /projects/dsn011/ppss_publication/2026-10-27_22genome_reanalysis/07_pangenome_synteny_SNP/synteny_4reps
nucmer --version 2>&1 | head -n 2; dnadiff --version 2>&1 | head -n 2
for r in CanS4 CX CicWB-2022 LWB; do for f in $G/*.fna; do g=$(basename $f .fna); [ $g = $r ] || echo "$r $g"; done; done > pairs.txt
wc -l pairs.txt
xargs -P 16 -L 1 bash -c 'dnadiff -p dnadiff/$0_vs_$1 $G/$0.fna $G/$1.fna > dnadiff/$0_vs_$1.log 2>&1' < pairs.txt
for f in $G/*.fna; do g=$(basename $f .fna); awk -v g=$g '/^>/{if(n!="")print g"\t"n"\t"l; n=substr($1,2); l=0; next}{l+=length($0)}END{print g"\t"n"\t"l}' $f; done > genome_contig_lengths.tsv
echo done; ls dnadiff/*.report | wc -l
