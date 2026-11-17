#!/bin/bash
# ANALYSIS: Average nucleotide identity (FastANI, 83 curated and 299 genomes; linkage sensitivity)
# MANUSCRIPT SECTION: Average nucleotide identity
# USED FOR: ANI levels n=83 and n=299; linkage sensitivity
# FOLDER: 02_ANI (see README.md in the parent folder)
#SBATCH -J fastani_aclai
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 01:00:00
#SBATCH -o fastani_aclai_%j.out
export PATH=/projects/dsn011/urbello_conda/envs/fastani_env/bin:$PATH
D=/projects/dsn011/ppss_publication/2026-10-20_new_species_taxonomy
echo $D/00_inputs/aclai_checkm/GCF_000018785.1_ASM1878v1_genomic.fna > q.txt
fastANI --ql q.txt --rl all299_paths.txt -t 8 --fragLen 3000 --minFraction 0.2 -o fastani_aclai_vs_299.tsv
echo "hits: $(wc -l < fastani_aclai_vs_299.tsv)"
