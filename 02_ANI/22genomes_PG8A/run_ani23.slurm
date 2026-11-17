#!/bin/bash
# ANALYSIS: Average nucleotide identity of the 22 genomes plus PG-8A (FastANI, ANIm, ANIb)
# MANUSCRIPT SECTION: Average nucleotide identity
# USED FOR: Fig. 2, Fig. S3.1, Tables S3.1 to S3.5
# FOLDER: 02_ANI/22genomes_PG8A (see README.md in the parent folder)
#SBATCH -J ani23
#SBATCH -p batch
#SBATCH -c 16
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -o /projects/dsn011/ppss_publication/ani_23_with_PG8A/ani23_%j.out
cd /projects/dsn011/ppss_publication/ani_23_with_PG8A
module load blast+/2.16.0
export PATH=/projects/dsn011/urbello_conda/pyani_venv/bin:/opt/anaconda3/2021.05/bin:$PATH
which nucmer delta-filter blastn makeblastdb
nucmer --version
echo ANIm start $(date)
average_nucleotide_identity.py -i genomes -o pyani_anim -m ANIm --scheduler multiprocessing --workers 16 --force -v
echo ANIm end $(date)
echo ANIb start $(date)
average_nucleotide_identity.py -i genomes -o pyani_anib -m ANIb --scheduler multiprocessing --workers 16 --force -v
echo ANIb end $(date)
ls -d $PWD/genomes/*.fna > genome_list.txt
echo "genomes: $(wc -l < genome_list.txt)"
source /opt/anaconda3/2021.05/etc/profile.d/conda.sh
conda activate /projects/dsn011/urbello_conda/envs/fastani_env
fastANI --ql genome_list.txt --rl genome_list.txt -t ${SLURM_CPUS_PER_TASK} --fragLen 3000 --minFraction 0.2 -o fastani_23.tsv
echo done $(date)
